################################################################
################################################################
################### Children's model ###########################
################################################################
################################################################
set.seed(1234)
  library(maxLik)
  library(readstata13)
  library(numDeriv)
  dta <- read.dta13("/home/vincent/Desktop/tmpah/intercultural/newdata.dta") # read formatted data
  dta$dat <- NULL
  dta[dta==(-1)] <- NA #recode missing
  
  source("allfunctions.R")
  
  ################################################################
  ###################### Format data #############################
  ################################################################
  
  ## get formatted data
  outbiglist <- buildprobadta()
  outdta <- outbiglist[[1]]
  X <- outbiglist[[2]]
  D <- outbiglist[[3]]
  S <- outbiglist[[4]]
  G <- outbiglist[[5]]
  Nvec <- outbiglist[[6]]
  sid <- outbiglist[[7]]
  rm(outbiglist)
  md <- max(outdta$dist)
  outdta$dist <- outdta$dist/md # rescale distance
  for (i in 1:length(X)){
    D[[i]][[8]] <- D[[i]][[8]]/md
  }
  outdta <- outdta[outdta$self==0,] # remove self links
  outdta$self <- NULL # remove "self" variable
  
  outbiglist <- droppsmallgroups(10) # drop groups of less than 10 students
  outdta <- outbiglist[[1]]
  X <- outbiglist[[2]]
  D <- outbiglist[[3]]
  S <- outbiglist[[4]]
  G <- outbiglist[[5]]
  Nvec <- outbiglist[[6]]
  sid <- outbiglist[[7]]
  rm(outbiglist)
  
  N <- nrow(outdta) # number of pairs
  ## school dummies matrix for pairwise regression
  schdummy <- matrix(0,N,length(unique(outdta$school)))
  for (i in 1:ncol(schdummy)){
    schdummy[,i] <- as.numeric(outdta$school==i)
  }
  
  
  ################################################################
  ####################### Estimation #############################
  ################################################################
  
  ## factor variable for LPM startup
  outdta$school.f <- factor(outdta$school)
  outdta$socialols <- outdta$social #/outdta$n
  
  ## initial value from LPM
  olsout <- lm(g ~ 0 + white:socialols + black:socialols + hisp:socialols + asian:socialols + mwork:socialols + gender:socialols + age:socialols + dist:socialols +
                 typeLH:socialols + typeHL:socialols + school.f:socialols ,data = outdta)
  
  ## initial value updated using conditional pairwise regression P(G|s)
  theta0 <- as.numeric(olsout$coefficients)
  fstry <- maxLik(logLik=objproba,grad=gradlik,start = theta0)
    
  ## Joint likelihood estimation
  fullmax <- optim(as.numeric(fstry$estimate),jlik)
  
  nproba <- length(theta0)
  nbeta <- length(c(estbeta))
  
  jlik(fullmax$par) ## gets optimum values for lambda and s2 (saved as global variables for the optimal parameter values)
  fulltheta <- c(fullmax$par,estbeta,lambdaest,s2est) # full model's parameters
  
  VV <- hessian(jlikhes,fulltheta) # numerical hessian
  VC <- solve(-VV) # variance-covariance matrix
  jlik(fullmax$par) ## gets optimum values for lambda and s2 (saved as global variables for the optimal parameter values)
  
rm(D,dta,outdta,S,schdummy,X,G,P,parentsdta,parentsdta0,parentsdta1) # remove confidential data
save.image("outestim_partial.RData") # save estimation results




