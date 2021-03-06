################################################################
################################################################
################### Children's model ###########################
################################################################
################################################################
set.seed(1234)
  library(maxLik)
  library(readstata13)
  library(numDeriv)
  setwd("~/Dropbox/intertransmission/revisionJOLE")
  dta <- read.dta13("/home/vincent/Desktop/tmpah/intercultural/newdata.dta") # read formatted data of secure server
  dta$dat <- NULL
  dta[dta==(-1)] <- NA #recode missing
  
  source("allfunctions.R")
  
  ################################################################
  ###################### Format data #############################
  ################################################################
  
  outbiglist <- buildprobadta() # get formatted data
  outdta <- outbiglist[[1]] # pair-level database
  X <- outbiglist[[2]] # individual characteristics (list)
  D <- outbiglist[[3]] # pairs' characteristics (list)
  S <- outbiglist[[4]] # socialization effort (list)
  G <- outbiglist[[5]] # observed network (list)
  Nvec <- outbiglist[[6]] # vector or group sizes
  sid <- outbiglist[[7]] # vector of id school numbers
  rm(outbiglist)
  
  ## rescale distance
  md <- max(outdta$dist)
  outdta$dist <- outdta$dist/md # rescale distance
  for (i in 1:length(X)){
    D[[i]][[8]] <- D[[i]][[8]]/md
  }
  
  #clean the data further
  outdta <- outdta[outdta$self==0,] # remove self links
  outdta$self <- NULL # remove "self" variable
  
  ## remove small groups and update variables
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
  
  ## initial value updated using conditional pairwise regression P(G|s)
  theta0 <- runif(25)
  fstry <- maxLik(logLik=objproba,grad=gradlik,start = theta0)
    
  ## Joint likelihood estimation
  fullmax <- optim(as.numeric(fstry$estimate),jlik)
  
  nproba <- length(theta0)
  nbeta <- length(c(estbeta))
  
  jlik(fullmax$par) ## gets optimum values for lambda and s2 (saved as global variables for the optimal parameter values)
  fulltheta <- c(fullmax$par,estbeta,lambdaest,s2est) # full model's parameters
  
  grplist <- unique(outdta$group)
  
  VV <- hessian(jlikhes,fulltheta) # numerical hessian
  VV <- solve(VV)
  outerproduct <- matrix(NA,length(grplist),length(fulltheta))
  for (grp in 1:length(grplist)){
    lfct <- function(ltheta) jlikhes_level(ltheta,grp)
    outerproduct[grp,] <- grad(lfct,fulltheta)
  }
  Outmat <- matrix(0,length(fulltheta),length(fulltheta))
  for (grp in 1:length(grplist)){
    Outmat <- Outmat + matrix(outerproduct[grp,],length(fulltheta),1)%*%matrix(outerproduct[grp,],1,length(fulltheta))
  }
  Vtheta <- VV%*%Outmat%*%VV
  
  stheta <- sqrt(diag(Vtheta))
  
  jlikhes(fulltheta) ## gets optimum values for lambda and s2 (saved as global variables for the optimal parameter values) and P
  
  prtheta <- c(fulltheta[1:11],fulltheta[26:33],fulltheta[48:49])
  prstheta <- c(stheta[1:11],stheta[26:33],stheta[48:49]) 
  listvars <- c("Both White", "Both Black", "Both Hisp.", "Both Asian", "Both Mother Work", "Same Gender", "Age Difference", "Geographical distance",
                "$d_{ij}(H,L)", "$d_{ij}(L,H)", "$d_{ij}(H,H)",
                "White", "Back", "Hisp.", "Asian", "Mother works", "Female", "Age", "Type H", "$\\phi$", "$\\sigma^2")
  for (i in 1:length(listvars)){
    cat(paste(listvars[i], " & ", format(round(prtheta[i], 3), nsmall=3), " & (", format(round(prstheta[i], 3), nsmall = 3), ") \\\\ \n", sep=""))
  }
  

  
rm(D,dta,outdta,S,schdummy,X,G,P,parentsdta,parentsdta0,parentsdta1) # remove confidential data
save.image("outestim_partial.RData") # save estimation results

