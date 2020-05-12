  ############# ESTIMATE THE LINKING PROBABILITIES ###############
  library(maxLik)
  library(readstata13)
  library(numDeriv)
  dta <- read.dta13("/home/vincent/Desktop/tmpah/intercultural/newdata.dta") # read formatted data
  dta$dat <- NULL
  dta[dta==(-1)] <- NA #recode missing
  
  source("allfunctions.R")
  
  ################
  ### RUN CODE
  ################
  
  
  ## get formatted data
  outbiglist <- buildprobadta()
  outdta <- outbiglist[[1]]
  X <- outbiglist[[2]]
  D <- outbiglist[[3]]
  S <- outbiglist[[4]]
  G <- outbiglist[[5]]
  rm(outbiglist)
  
  outdta$dist <- outdta$dist/max(outdta$dist) # rescale distance 
  outdta <- outdta[outdta$self==0,] # remove self links
  outdta$self <- NULL # remove "self" variable
  N <- nrow(outdta) # number of pairs
  
  ## school dummies matrix for pairwise regression
  schdummy <- matrix(0,N,length(unique(outdta$school)))
  for (i in 1:ncol(schdummy)){
    schdummy[,i] <- as.numeric(outdta$school==i)
  }
  
  ## factor variable for LPM startup
  outdta$school.f <- factor(outdta$school)
  
  ## initial value from LPM
  olsout <- lm(g ~ 0 + white:social + black:social + hisp:social + asian:social + mwork:social + gender:social + age:social + dist:social + typeLH:social + typeHL:social + school.f:social,data = outdta)
  
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
    
rm(D,dta,outdta,S,schdummy,X,G,P,parentsdta,parentsdta0,parentsdta1) # remove confidential data
save.image("outestim_partial.RData") # save estimation results


#### Not used

#  hld <- 99 # initialize bound on lambda
#  for (i in 1:length(D)){
    
#    ## build probability matrix
#    nt <- nrow(X[[i]])
#    Pt <- G[[i]]
#    diag(Pt) <- 0
#    hld <- min(hld,(1/norm(Pt,"2"))) # update bound (tighten)
#  }
#  ls <- optim(0,liksar_net,method="Brent",lower=(-hld+1e-8),upper=(hld-1e-8)) # optimize SAR likelihood G
#  lic <- function(lbda) liksar(lbda,fullmax$par) # creates a function of lambda
#  lsc <- optim(0,lic,method="Brent",lower=(-hld+1e-8),upper=(hld-1e-8)) # optimize SAR likelihood D


