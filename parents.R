################################################################
################################################################
#################### Parents' model ############################
################################################################
################################################################

load("~/Dropbox/intertransmission/newcodes/grade-school/outestim_partial.RData")
library(mvtnorm)
library(maxLik)
library(readstata13)
library(numDeriv)
dta <- read.dta13("/home/vincent/Desktop/tmpah/intercultural/newdata.dta") # read formatted data
dta$dat <- NULL
dta[dta==(-1)] <- NA #recode missing

nsim <- 100 # number of simulations for Eh_i^t
bsim <- 500 # number of bootstrap replications

source("allfunctions.R")


################################################################
###################### Format data #############################
################################################################

  ## get formatted data (children's model)
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

  ### get formatted data (parents' model)

  parentsdta <- builddta() # get formatted data
  parentsdta$school <- factor(parentsdta$scid) # schools as factor var.
  
  ## keep only if present in the children's data
  parentsdta$scid_bkp <- parentsdta$scid
  sl <- unique(parentsdta$scid_bkp)
  rsl <- 1:length(sl)
  for (i in 1:length(sl)){
    parentsdta$scid_bkp[parentsdta$scid==sl[i]] <- rsl[i]
  }
  parentsdta$scid <- parentsdta$scid_bkp*100 + parentsdta$h1gi20
  parentsdta <- parentsdta[sapply(parentsdta$scid, function(x) is.element(x,outdta$group)),]
  
  parentsdta0 <- parentsdta[parentsdta$type==0,] # lowly educated
  parentsdta1 <- parentsdta[parentsdta$type==1,] # higly educated

  ################################################################
  ####################### Estimation #############################
  ################################################################

  result0 <- result1 <- matrix(NA,bsim,23) # initialize results matrices
  for (s in 1:bsim){
    gammatilde <- rmvnorm(1,fulltheta,VC) # draw coefficient
    parentsdta$hom <- buildEhomophily(gammatilde) # get Eh_i^t
    parentsdta0$hom <- parentsdta[parentsdta$type==0,"hom"]
    parentsdta1$hom <- parentsdta[parentsdta$type==1,"hom"]
    parentsdta0$het <- 1-parentsdta0$hom
    parentsdta1$het <- 1-parentsdta1$hom

    obser0 <- sample(1:nrow(parentsdta0),replace=T) # bootstrap observations
    obser1 <- sample(1:nrow(parentsdta1),replace=T) # bootstrap observations
    out0 <- lm(peffort ~ het:white + het:black + het:hisp + het:asian + het:female + het:age + hom + het:school, data=parentsdta0[obser0,]) # regression t=L
    out1 <- lm(peffort ~ hom:white + hom:black + hom:hisp + hom:asian + hom:female + hom:age + hom + hom:school, data=parentsdta1[obser1,]) # regression t=H
    result0[s,1:length(out0$coefficients)] <- out0$coefficients
    result1[s,1:length(out1$coefficients)] <- out1$coefficients
    print(s)
  }

  ## point estimates
  parentsdta$hom <- buildEhomophily(fulltheta) # get Eh_i^t
  out0 <- lm(peffort ~ het:white + het:black + het:hisp + het:asian + het:female + het:age + hom + het:school, data=parentsdta0)
  out1 <- lm(peffort ~ hom:white + hom:black + hom:hisp + hom:asian + hom:female + hom:age + hom + hom:school, data=parentsdta1)

  rm(D,dta,outdta,S,schdummy,X,G,P,parentsdta,parentsdta0,parentsdta1) # remove confidential data
  save.image("outestim_final.RData") # save estimation results
