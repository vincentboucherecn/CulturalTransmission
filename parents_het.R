################################################################
################################################################
#################### Parents' model ############################
################################################################
################################################################

set.seed(1234)
setwd("~/Dropbox/intertransmission/revisionJOLE")
load("outestim_partial.RData")
library(mvtnorm)
library(maxLik)
library(flexmix)
library(readstata13)
library(numDeriv)
dta <- read.dta13("/home/vincent/Desktop/tmpah/intercultural/newdata.dta") # read formatted data from secure server
dta$dat <- NULL
dta[dta==(-1)] <- NA #recode missing

nsim <- 500 # number of simulations for Eh_i^t
bsim <- 500 # number of bootstrap replications

source("allfunctions.R")


################################################################
###################### Format data #############################
################################################################

  ### reload data from children model

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
numberschools <- ncol(schdummy)
schdummy[,numberschools] <- 1 # lastschool is the constant

  ### get formatted data (parents' model)

  parentsdta <- builddta() # get formatted data
  parentsdta$school <- parentsdta$scid # get school id
  parentsdta$schoolf <- factor(parentsdta$scid) # schools as factor var.
  
  ## keep only the data present in the children's data
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

  ## point estimates
  parentsdta$hom <- buildEhomophily(fulltheta) # get Eh_i^t
  parentsdta0$hom <- parentsdta[parentsdta$type==0,"hom"]
  parentsdta1$hom <- parentsdta[parentsdta$type==1,"hom"]
  parentsdta0$het <- 1-parentsdta0$hom
  parentsdta1$het <- 1-parentsdta1$hom
  out0 <- flexmix(peffort ~ 0 + I(hom*white) + I(hom*black) + I(hom*hisp) + I(hom*asian) + I(hom*momworks) + I(hom*female) + I(hom*age) + hom:schoolf + het, data=parentsdta0, k=2)
  out1 <- flexmix(peffort ~ 0 + I(het*white) + I(het*black) + I(het*hisp) + I(het*asian) + I(het*momworks) + I(het*female) + I(het*age) + het:schoolf + hom, data=parentsdta1, k=2)
  PEout0 <- out0 # backup estimate
  PEout1 <- out1 # backup estimate

  ######
  #Bootstrap SE
  ######
   
  result01 <- result02 <- result11 <- result12 <- matrix(NA,bsim,30) # initialize results matrices
  for (s in 1:bsim){
    gammatilde <- rmvnorm(1,fulltheta,Vtheta) # draw coefficient
    parentsdta$hom <- buildEhomophily(gammatilde) # get Eh_i^t
    parentsdta0$hom <- parentsdta[parentsdta$type==0,"hom"]
    parentsdta1$hom <- parentsdta[parentsdta$type==1,"hom"]
    parentsdta0$het <- 1-parentsdta0$hom
    parentsdta1$het <- 1-parentsdta1$hom

    obser0 <- sample(1:nrow(parentsdta0),replace=T) # bootstrap observations
    obser1 <- sample(1:nrow(parentsdta1),replace=T) # bootstrap observations
    out0 <- flexmix(peffort ~ 0 + I(hom*white) + I(hom*black) + I(hom*hisp) + I(hom*asian) + I(hom*momworks) + I(hom*female) + I(hom*age) + hom:schoolf + het, data=parentsdta0[obser0,], k=2) # regression t=L
    out1 <- flexmix(peffort ~ 0 + I(het*white) + I(het*black) + I(het*hisp) + I(het*asian) + I(het*momworks) + I(het*female) + I(het*age) + het:schoolf + hom, data=parentsdta1[obser1,], k=2) # regression t=H
    coef01 <- parameters(out0,component=1)[,1]
    coef02 <- parameters(out0,component=2)[,1]
    coef11 <- parameters(out1,component=1)[,1]
    coef12 <- parameters(out1,component=2)[,1]
    
    result01[s,1:length(coef01)] <- coef01
    result02[s,1:length(coef02)] <- coef02
    result11[s,1:length(coef11)] <- coef11
    result12[s,1:length(coef12)] <- coef12
    print(s)
  }

  ## identifiability
  wmin <- c(result01[,8]) <= c(result02[,8])
  r01 <- result01
  r01[wmin==F,] <- result02[wmin==F,]
  r02 <- result02
  r02[wmin==F,] <- result01[wmin==F,]
  wmin <- c(result11[,8]) <= c(result12[,8])
  r11 <- result11
  r11[wmin==F,] <- result12[wmin==F,]
  r12 <- result12
  r12[wmin==F,] <- result11[wmin==F,]
  
  rm(D,dta,outdta,S,schdummy,X,G,P,parentsdta,parentsdta0,parentsdta1) # remove confidential data
  save.image("outestim_final_het.RData") # save estimation results
  

  
  
