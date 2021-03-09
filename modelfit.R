
set.seed(1234)
setwd("~/Dropbox/intertransmission/revisionJOLE")
load("outestim_partial.RData")
library(mvtnorm)
library(maxLik)
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

fith <- matrix(NA,nrow(parentsdta),nsim)
for (sim in 1:nsim){
  fith[,sim] <- buildEhomophily_single(fulltheta)
}
rh <- buildrealh()
tst0 <- tst1 <- tst <- ttst <- rep(NA,nsim)

for (sim in 1:nsim){
  ch <- fith[,sim]
  ttst[sim] <- mean(ch)
  tst[sim] <- mean(ch[ch>0 & ch<1])
  tst0[sim] <- mean(ch==0)
  tst1[sim] <- mean(ch==1)
}
mean(tst)
sd(tst)
mean(rh[rh>0 & rh<1])
mean(tst0)
sd(tst0)
mean(tst1)
sd(tst1)
mean(rh==0)
mean(rh==1)

mean(ttst)
sd(ttst)
mean(rh)
