################################################################
################################################################
################### Children's model ###########################
################################################################
################################################################
set.seed(1234)
library(maxLik)
library(readstata13)
library(numDeriv)
library(ggplot2)
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

fractypes <- rep(NA,length(X))
for (i in 1:length(X)){
  fractypes[i] <- mean(X[[i]][,8])
}
fractypes <- as.data.frame(fractypes)

ggplot(fractypes, aes(fractypes)) + geom_histogram(fill="#69b3a2",col="white") + theme_minimal() + xlab("Fraction of High type") + ylab("Number of groups") + labs(fill="")

ggplot(outdta, aes(si)) + geom_histogram(fill="#69b3a2",col="white") + theme_minimal() + xlab("Socialisation Levels") + ylab("") + labs(fill="")
