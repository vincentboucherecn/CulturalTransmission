
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


###############################################
############ Initial graphs ###################
###############################################

library(ggplot2)
library(tidyverse)

parentsdta$h <- rh

theme_set(theme_minimal() + theme(legend.position = c(0.8,0.8)))

p <- ggplot(parentsdta, aes(x=social, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=20) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") +
  scale_x_continuous(name = "Socialization efforts (data)") +
  scale_y_continuous(name = "Number of Students") + theme(text = element_text(size = 15))  
plot(p)

p <- ggplot(parentsdta, aes(x=peffort, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=20) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") +
  scale_x_continuous(name = "Education effort (data)") +
  scale_y_continuous(name = "Number of Parents") + theme(text = element_text(size = 15))  
plot(p)

p <- ggplot(parentsdta, aes(x=h, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=20) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") +
  scale_x_continuous(name = "Fraction of same-type links (data)") +
  scale_y_continuous(name = "Number of Students") + theme(text = element_text(size = 15))  
plot(p)

####peduc

parentsdta$pe <- NA
parentsdta$pe[parentsdta$type==0] <- parentsdta$peffort[parentsdta$type==0] + (1-parentsdta$peffort[parentsdta$type==0])*(1-parentsdta$h[parentsdta$type==0])
parentsdta$pe[parentsdta$type==1] <- parentsdta$peffort[parentsdta$type==1] + (1-parentsdta$peffort[parentsdta$type==1])*parentsdta$h[parentsdta$type==1]

p <- ggplot(parentsdta, aes(x=pe, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=20) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") +
  scale_x_continuous(name = "Probabilities that the student becomes educated (simulated)") +
  scale_y_continuous(name = "Number of Students") + theme(text = element_text(size = 15))  
plot(p)

creategraph <- as.data.frame(c("s","tau","h","piH"))
colnames(creategraph) <- "Variable"
creategraph$type <- 0
creategraph$mean <- colMeans(parentsdta[parentsdta$type==0,c("social","peffort","h","pe")])
creategraph$sd <- apply(parentsdta[parentsdta$type==0,c("social","peffort","h","pe")],2,sd)
creategraph0 <- creategraph

creategraph <- as.data.frame(c("s","tau","h","piH"))
colnames(creategraph) <- "Variable"
creategraph$type <- 1
creategraph$mean <- colMeans(parentsdta[parentsdta$type==1,c("social","peffort","h","pe")])
creategraph$sd <- apply(parentsdta[parentsdta$type==1,c("social","peffort","h","pe")],2,sd)

creategraph <- rbind(creategraph0,creategraph)
creategraph$data <- 1
creategraph1 <- creategraph

load("total_withpolicy.RData")
simdata <- collecttot1[collecttot1$strength==0,]

creategraph <- as.data.frame(c("s","tau","h","piH"))
colnames(creategraph) <- "Variable"
creategraph$type <- 0
creategraph$mean <- colMeans(simdata[simdata$type==0,c("s","tau","h","peducated")])
creategraph$sd <- apply(simdata[simdata$type==0,c("s","tau","h","peducated")],2,sd)
creategraph0 <- creategraph

creategraph <- as.data.frame(c("s","tau","h","piH"))
colnames(creategraph) <- "Variable"
creategraph$type <- 1
creategraph$mean <- colMeans(simdata[simdata$type==1,c("s","tau","h","peducated")])
creategraph$sd <- apply(simdata[simdata$type==1,c("s","tau","h","peducated")],2,sd)
creategraph <- rbind(creategraph0,creategraph)
creategraph$data <- 0

creategraph <- rbind(creategraph1,creategraph)

ggplot(creategraph[creategraph$data==1,], aes(x=factor(Variable),y=mean,group=factor(type),color=factor(type))) + geom_point(position=position_dodge(width=0.5)) +
  geom_pointrange(aes(x=factor(Variable),ymax=mean+sd,ymin=mean-sd),position=position_dodge(width=0.5))

ggplot(creategraph[creategraph$data==0,], aes(x=factor(Variable),y=mean,group=factor(type),color=factor(type))) + geom_point(position=position_dodge(width=0.5)) +
  geom_pointrange(aes(x=factor(Variable),ymax=mean+sd,ymin=mean-sd),position=position_dodge(width=0.5))

creategraph$type <- factor(creategraph$type, labels=c('Lowly educated','Highly educated'))
creategraph$data <- factor(creategraph$data, labels=c('Simulated','Data'))

ggplot(creategraph[creategraph$Variable=="s",], aes(x=(data),y=mean,group=type,color=type)) + geom_point(position=position_dodge(width=0.5)) +
  geom_pointrange(aes(x=(data),ymax=mean+sd,ymin=mean-sd),position=position_dodge(width=0.5)) + xlab("Simulated values vs Data") +
  ylab("Socialization effort (child)") + scale_color_manual(values=c("#69b3a2", "#404080")) + theme_minimal() + theme(text = element_text(size = 15))  
ggsave("~/Dropbox/intertransmission/2ndRevisionJOLE/fits.pdf")

ggplot(creategraph[creategraph$Variable=="tau",], aes(x=(data),y=mean,group=type,color=type)) + geom_point(position=position_dodge(width=0.5)) +
  geom_pointrange(aes(x=(data),ymax=mean+sd,ymin=mean-sd),position=position_dodge(width=0.5)) + xlab("Simulated values vs Data") +
  ylab("Education effort (parent)") + scale_color_manual(values=c("#69b3a2", "#404080")) + theme_minimal() + theme(text = element_text(size = 15))  
ggsave("~/Dropbox/intertransmission/2ndRevisionJOLE/fittau.pdf")

ggplot(creategraph[creategraph$Variable=="h",], aes(x=(data),y=mean,group=type,color=type)) + geom_point(position=position_dodge(width=0.5)) +
  geom_pointrange(aes(x=(data),ymax=mean+sd,ymin=mean-sd),position=position_dodge(width=0.5)) + xlab("Simulated values vs Data") +
  ylab("Fraction of same-type links") + scale_color_manual(values=c("#69b3a2", "#404080")) + theme_minimal() + theme(text = element_text(size = 15))  
ggsave("~/Dropbox/intertransmission/2ndRevisionJOLE/fith.pdf")

ggplot(creategraph[creategraph$Variable=="piH",], aes(x=(data),y=mean,group=type,color=type)) + geom_point(position=position_dodge(width=0.5)) +
  geom_pointrange(aes(x=(data),ymax=mean+sd,ymin=mean-sd),position=position_dodge(width=0.5)) + xlab("Simulated values vs Data") +
  ylab("Probability of high education") + scale_color_manual(values=c("#69b3a2", "#404080")) + theme_minimal() + theme(text = element_text(size = 15))  
ggsave("~/Dropbox/intertransmission/2ndRevisionJOLE/fitpi.pdf")
