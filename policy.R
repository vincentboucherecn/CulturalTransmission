################################################################
################################################################
####################  First Policy  ############################
################################################################
################################################################

set.seed(1234)
setwd("~/Dropbox/intertransmission/revisionJOLE")
load("outestim_final_cluster.RData")
library(mvtnorm)
library(maxLik)
library(readstata13)
library(numDeriv)
library(ggplot2)
library(tidyverse)
dta <- read.dta13("/home/vincent/Desktop/tmpah/intercultural/newdata.dta") # read formatted data
dta$dat <- NULL
dta[dta==(-1)] <- NA #recode missing

nsim <- 500 # number of simulations for each swap


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

jlikhes(fulltheta) ## gets optimum values for lambda and s2 (saved as global variables for the optimal parameter values) and P

## re-arrange new vector of parameter values (works even if some groups are dropped (this is why there's is a loop))

t1 <- fulltheta[1:11] # network para [both white.... typeH with typeH]
t2 <- fulltheta[26:33] # socialization para [white....Type]
tL <- as.numeric(PEout0$coefficients[1:8]) # type L parents
tH <- as.numeric(PEout1$coefficients[1:8]) # type H parents
tLt <- tHt <- t1t <- t2t <- rep(0,length(X))
for (i in 1:length(X)){
  t1t[i] <- fulltheta[(11+sid[i])]
  t2t[i] <- fulltheta[(33+sid[i])]
  tLt[i] <- as.numeric(PEout0$coefficients[(8+sid[i])])
  tHt[i] <- as.numeric(PEout1$coefficients[(8+sid[i])])
}
t1 <- c(t1,t1t)
t2 <- c(t2,t2t)
tL <- c(tL,tLt)
tH <- c(tH,tHt)

#########################################3
############ uniform policy #############
#########################################3

collecttot <- NULL
sds <- sd(parentsdta$social) # standard dev of observed socialization efforts

for (str in 0:10){
  print(str)
  ## compute policy shift
  bp <- vector("list",length(X))
  for (i in 1:length(X)){
    bp[[i]] <- matrix((str*sds/10),nrow(X[[i]]),1)
  }
  ## simul counterfactual
  collectt <- simulmore(bp,t2,tL,tH)
  collectt$strength <- str
  collecttot <- rbind(collecttot,collectt)
}
collecttot$strfactor <- factor(collecttot$strength, labels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))

collecttot1 <- collecttot
pdf(file="subsidy.pdf")

collecttot1 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=tau, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Subsidy") + ylab("Equilibrium education effort") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot1 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=h, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Subsidy") + ylab("Equilibrium fraction of same-type friends") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot1 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=s, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Subsidy") + ylab("Equilibrium socialization effort") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot1 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=degree, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Subsidy") + ylab("Equilibrium number of friends") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot1 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=peducated, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Subsidy") + ylab("Probability that the child becomes educated") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

dev.off()
#########################################3
############ more heterophily #############
#########################################3


collecttot <- NULL
split <- (rep(fulltheta[11],2) - fulltheta[9:10])/10
for (str in 0:10){
  print(str)
  bp <- vector("list",length(X))
  for (i in 1:length(X)){
    bp[[i]] <- matrix(0,nrow(X[[i]]),1) #no shift of b
  }
  thetatmp <- fulltheta
  thetatmp[9:10] <- fulltheta[9:10] + str*split # increase heterophily
  jlikhes(thetatmp) ## gets optimum values for lambda and s2 (saved as global variables for the optimal parameter values) and P (preference bias!)
  
  collectt <- simulmore(bp,t2,tL,tH)
  collectt$strength <- str
  collecttot <- rbind(collecttot,collectt)
}
jlikhes(fulltheta) # reset values of global variables
collecttot$strfactor <- factor(collecttot$strength, labels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))

collecttot2 <- collecttot
pdf(file="lessHH.pdf")
collecttot2 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=tau, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Increasing Heterophily") + ylab("Equilibrium education effort") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot2 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=h, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Increasing Heterophily") + ylab("Equilibrium fraction of same-type friends") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot2 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=s, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Increasing Heterophily") + ylab("Equilibrium socialization effort") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot2 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=degree, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Increasing Heterophily") + ylab("Equilibrium number of friends") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot2 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=peducated, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Increasing Heterophily") + ylab("Probability that the child becomes educated") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

dev.off()
#########################################3
################ type=0 #################3
#########################################3

collecttot <- NULL
sds <- sd(parentsdta$social)
for (str in 0:10){
  print(str)
  bp <- vector("list",length(X))
  for (i in 1:length(X)){
    bp[[i]] <- matrix((str*sds/10),nrow(X[[i]]),1)*(1-X[[i]][,8])
  }
  collectt <- simulmore(bp,t2,tL,tH)
  collectt$strength <- str
  collecttot <- rbind(collecttot,collectt)
}
collecttot$strfactor <- factor(collecttot$strength, labels=c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"))

collecttot3 <- collecttot

pdf(file="subsidy0.pdf")

collecttot3 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=tau, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Subsidy - low income") + ylab("Equilibrium education effort") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot3 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=h, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Subsidy - low type") + ylab("Equilibrium fraction of same-type friends") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot3 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=s, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Subsidy - low type") + ylab("Equilibrium socialization effort") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot3 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=degree, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Subsidy - low type") + ylab("Equilibrium number of friends") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

collecttot3 %>% 
  dplyr::filter(type %in% c(0,1)) %>%
  ggplot(aes(x=strfactor, y=peducated, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_boxplot() + theme_minimal() + xlab("Subsidy - low type") + ylab("Probability that the child becomes educated") +scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") + theme(text = element_text(size = 15))  

dev.off()


###############################################
############# Regressions #####################
###############################################

rs1 <- lm(s ~ 0 + type + I(1-type) + I(strength*type) + I(strength*(1-type)), data=collecttot1)
rh1 <- lm(h ~ 0 + type + I(1-type) + I(strength*type) + I(strength*(1-type)), data=collecttot1)
red1 <- lm(peducated ~ 0 + type + I(1-type) + I(strength*type) + I(strength*(1-type)), data=collecttot1)
rtau1 <- lm(tau ~ 0 + type + I(1-type) + I(strength*type) + I(strength*(1-type)), data=collecttot1)

rs2 <- lm(s ~ 0 + type + I(1-type) + I(strength*type) + I(strength*(1-type)), data=collecttot3)
rh2 <- lm(h ~ 0 + type + I(1-type) + I(strength*type) + I(strength*(1-type)), data=collecttot3)
red2 <- lm(peducated ~ 0 + type + I(1-type) + I(strength*type) + I(strength*(1-type)), data=collecttot3)
rtau2 <- lm(tau ~ 0 + type + I(1-type) + I(strength*type) + I(strength*(1-type)), data=collecttot3)


###############################################
############ Initial graphs ###################
###############################################

theme_set(theme_minimal() + theme(legend.position = c(0.8,0.8)))
p <- ggplot(collecttot1[collecttot1$strength==0,], aes(x=s, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=20) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") +
  scale_x_continuous(name = "Ex-ante simulated socialization efforts") +
  scale_y_continuous(name = "Number of Students") + theme(text = element_text(size = 15))  
plot(p)

p <- ggplot(collecttot1[collecttot1$strength==0,], aes(x=h, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=20) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") +
  scale_x_continuous(name = "Ex-ante simulated fraction of same-type links") +
  scale_y_continuous(name = "Number of Students") + theme(text = element_text(size = 15))  
plot(p)


p <- ggplot(collecttot1[collecttot1$strength==0,], aes(x=tau, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=20) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") +
  scale_x_continuous(name = "Ex-ante simulated education effort") +
  scale_y_continuous(name = "Number of Parents") + theme(text = element_text(size = 15))  
plot(p)

p <- ggplot(collecttot1[collecttot1$strength==0,], aes(x=peducated, fill=factor(type, labels=c('Lowly educated','Highly educated')))) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=20) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) + labs(fill="") +
  scale_x_continuous(name = "Ex-ante simulated probabilities that the student becomes educated") +
  scale_y_continuous(name = "Number of Students") + theme(text = element_text(size = 15))  
plot(p)


rm(D,dta,outdta,S,schdummy,X,G,P,parentsdta,parentsdta0,parentsdta1) # remove confidential data
save.image("total_withpolicy.RData")

