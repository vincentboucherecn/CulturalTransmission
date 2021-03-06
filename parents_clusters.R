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
  
  ggplot(parentsdta, aes(peffort)) + geom_histogram(fill="#69b3a2",col="white") + theme_minimal() + xlab("Education Effort") + ylab("") + labs(fill="")
  
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
  out0 <- lm(peffort ~ 0 + I(hom*white) + I(hom*black) + I(hom*hisp) + I(hom*asian) + I(hom*momworks) + I(hom*female) + I(hom*age) + hom:schoolf + het, data=parentsdta0)
  out1 <- lm(peffort ~ 0 + I(het*white) + I(het*black) + I(het*hisp) + I(het*asian) + I(het*momworks) + I(het*female) + I(het*age) + het:schoolf + hom, data=parentsdta1)
  PEout0 <- out0 # backup estimate
  PEout1 <- out1 # backup estimate

    
  ## compute the predicted derivative of tau with respect to h
  margin0 <- out0$coefficients["I(hom * white)"]*parentsdta0$white + out0$coefficients["I(hom * black)"]*parentsdta0$black +
    out0$coefficients["I(hom * hisp)"]*parentsdta0$hisp + out0$coefficients["I(hom * asian)"]*parentsdta0$asian + out0$coefficients["I(hom * momworks)"]*parentsdta0$momworks + out0$coefficients["I(hom * female)"]*parentsdta0$female +
    out0$coefficients["I(hom * age)"]*parentsdta0$age + out0$coefficients["hom:schoolf2"]*as.numeric(parentsdta0$school==2) +
    out0$coefficients["hom:schoolf3"]*as.numeric(parentsdta0$school==3) + out0$coefficients["hom:schoolf7"]*as.numeric(parentsdta0$school==7) +
    out0$coefficients["hom:schoolf8"]*as.numeric(parentsdta0$school==8) + out0$coefficients["hom:schoolf28"]*as.numeric(parentsdta0$school==28) +
    out0$coefficients["hom:schoolf58"]*as.numeric(parentsdta0$school==58) + out0$coefficients["hom:schoolf77"]*as.numeric(parentsdta0$school==77) +
    out0$coefficients["hom:schoolf81"]*as.numeric(parentsdta0$school==81) + out0$coefficients["hom:schoolf88"]*as.numeric(parentsdta0$school==88) +
    out0$coefficients["hom:schoolf106"]*as.numeric(parentsdta0$school==106) + out0$coefficients["hom:schoolf126"]*as.numeric(parentsdta0$school==126) +
    out0$coefficients["hom:schoolf175"]*as.numeric(parentsdta0$school==175) + out0$coefficients["hom:schoolf194"]*as.numeric(parentsdta0$school==194) +
    out0$coefficients["hom:schoolf369"]*as.numeric(parentsdta0$school==369) - out0$coefficients["het"]
  
  margin1 <- - out1$coefficients["I(het * white)"]*parentsdta1$white - out1$coefficients["I(het * black)"]*parentsdta1$black -
    out1$coefficients["I(het * hisp)"]*parentsdta1$hisp - out1$coefficients["I(het * asian)"]*parentsdta1$asian - out1$coefficients["I(het * momworks)"]*parentsdta1$momworks - out1$coefficients["I(het * female)"]*parentsdta1$female -
    out1$coefficients["I(het * age)"]*parentsdta1$age - out1$coefficients["het:schoolf2"]*as.numeric(parentsdta1$school==2) -
    out1$coefficients["het:schoolf3"]*as.numeric(parentsdta1$school==3) - out1$coefficients["het:schoolf7"]*as.numeric(parentsdta1$school==7) -
    out1$coefficients["het:schoolf8"]*as.numeric(parentsdta1$school==8) - out1$coefficients["het:schoolf28"]*as.numeric(parentsdta1$school==28) -
    out1$coefficients["het:schoolf58"]*as.numeric(parentsdta1$school==58) - out1$coefficients["het:schoolf77"]*as.numeric(parentsdta1$school==77) -
    out1$coefficients["het:schoolf81"]*as.numeric(parentsdta1$school==81) - out1$coefficients["het:schoolf88"]*as.numeric(parentsdta1$school==88) -
    out1$coefficients["het:schoolf106"]*as.numeric(parentsdta1$school==106) - out1$coefficients["het:schoolf126"]*as.numeric(parentsdta1$school==126) -
    out1$coefficients["het:schoolf175"]*as.numeric(parentsdta1$school==175) - out1$coefficients["het:schoolf194"]*as.numeric(parentsdta1$school==194) -
    out1$coefficients["het:schoolf369"]*as.numeric(parentsdta1$school==369) + out1$coefficients["hom"]
  
   ######
    #Bootstrap SE
   ######
  
  id0 <- cbind(1:nrow(parentsdta0),parentsdta0$scid)
  lst <- unique(id0[,2])
  s0 <- rep(NA,length(D))
  for (i in 1:length(D)){
    s0[i] <- length(id0[id0[,2]==lst[i],2])
  }
  id1 <- cbind(1:nrow(parentsdta1),parentsdta1$scid)
  s1 <- rep(NA,length(D))
  for (i in 1:length(D)){
    s1[i] <- length(id1[id1[,2]==lst[i],2])
  }
  

  result0 <- result1 <- matrix(NA,bsim,30) # initialize results matrices
  for (s in 1:bsim){
    gammatilde <- rmvnorm(1,fulltheta,Vtheta) # draw coefficient
    parentsdta$hom <- buildEhomophily(gammatilde) # get Eh_i^t
    parentsdta0$hom <- parentsdta[parentsdta$type==0,"hom"]
    parentsdta1$hom <- parentsdta[parentsdta$type==1,"hom"]
    parentsdta0$het <- 1-parentsdta0$hom
    parentsdta1$het <- 1-parentsdta1$hom

    obser0 <- NULL
    obser1 <- NULL
    for (i in 1:length(D)){
      gr0 <- sample(1:length(D),1) # sample a group
      grid0 <- lst[gr0] # group id
      grs0 <- s0[gr0] # group size
      if (grs0==0){
        ind0 <- NULL
      } else if (grs0==1) {
        ind0 <- id0[id0[,2]==grid0,1] #
      } else {
        ind0 <- sample(id0[id0[,2]==grid0,1],grs0) # sample with replacement
      }
      obser0 <- c(obser0,ind0)
      
      gr1 <- sample(1:length(D),1) # sample a group
      grid1 <- lst[gr1] # group id
      grs1 <- s1[gr1] # group size
      if (grs1==0){
        ind1 <- NULL
      } else if (grs1==1) {
        ind1 <- id1[id1[,2]==grid1,1] #
      } else {
        ind1 <- sample(id1[id1[,2]==grid1,1],grs1) # sample with replacement
      }
      obser1 <- c(obser1,ind1)
    }
    
    
    out0 <- lm(peffort ~ 0 + I(hom*white) + I(hom*black) + I(hom*hisp) + I(hom*asian) + I(hom*momworks) + I(hom*female) + I(hom*age) + hom:schoolf + het, data=parentsdta0[obser0,]) # regression t=L
    out1 <- lm(peffort ~ 0 + I(het*white) + I(het*black) + I(het*hisp) + I(het*asian) + I(het*momworks) + I(het*female) + I(het*age) + het:schoolf + hom, data=parentsdta1[obser1,]) # regression t=H
    result0[s,1:length(out0$coefficients)] <- out0$coefficients
    result1[s,1:length(out1$coefficients)] <- out1$coefficients
    print(s)
  }

  rm(D,dta,outdta,S,schdummy,X,G,P,parentsdta,parentsdta0,parentsdta1) # remove confidential data
  save.image("outestim_final_cluster.RData") # save estimation results
  
  ##### graphics #####
  
  library(ggplot2)
  dtaggplot <- as.data.frame(rbind(cbind(margin0,rep(0,length(margin0))), cbind(margin1,rep(1,length(margin1)))))
  dtaggplot$V2 <- factor(dtaggplot$V2, labels=c('Lowly educated','Highly educated'))
  l0 <- c(quantile(dtaggplot[dtaggplot$V2=="Lowly educated","margin0"])[2])
  l1 <- c(quantile(dtaggplot[dtaggplot$V2=="Lowly educated","margin0"])[4])
  h0 <- c(quantile(dtaggplot[dtaggplot$V2=="Highly educated","margin0"])[2])
  h1 <- c(quantile(dtaggplot[dtaggplot$V2=="Highly educated","margin0"])[4])
  dtaggplotI <- dtaggplot
  dtaggplotI$V2 <- as.character(dtaggplotI$V2)
  dtaggplotI[dtaggplotI$margin0>l0 & dtaggplotI$margin0<l1 & dtaggplotI$V2=="Lowly educated","V2"] <- "Lowly educated - IQ"
  dtaggplotI[dtaggplotI$margin0>h0 & dtaggplotI$margin0<h1 & dtaggplotI$V2=="Highly educated","V2"] <- "Highly educated - IQ"
  dtaggplotI$V2 <- as.factor(dtaggplotI$V2)
  dtaggplotI$V2 <- factor(dtaggplotI$V2, levels=c("Lowly educated"    ,   "Lowly educated - IQ","Highly educated"   ,   "Highly educated - IQ" ))
  theme_set(theme_minimal() + theme(legend.position = c(0.8,0.8)))
  p<- ggplot(dtaggplotI, aes(x=margin0, fill=V2)) +
    geom_histogram( color="#e9ecef", position = 'identity',bins=20) +
    scale_fill_manual(values=c(alpha("#69b3a2",0.3),alpha("#69b3a2",0.8),alpha("#404080",0.3),alpha("#404080",0.8))) + labs(fill="") +
    scale_x_continuous(name = "Derivative of effort w.r.t. child's homophily") +
    scale_y_continuous(name = "Number of Parents") + geom_vline(xintercept=0)
    plot(p)
    p<- ggplot(dtaggplot, aes(x=margin0, fill=V2)) +
      geom_histogram( color="#e9ecef", position = 'identity',bins=20,alpha=0.6) +
      scale_fill_manual(values=c("#69b3a2","#404080")) + labs(fill="") +
      scale_x_continuous(name = "Derivative of effort w.r.t. child's homophily") +
      scale_y_continuous(name = "Number of Parents") + geom_vline(xintercept=0) + theme(text = element_text(size = 15))  
    plot(p)
  
