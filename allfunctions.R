################################################################
################################################################
################### Children's model ###########################
################################################################
################################################################


################################################################
###################### Format data #############################
################################################################

buildprobadta <- function(){
  ### function to format and organize data for the joint likelihood P(G,s)
  
  ### keep only relevant variables and clean data
  dtaProba <- dta[,c(which(colnames(dta)=="aid"):which(colnames(dta)=="h1rm5"),which(colnames(dta)=="h1rf1"),which(colnames(dta)=="h1nb1"):which(colnames(dta)=="h1nb7"),
                     which(colnames(dta)=="s44a1"):which(colnames(dta)=="s44"))]
  dtaProba <- dtaProba[dtaProba$h1gi20<=12,] # drop individuals with no grade levels
  dtaProba <- dtaProba[((dtaProba$scid != 1) & (dtaProba$scid != 115)),]
  ## recode/clean variables, missing values, individual characteristics
  dtaProba$female <- dtaProba$bio_sex-1 #female
  dtaProba$bio_sex <- NULL
  dtaProba$age <- 95 - (dtaProba$h1gi1y+(dtaProba$h1gi1m-1)/12) #age
  dtaProba$age <- pmin(dtaProba$age,18)
  dtaProba$age <- pmax(dtaProba$age,12)
  dtaProba$h1gi1m <- NULL
  dtaProba$h1gi1y <- NULL
  dtaProba$h1gi3 <- NULL
  dtaProba[dtaProba$h1gi4>1,"h1gi4"] <- 0
  colnames(dtaProba)[which(colnames(dtaProba)=="h1gi4")] <- "hisp"
  dtaProba[dtaProba$h1gi6a>1,"h1gi6a"] <- 0
  colnames(dtaProba)[which(colnames(dtaProba)=="h1gi6a")] <- "white"
  dtaProba[dtaProba$h1gi6b>1,"h1gi6b"] <- 0
  colnames(dtaProba)[which(colnames(dtaProba)=="h1gi6b")] <- "black"
  dtaProba$h1gi6c <- NULL
  dtaProba[dtaProba$h1gi6d>1,"h1gi6d"] <- 0
  colnames(dtaProba)[which(colnames(dtaProba)=="h1gi6d")] <- "asian"
  dtaProba$h1gi6e <- NULL
  dtaProba$h1rm4 <- NULL
  dtaProba[dtaProba$h1rm5>1,"h1rm5"] <- 0
  colnames(dtaProba)[which(colnames(dtaProba)=="h1rm5")] <- "momworks"
  dtaProba$type <- as.numeric((dtaProba$h1rf1>=8 & dtaProba$h1rf1<=9) | (dtaProba$h1rm1>=8 & dtaProba$h1rm1<=9 ) ) # type 1 if either parent is college educated
  dtaProba$h1rf1 <- NULL
  dtaProba$h1rm1 <- NULL
  
  ## socialization effort
  dtaProba[dtaProba$h1da4>3,"h1da4"] <- NA
  dtaProba[dtaProba$h1da5>3,"h1da5"] <- NA
  dtaProba[dtaProba$h1da6>3,"h1da6"] <- NA
  dtaProba[dtaProba$h1da7>3,"h1da7"] <- NA
  dtaProba[dtaProba$h1nb1>1,"h1nb1"] <- 0
  dtaProba[dtaProba$h1nb4>1,"h1nb4"] <- 0
  dtaProba$schoolact <- rowSums(dtaProba[,c(which(colnames(dtaProba)=="s44a1"):which(colnames(dtaProba)=="s44a33"))],na.rm=T) # number of activities
  
  dtaProba$schoolact <- (pmin(dtaProba$schoolact,5)+1)/6
  dtaProba[,c(which(colnames(dtaProba)=="s44a1"):which(colnames(dtaProba)=="s44"))] <- NULL
  dtaProba$dailyact <- rowMeans(dtaProba[,c(which(colnames(dtaProba)=="h1da4"):which(colnames(dtaProba)=="h1da7"))],na.rm=T) # average daily activities
  dtaProba[,c(which(colnames(dtaProba)=="h1da4"):which(colnames(dtaProba)=="h1da7"))] <- NULL
  dtaProba$neighact <- rowMeans(dtaProba[,c("h1nb1","h1nb4")],na.rm=T) # average neighbourhood participation
  dtaProba[,c(which(colnames(dtaProba)=="h1nb1"):which(colnames(dtaProba)=="h1nb7"))] <- NULL
  
  dtaProba$social <- rowMeans(dtaProba[,c("schoolact","dailyact","neighact")],na.rm=T) # socialization index
  dtaProba$social <- pmin(dtaProba$social,quantile(dtaProba$social,0.99))/quantile(dtaProba$social,0.99) # trim top 1% and normalize between 0 and 1
  dtaProba$social <- log(dtaProba$social+1)/max(log(dtaProba$social+1)) # log scale
  #dtaProba[,c("schoolact","dailyact","neighact")] <- NULL
  
  
  ## build data for estimation procedures
  
  ## recode school names to 1:nschools
  dtaProba$scid_bkp <- dtaProba$scid
  sl <- unique(dtaProba$scid_bkp)
  rsl <- 1:length(sl)
  for (i in 1:length(sl)){
    dtaProba$scid_bkp[dtaProba$scid==sl[i]] <- rsl[i]
  }
  
  ## create school-grade index
  dtaProba$scid <- dtaProba$scid_bkp*100 + dtaProba$h1gi20
  listsch <- unique(dtaProba$scid) # list of schools
  lM <- length(listsch) # number of schools/grades
  
  # initialize lists
  Dlist <- vector("list",lM)
  Xlist <- vector("list",lM)
  Slist <- vector("list",lM)
  Glist <- vector("list",lM)
  Poslist <- vector("list",lM)
  nlist <- rep(0,lM)
  sidlist <- rep(0,lM)
  lstfr <- c(which(colnames(dtaProba)=="mf_aid1"):which(colnames(dtaProba)=="ff_aid5")) # list of friends' id
  
  # compute total number of pairs of students (including self links, which will be removed later)
  nobs <- 0
  for (i in 1:lM){
    nobs <- nobs + sum(as.numeric(dtaProba$scid==listsch[i]))^2
  }
  
  dtaset <- as.data.frame(matrix(NA,nobs,18)) # initialize dataset of correct size for dyadic regression
  colnames(dtaset) <- c("g","white","black","hisp","asian","mwork","gender","age","dist","typeLH","typeHL","typeHH","social","school","self","si","sj","n")
  
  ## for each school
  pos <- 1
  for (i in 1:lM){
    print(i)
    
    ## keep only students from this school
    tdta <- dtaProba[dtaProba$scid==listsch[i],] # temporary dataset for school i
    nt <- nrow(tdta) # number of students
    if (nt<10){
      print("**** check n<10 ****")
    }
    nlist[i] <- nt
    sidlist[i] <- tdta[1,"scid_bkp"]
    Xlist[[i]] <- as.matrix(tdta[,c("white","black","hisp","asian","momworks","female","age","type")]) # store expl var for school i
    Slist[[i]] <- matrix(tdta$social,nt,1) # store socialization efforts
    Poslist[[i]] <- as.matrix(tdta[,c("x","y","h1gi20")])
    ## generate network
    Gt <- matrix(0,nt,nt)
    for (j in 1:nt){
      for (k in lstfr){
        if (!is.na(tdta[j,k])){
          fr <- which(tdta$aid==tdta[j,k])
          if (length(fr)>0){
            Gt[j,fr[1]] <- 1
          }
        }
      }
    }
    diag(Gt) <- 0 # zero on diagonal
    Glist[[i]] <- Gt # store observed network
    print(mean(rowSums(Gt)))
    
    tDlist <- list("vector",11) # for each explanatory variable
    dtaset[pos:(pos+nt^2-1),"g"] <- c(Gt)
    
    tDlist[[1]] <- (matrix(rep(tdta[,"white"],nt),nt,nt)*t(matrix(rep(tdta[,"white"],nt),nt,nt))) ## i*j = 1 if i==j==1 "Both white"
    dtaset[pos:(pos+nt^2-1),"white"] <- c(tDlist[[1]])
    tDlist[[2]] <- (matrix(rep(tdta[,"black"],nt),nt,nt)*t(matrix(rep(tdta[,"black"],nt),nt,nt))) ## i*j = 1 if i==j==1 "Both black"
    dtaset[pos:(pos+nt^2-1),"black"] <- c(tDlist[[2]])
    tDlist[[3]] <- (matrix(rep(tdta[,"hisp"],nt),nt,nt)*t(matrix(rep(tdta[,"hisp"],nt),nt,nt))) ## i*j = 1 if i==j==1 "Both hisp"
    dtaset[pos:(pos+nt^2-1),"hisp"] <- c(tDlist[[3]])
    tDlist[[4]] <- (matrix(rep(tdta[,"asian"],nt),nt,nt)*t(matrix(rep(tdta[,"asian"],nt),nt,nt))) ## i*j = 1 if i==j==1 "Both asian"
    dtaset[pos:(pos+nt^2-1),"asian"] <- c(tDlist[[4]])
    tDlist[[5]] <- (matrix(rep(tdta[,"momworks"],nt),nt,nt)*t(matrix(rep(tdta[,"momworks"],nt),nt,nt))) ## i*j = 1 if i==j==1 "Both moms work"
    dtaset[pos:(pos+nt^2-1),"mwork"] <- c(tDlist[[5]])
    tDlist[[6]] <- matrix(as.numeric(matrix(rep(tdta[,"female"],nt),nt,nt)==t(matrix(rep(tdta[,"female"],nt),nt,nt))),nt,nt) ## "same gender"
    dtaset[pos:(pos+nt^2-1),"gender"] <- c(tDlist[[6]])
    tDlist[[7]] <- (abs(matrix(rep(tdta[,"age"],nt),nt,nt)-t(matrix(rep(tdta[,"age"],nt),nt,nt)))) ## "age difference"
    dtaset[pos:(pos+nt^2-1),"age"] <- c(tDlist[[7]])
    tDlist[[8]] <- (  sqrt( (matrix(rep(tdta[,"x"],nt),nt,nt)-t(matrix(rep(tdta[,"x"],nt),nt,nt)))^2 +
                              (matrix(rep(tdta[,"y"],nt),nt,nt)-t(matrix(rep(tdta[,"y"],nt),nt,nt)))^2 ) ) ## "geo distance"
    dtaset[pos:(pos+nt^2-1),"dist"] <- c(tDlist[[8]])
    tDlist[[9]] <- matrix(as.numeric(matrix(rep(tdta[,"type"],nt),nt,nt)==0)*as.numeric(t(matrix(rep(tdta[,"type"],nt),nt,nt))==1),nt,nt) # type LH
    dtaset[pos:(pos+nt^2-1),"typeLH"] <- c(tDlist[[9]])
    tDlist[[10]] <- matrix(as.numeric(matrix(rep(tdta[,"type"],nt),nt,nt)==1)*as.numeric(t(matrix(rep(tdta[,"type"],nt),nt,nt))==0),nt,nt) # type HL
    dtaset[pos:(pos+nt^2-1),"typeHL"] <- c(tDlist[[10]])
    tDlist[[11]] <- matrix(as.numeric(matrix(rep(tdta[,"type"],nt),nt,nt)==1)*as.numeric(t(matrix(rep(tdta[,"type"],nt),nt,nt))==1),nt,nt) # type HH
    dtaset[pos:(pos+nt^2-1),"typeHH"] <- c(tDlist[[11]])
    
    Dlist[[i]] <- tDlist # store matrices of explanatory variables
    dtaset[pos:(pos+nt^2-1),"social"] <- c(matrix(rep(tdta[,"social"],nt),nt,nt)*t(matrix(rep(tdta[,"social"],nt),nt,nt))) ## "social"
    dtaset[pos:(pos+nt^2-1),"school"] <- c(matrix(rep(tdta[,"scid_bkp"],nt),nt,nt))
    dtaset[pos:(pos+nt^2-1),"group"] <- c(matrix(rep(tdta[,"scid"],nt),nt,nt))
    dtaset[pos:(pos+nt^2-1),"self"] <- c(diag(nt)) # flag ==1 if self link
    dtaset[pos:(pos+nt^2-1),"si"] <- c(matrix(rep(tdta[,"social"],nt),nt,nt)) ## "social"
    dtaset[pos:(pos+nt^2-1),"sj"] <- c(t(matrix(rep(tdta[,"social"],nt),nt,nt))) ## "social"
    dtaset[pos:(pos+nt^2-1),"n"] <- nt
    pos <- pos + nt^2
  }
  return(list(dtaset,Xlist,Dlist,Slist,Glist,nlist,sidlist,Poslist))
}

droppsmallgroups <- function(small){
  
  ## this function removes small groups
  lM <- sum(as.numeric(Nvec>=small))
  Dlist <- vector("list",lM)
  Xlist <- vector("list",lM)
  Slist <- vector("list",lM)
  Glist <- vector("list",lM)
  Poslist <- vector("list",lM)
  nlist <- rep(0,lM)
  sidlist <- rep(0,lM)
  pos <- 1
  for (i in 1:length(D)){
    if (Nvec[i]>=small){
      Dlist[[pos]] <- D[[i]]
      Xlist[[pos]] <- X[[i]]
      Glist[[pos]] <- G[[i]]
      Slist[[pos]] <- S[[i]]
      nlist[[pos]] <- Nvec[[i]]
      sidlist[[pos]] <- sid[[i]]
      if (exists("Posdta")){
        Poslist[[pos]] <- Posdta[[i]]
      }
      pos <- pos + 1
    }
  }
  dtaset <- outdta[outdta$n>=small,]
  return(list(dtaset,Xlist,Dlist,Slist,Glist,nlist,sidlist,Poslist))
}


################################################################
###################### Likelihoods #############################
################################################################

################
### likelihood for the pairwise regression P(G|s)
################

objproba <- function(theta){
  ## conditional likelihood P(G|s)
  print(theta)
  rho <- pnorm(as.matrix(outdta[,2:12])%*%matrix(theta[1:11],11,1)+schdummy%*%matrix(theta[12:(12+ncol(schdummy)-1)],ncol(schdummy),1))
  proba <- outdta[,"social"]*rho
  lik <- outdta[,"g"]*log(proba) + (1-outdta[,"g"])*log(1-proba)
  return(sum(lik,na.rm = T))
}


gradlik <- function(theta){
  ## gradient of the conditional likelihood P(G|s)
  rho <- as.matrix(outdta[,2:12])%*%matrix(theta[1:11],11,1)+schdummy%*%matrix(theta[12:(12+ncol(schdummy)-1)],ncol(schdummy),1)
  rhop <- dnorm(rho)
  rho <- pnorm(rho)
  proba <- outdta[,"social"]*rho
  mbeta <- ((outdta[,"g"]-proba)/(proba*(1-proba)))*outdta[,"social"]*rhop
  outgrad <- matrix(NA,N,length(theta))
  for (i in 1:11){
    outgrad[,i] <- mbeta*outdta[,(i+1)]
  }
  for (i in 1:ncol(schdummy)){
    outgrad[,(i+11)] <- mbeta*schdummy[,i]
  }
  return(colSums(outgrad,na.rm = T))
}


################
### likelihood for the SAR P(s)
################

liksar <- function(lambda,theta){
  ## likelihood P(s), concentrated around lambda (=phi in the paper)
  thetaprob <- theta # parameters for P(G|s)
  lambda <- lambda

  ## initialize variables  
  P <<- list("vector",length(D))
  kx <- ncol(X[[1]]) # number of explanatory variables
  XX <- matrix(0,(kx+ncol(schdummy)),(kx+ncol(schdummy)))
  XMS <-  matrix(0,(kx+ncol(schdummy)),1)

  ## for each school
  for (i in 1:length(D)){
    
    ## build probability matrix
    nt <- nrow(X[[i]])
    Pt <- matrix(0,nt,nt)
    
    ## build expected network
    for (j in 1:length(D[[i]])){
      Pt <- Pt + thetaprob[j]*D[[i]][[j]] # add contribution of expl var j
    }
    Pt <- pnorm(Pt + thetaprob[(length(D[[i]])+sid[i])]) # add school dummy and compute proba
    diag(Pt) <- 0
    P[[i]] <<- Pt # store proba (=d_ij in paper's notation)
    
    ## temporary matrix to compute beta
    Mt <- diag(nt)-lambda*Pt
    Xt <- cbind(X[[i]],matrix(0,nt,ncol(schdummy))) # matrix of expl var, including school dummies
    Xt[,(kx + sid[i])] <- 1 # put 1 on dummy school sid[i]
    
    XX <- XX + t(Xt) %*% Xt
    XMS <- XMS + t(Xt) %*% Mt %*% S[[i]]
    
  }
  estbeta <<- solve(XX) %*% XMS ## save estimated beta as global variable
  s2 <- 0
  ldetm <- 0
  nn <- 0
  for (i in 1:length(D)){
    nt <- nrow(X[[i]])
    Mt <- diag(nt)-lambda*P[[i]]
    Xt <- cbind(X[[i]],matrix(0,nt,ncol(schdummy))) # matrix of expl var, including school dummies
    Xt[,(kx + sid[i])] <- 1 # put 1 on dummy i
    
    eet <- t( Mt %*% S[[i]] - Xt %*% estbeta ) %*% ( Mt %*% S[[i]] - Xt %*% estbeta )
    s2 <- s2 + c(eet)
    
    ldetm <- ldetm + log(det(Mt)) 
    nn <- nn + nt
  }
  s2est <<- s2/nn # save estimated sigma2 as global variable
  likout <- -0.5*nn*log(s2/nn)+ldetm -0.5*nn*(log(2*pi)+1) # concentrated likelihood
  return(-likout ) # the objective function will be minimized
}

################
### JOINT LIKELIHOOD
################

jlik <- function(theta){
  ## the likelihood uses a two step procedure. Conditional on theta, it computes P(G|s) and the bounds for lambda (=phi in paper)
  ## then, it optimizes the likelihood P(s) at the maximized value for lambda.
  
  hld <- 99 # initialize bound on lambda
  for (i in 1:length(D)){
    
    ## build probability matrix
    nt <- nrow(X[[i]])
    Pt <- matrix(0,nt,nt)
    for (j in 1:length(D[[i]])){
      Pt <- Pt + theta[j]*D[[i]][[j]] # add contribution of expl var j
    }
    Pt <- pnorm(Pt + theta[(length(D[[i]])+sid[i])]) # add school dummy and compute proba
    diag(Pt) <- 0
    hld <- min(hld,(1/norm(Pt,"2"))) # update bound (tighten)
  }
  
  likt <- function(lbda) liksar(lbda,theta) # creates a function of lambda
  ls <- optim(0,likt,method="Brent",lower=(-hld+1e-8),upper=(hld-1e-8)) # optimize SAR likelihood
  lambdaest <<- ls$par ## save lambda as global variable
  return(-objproba(theta)+ls$value) # obj function will be minimized
}

jlikhes <- function(longtheta){
  
  ## joint likelihood, not concentrated. Is used for the computation of the numerical hessian
  
  thetaprob <- longtheta[1:nproba]
  thetabeta <- longtheta[(nproba+1):(nproba+nbeta)]
  thetal <- longtheta[(nproba+nbeta+1)]
  thetas <- longtheta[(nproba+nbeta+2)]
  
  likproba <- objproba(thetaprob) # gets P(G|s)
  
  ## for each school
  kx <- ncol(X[[1]])
  liksarout <- 0
  nn <- 0
  P <<- list("vector",length(D))
  Peq <<- list("vector",length(D))
  for (i in 1:length(D)){
    
    ## build probability matrix
    nt <- nrow(X[[i]])
    Pt <- matrix(0,nt,nt)
    for (j in 1:length(D[[i]])){
      Pt <- Pt + thetaprob[j]*D[[i]][[j]] # add contribution of expl var j
    }
    Pt <- pnorm(Pt + thetaprob[(length(D[[i]])+sid[i])])#/Nvec[i] # add school dummy and compute proba
    diag(Pt) <- 0
    P[[i]] <<- Pt
    ## temporary matrix to compute beta
    Mt <- diag(nt)-thetal*Pt
    Xt <- cbind(X[[i]],matrix(0,nt,ncol(schdummy))) # matrix of expl var, including school dummies
    Xt[,(kx + sid[i])] <- 1 # put 1 on dummy i
    ept <- Mt%*%S[[i]] - Xt%*%matrix(thetabeta,nbeta,1)
    #seq <- solve(Mt)%*%Xt%*%matrix(thetabeta,nbeta,1)
    #Peq[[i]] <<- Pt*(seq%*%t(seq))
    liksarout <- liksarout + log(det(Mt)) - c((0.5/thetas)*t(ept)%*%ept)
    nn <- nn + nt
  }
  liksarout <- liksarout - 0.5*nn*log(thetas) - 0.5*nn*log(2*pi)  # P(s)
  return(likproba + liksarout)
}


################
### JOINT LIKELIHOOD -- restricted to single school-grade
################

objproba_level <- function(theta,lev){
  ## conditional likelihood P(G|s)
  print(theta)
  flag_lev <- outdta$group==lev
  rho <- pnorm(as.matrix(outdta[flag_lev,2:12])%*%matrix(theta[1:11],11,1)+schdummy[flag_lev,]%*%matrix(theta[12:(12+ncol(schdummy)-1)],ncol(schdummy),1))
  proba <- outdta[flag_lev,"social"]*rho
  lik <- outdta[flag_lev,"g"]*log(proba) + (1-outdta[flag_lev,"g"])*log(1-proba)
  return(sum(lik,na.rm = T))
}

jlikhes_level <- function(longtheta,levnum){
  
  ## joint likelihood, not concentrated. Is used for the computation of the numerical hessian
  
  thetaprob <- longtheta[1:nproba]
  thetabeta <- longtheta[(nproba+1):(nproba+nbeta)]
  thetal <- longtheta[(nproba+nbeta+1)]
  thetas <- longtheta[(nproba+nbeta+2)]
  
  likproba <- objproba_level(thetaprob,grplist[levnum]) # gets P(G|s)
  
  ## for each school
  kx <- ncol(X[[1]])
  i <- levnum
  ## build probability matrix
  nt <- nrow(X[[i]])
  Pt <- matrix(0,nt,nt)
  for (j in 1:length(D[[i]])){
    Pt <- Pt + thetaprob[j]*D[[i]][[j]] # add contribution of expl var j
  }
  Pt <- pnorm(Pt + thetaprob[(length(D[[i]])+sid[i])]) # add school dummy and compute proba
  diag(Pt) <- 0
  
  ## temporary matrix to compute beta
  Mt <- diag(nt)-thetal*Pt
  Xt <- cbind(X[[i]],matrix(0,nt,ncol(schdummy))) # matrix of expl var, including school dummies
  Xt[,(kx + sid[i])] <- 1 # put 1 on dummy i
  ept <- Mt%*%S[[i]] - Xt%*%matrix(thetabeta,nbeta,1)
  liksarout <- log(det(Mt)) - c((0.5/thetas)*t(ept)%*%ept)
  liksarout <- liksarout - 0.5*nt*log(thetas) - 0.5*nt*log(2*pi)  # P(s)|lev
  return(likproba + liksarout)
}


################################################################
################################################################
#################### Parents' model ############################
################################################################
################################################################


################################################################
###################### Format data #############################
################################################################

builddta <- function(){
  ##  Format data for the parents' model
  
  ### keep only relevant variables and clean data
  dtaProba <- dta
  
  dtaProba <- dtaProba[dtaProba$h1gi20<=12,] # drop individuals with no grade levels
  dtaProba <- dtaProba[((dtaProba$scid != 1) & (dtaProba$scid != 115)),]
  
  ## recode/clean variables, missing values, individual characteristics
  dtaProba$female <- dtaProba$bio_sex-1 #female
  dtaProba$bio_sex <- NULL
  dtaProba$age <- 95 - (dtaProba$h1gi1y+(dtaProba$h1gi1m-1)/12) #age
  dtaProba$age <- pmin(dtaProba$age,18)
  dtaProba$age <- pmax(dtaProba$age,12)
  dtaProba$h1gi1m <- NULL
  dtaProba$h1gi1y <- NULL
  dtaProba$h1gi3 <- NULL
  dtaProba[dtaProba$h1gi4>1,"h1gi4"] <- 0
  colnames(dtaProba)[which(colnames(dtaProba)=="h1gi4")] <- "hisp"
  dtaProba[dtaProba$h1gi6a>1,"h1gi6a"] <- 0
  colnames(dtaProba)[which(colnames(dtaProba)=="h1gi6a")] <- "white"
  dtaProba[dtaProba$h1gi6b>1,"h1gi6b"] <- 0
  colnames(dtaProba)[which(colnames(dtaProba)=="h1gi6b")] <- "black"
  dtaProba$h1gi6c <- NULL
  dtaProba[dtaProba$h1gi6d>1,"h1gi6d"] <- 0
  colnames(dtaProba)[which(colnames(dtaProba)=="h1gi6d")] <- "asian"
  dtaProba$h1gi6e <- NULL
  dtaProba$h1rm4 <- NULL
  dtaProba[dtaProba$h1rm5>1,"h1rm5"] <- 0
  colnames(dtaProba)[which(colnames(dtaProba)=="h1rm5")] <- "momworks"
  dtaProba$type <- as.numeric((dtaProba$h1rf1>=8 & dtaProba$h1rf1<=9) | (dtaProba$h1rm1>=8 & dtaProba$h1rm1<=9 ) ) # type 1 if either parent is college educated
  dtaProba$h1rf1 <- NULL
  dtaProba$h1rm1 <- NULL
  
  
  ## socialization effort
  dtaProba[dtaProba$h1da4>3,"h1da4"] <- NA
  dtaProba[dtaProba$h1da5>3,"h1da5"] <- NA
  dtaProba[dtaProba$h1da6>3,"h1da6"] <- NA
  dtaProba[dtaProba$h1da7>3,"h1da7"] <- NA
  dtaProba[dtaProba$h1nb1>1,"h1nb1"] <- 0
  dtaProba[dtaProba$h1nb4>1,"h1nb4"] <- 0
  dtaProba$schoolact <- rowSums(dtaProba[,c(which(colnames(dtaProba)=="s44a1"):which(colnames(dtaProba)=="s44a33"))],na.rm=T) # number of activities
  
  dtaProba$schoolact <- (pmin(dtaProba$schoolact,5)+1)/6
  dtaProba[,c(which(colnames(dtaProba)=="s44a1"):which(colnames(dtaProba)=="s44"))] <- NULL
  dtaProba$dailyact <- rowMeans(dtaProba[,c(which(colnames(dtaProba)=="h1da4"):which(colnames(dtaProba)=="h1da7"))],na.rm=T) # average daily activities
  dtaProba[,c(which(colnames(dtaProba)=="h1da4"):which(colnames(dtaProba)=="h1da7"))] <- NULL
  dtaProba$neighact <- rowMeans(dtaProba[,c("h1nb1","h1nb4")],na.rm=T) # average neighbourhood participation
  dtaProba[,c(which(colnames(dtaProba)=="h1nb1"):which(colnames(dtaProba)=="h1nb7"))] <- NULL
  
  dtaProba$social <- rowMeans(dtaProba[,c("schoolact","dailyact","neighact")],na.rm=T) # socialization index
  dtaProba$social <- pmin(dtaProba$social,quantile(dtaProba$social,0.99))/quantile(dtaProba$social,0.99) # trim top 1% and normalize between 0 and 1
  dtaProba$social <- log(dtaProba$social+1)/max(log(dtaProba$social+1)) # log scale
  #dtaProba[,c("schoolact","dailyact","neighact")] <- NULL
  
  ## parent educational effort
  ## decisions
  dtaProba[dtaProba$h1wp1>1,"h1wp1"] <- 1 # missing values
  dtaProba[dtaProba$h1wp2>1,"h1wp2"] <- 1 # missing values
  dtaProba[dtaProba$h1wp3>1,"h1wp3"] <- 1 # missing values
  dtaProba[dtaProba$h1wp4>1,"h1wp4"] <- 1 # missing values
  dtaProba[dtaProba$h1wp5>1,"h1wp5"] <- 1 # missing values
  dtaProba[dtaProba$h1wp6>1,"h1wp6"] <- 1 # missing values
  dtaProba[dtaProba$h1wp7>1,"h1wp7"] <- 1 # missing values
  dtaProba$decision <- 1- (dtaProba$h1wp1 + dtaProba$h1wp2 + dtaProba$h1wp3 + dtaProba$h1wp4 + dtaProba$h1wp5 + dtaProba$h1wp6 + dtaProba$h1wp7)/7 # variables on how much students make their own decisions
  
  ## caring
  dtaProba[dtaProba$h1wp10>5,"h1wp10"] <- 0 # missing values
  dtaProba[dtaProba$h1wp14>5,"h1wp14"] <- 0 # missing values
  dtaProba$care <- (dtaProba$h1wp10 + dtaProba$h1wp14)/10 # variables on how much students care
  
  # activities
  dtaProba[dtaProba$h1wp17h>1,"h1wp17h"] <- 0 # missing values
  dtaProba[dtaProba$h1wp17i>1,"h1wp17i"] <- 0 # missing values
  dtaProba[dtaProba$h1wp17j>1,"h1wp17j"] <- 0 # missing values
  dtaProba[dtaProba$h1wp18h>1,"h1wp18h"] <- 0 # missing values
  dtaProba[dtaProba$h1wp18i>1,"h1wp18i"] <- 0 # missing values
  dtaProba[dtaProba$h1wp18j>1,"h1wp18j"] <- 0 # missing values
  dtaProba$activities <- (dtaProba$h1wp17h + dtaProba$h1wp17i + dtaProba$h1wp17j + dtaProba$h1wp18h + dtaProba$h1wp18i + dtaProba$h1wp18j)/6 # variables on parent-child activities related to school
  dtaProba$peffort <- (dtaProba$decision + dtaProba$care + dtaProba$activities)/3
  ### note on missing values in constructing peffort: if missing values are random, it is equivalent to classical measurement error on y so it does not affect the estimation
  
  return(dtaProba)
}


################################################################
####################### Simulate H #############################
################################################################

buildEhomophily <- function(gammatilde){
  
  ## simulate Eh^t conditional on gammatilde
  jlikhes(gammatilde) # computes probabilities: P is written as a global variable
  hsim <- rep(0,nrow(parentsdta))
  pos <- 1
  for (i in 1:length(D)){
    Pt <- P[[i]]*(S[[i]]%*%t(S[[i]]))
    nt <- nrow(Pt)
    types <- X[[i]][,8]
    ## matrix = 1 if same type and 0 otherwise
    sametype <- matrix(as.numeric(matrix(rep(types,nt),nt,nt)==t(matrix(rep(types,nt),nt,nt))),nt,nt)
    for (s in 1:nsim){
      Gt <- matrix(rbinom((nt*nt),1,Pt),nt,nt) # draw random network given P
      hsim[pos:(pos+nt-1)] <- hsim[pos:(pos+nt-1)] + (rowSums(Gt*sametype)/pmax(rowSums(Gt),1))/nsim
    }
    pos <- pos + nt
  }
  return(hsim)
}

################################################################
################################################################
################## Policy simulations ##########################
################################################################
################################################################

simulmore <- function(bplus,thetain2,thetainL,thetainH){
  
  # compute sample size
  siz <- 0
  for (i in 1:length(X)){
    siz <- siz + nrow(X[[i]])
  }
  
  #Empty database
  collect <- as.data.frame(matrix(0,siz,8))
  colnames(collect) <- c("id","s","h","tau","school","type","degree","peducated")
  
  pos <- 1 # position to fill-in
  for (i in 1:length(X)){
    Xt <- X[[i]] # indiv. var
    nt <- nrow(Xt) # size of group i
    Dt <- P[[i]] # predicted preference bias (D using the paper's notation)
    bpt <- bplus[[i]] # policy shift
    hmat <- matrix(NA,nt,nsim)
    taumat <- matrix(NA,nt,nsim)
    for (sim in 1:nsim){
      et <- matrix(rnorm(nt,0,sqrt(s2est)),nt,1) # draw errors 
      bt <- Xt%*%matrix(thetain2[1:ncol(Xt)],ncol(Xt),1) + et + matrix(thetain2[(ncol(Xt)+i)],nt,1) + bpt # compute bt
      Mt <- diag(nt)-lambdaest*Dt
      st <- solve(Mt)%*%bt # counterfactual equilibrium socialization effort
      collect[pos:(pos+nt-1),"s"] <- collect[pos:(pos+nt-1),"s"] + st/nsim
      Pt <- Dt * (st%*%t(st)) # counterfactual equilibrium proba
      Pt <- pmin(pmax(Pt,0),1) # ensure all proba in 0,1 (not binding)
      Gt <- matrix(rbinom((nt*nt),1,Pt),nt,nt) # draw network
      sametype <- matrix(as.numeric(  matrix(rep(Xt[,8],nt),nt,nt)==t(matrix(rep(Xt[,8],nt),nt,nt)) ),nt,nt)
      ht <- (rowSums(Gt*sametype)/pmax(1,rowSums(Gt))) # homophily index
      hmat[,sim] <- ht
      collect[pos:(pos+nt-1),"h"] <- collect[pos:(pos+nt-1),"h"] + ht/nsim # fraction of same-type links
      collect[pos:(pos+nt-1),"degree"] <- collect[pos:(pos+nt-1),"degree"] + rowSums(Gt)/nsim
    }
    for (sim in 1:nsim){
      ## compute parents' counter factual equilibrium effort
      ht <- collect[pos:(pos+nt-1),"h"] # expected h
      tauL <- thetainL[8]*(1-ht) + thetainL[8+i]*ht + rnorm(length(ht),0,sqrt(mean(PEout0$residuals^2)))
      tauH <- thetainH[8]*ht + thetainH[8+i]*(1-ht) + rnorm(length(ht),0,sqrt(mean(PEout1$residuals^2)))
      for (j in 1:7){
        tauL <- tauL + thetainL[j]*(ht*Xt[,j])
        tauH <- tauH + thetainH[j]*((1-ht)*Xt[,j])
      }
      tau <- rep(NA,nt)
      tau[c(Xt[,8]==0)] <- tauL[c(Xt[,8]==0)]
      tau[c(Xt[,8]==1)] <- tauH[c(Xt[,8]==1)]
      collect[pos:(pos+nt-1),"tau"] <- collect[pos:(pos+nt-1),"tau"] + tau/nsim
      taumat[,sim] <- tau
    }
    collect[pos:(pos+nt-1),"type"] <- Xt[,8]
    tmat <- matrix(rep(Xt[,8],nsim),nt,nsim)
    collect[pos:(pos+nt-1),"peducated"] <- rowMeans(taumat + (1-taumat)*( hmat*tmat + (1-hmat)*(1-tmat) ))
    collect[pos:(pos+nt-1),"id"] <- i
    collect[pos:(pos+nt-1),"school"] <- sid[i]
    
    pos <- pos + nt
  }
  return(collect)
}

################################################################
################################################################
######################## Model Fit #############################
################################################################
################################################################


buildrealh <- function(){
  
  hsim <- rep(0,nrow(parentsdta))
  pos <- 1
  for (i in 1:length(D)){
    Pt <- P[[i]]*(S[[i]]%*%t(S[[i]]))
    nt <- nrow(Pt)
    types <- X[[i]][,8]
    ## matrix = 1 if same type and 0 otherwise
    sametype <- matrix(as.numeric(matrix(rep(types,nt),nt,nt)==t(matrix(rep(types,nt),nt,nt))),nt,nt)
    Gt <- G[[i]] # get real network
    hsim[pos:(pos+nt-1)] <- (rowSums(Gt*sametype)/pmax(rowSums(Gt),1))
    pos <- pos + nt
  }
  return(hsim)
}

buildEhomophily_single <- function(gammatilde){
  
  ## simulate Eh^t conditional on gammatilde
  jlikhes(gammatilde) # computes probabilities: P is written as a global variable
  hsim <- rep(0,nrow(parentsdta))
  pos <- 1
  for (i in 1:length(D)){
    Pt <- P[[i]]*(S[[i]]%*%t(S[[i]]))
    nt <- nrow(Pt)
    types <- X[[i]][,8]
    ## matrix = 1 if same type and 0 otherwise
    sametype <- matrix(as.numeric(matrix(rep(types,nt),nt,nt)==t(matrix(rep(types,nt),nt,nt))),nt,nt)
    Gt <- matrix(rbinom((nt*nt),1,Pt),nt,nt) # draw random network given P
    hsim[pos:(pos+nt-1)] <- + (rowSums(Gt*sametype)/pmax(rowSums(Gt),1))
    pos <- pos + nt
  }
  return(hsim)
}

objprobit <- function(theta){
  ## conditional likelihood P(G|s=1)
  print(theta)
  rho <- pnorm(as.matrix(outdta[,2:12])%*%matrix(theta[1:11],11,1)+schdummy%*%matrix(theta[12:(12+ncol(schdummy)-1)],ncol(schdummy),1))
  proba <- rho
  lik <- outdta[,"g"]*log(proba) + (1-outdta[,"g"])*log(1-proba)
  return(sum(lik,na.rm = T))
}

gradprobit <- function(theta){
  ## gradient of the conditional likelihood P(G|s)
  rho <- as.matrix(outdta[,2:12])%*%matrix(theta[1:11],11,1)+schdummy%*%matrix(theta[12:(12+ncol(schdummy)-1)],ncol(schdummy),1)
  rhop <- dnorm(rho)
  rho <- pnorm(rho)
  proba <- rho
  mbeta <- ((outdta[,"g"]-proba)/(proba*(1-proba)))*rhop
  outgrad <- matrix(NA,N,length(theta))
  for (i in 1:11){
    outgrad[,i] <- mbeta*outdta[,(i+1)]
  }
  for (i in 1:ncol(schdummy)){
    outgrad[,(i+11)] <- mbeta*schdummy[,i]
  }
  return(colSums(outgrad,na.rm = T))
}


predprobit <- function(theta){
  rho <- pnorm(as.matrix(outdta[,2:12])%*%matrix(theta[1:11],11,1)+schdummy%*%matrix(theta[12:(12+ncol(schdummy)-1)],ncol(schdummy),1))
  proba <- rho
  outint <- matrix(NA,nrow(outdta),3)
  outint[,1] <- outdta[,"g"]
  outint[,2] <- proba
  pos <- 1
  for (i in 1:length(D)){
    Pt <- P[[i]]*(S[[i]]%*%t(S[[i]]))
    diag(Pt) <- NA
    lPt <- c(Pt)
    lPt <- lPt[is.na(lPt)==F]
    outint[pos:(pos+length(lPt)-1),3] <- lPt
    pos <- pos + length(lPt)
  }
  outint <- as.data.frame(outint)
  colnames(outint) <- c("data","probit","model")
  return(outint)
}



