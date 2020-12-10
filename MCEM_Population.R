# library(dplyr)
# library(SNPknock)
# library(readxl)
# library(factoextra)
# library(FactoMineR)
# library(ggfortify)

simSNPCov_popl=function(n,Beta=log(rf(K,m,m)),MAF1,MAF2,beta_s){
  ncase <- ncon <- n/2
  conX <- matrix(NA,ncol=K+1,nrow=ncon)
  caseX <- matrix(NA,ncol=K+1,nrow=ncase)
  
  # Simulate population status first (population is the K+1 column)
  f_0 <- f_1 <- 0.5
  conX[,K+1] <- rbinom(ncon,1,f_0)
  caseX[,K+1] <- sample(0:1,size=ncase,replace=TRUE,prob=c(f_0,f_1*exp(beta_s)))
  
  # Simulate SNVs conditional on population status
  for (i in 1:K){ 
    maf0 <- MAF1[i]
    maf1 <- MAF2[i]
    beta <- Beta[i]
    
    ## control
    con_pop0 <- which(conX[,K+1]==0)
    con_pop1 <- which(conX[,K+1]==1)
    
    conX[con_pop0,i] <- rbinom(length(con_pop0),size=2,prob=maf0)
    conX[con_pop1,i] <- rbinom(length(con_pop1),size=2,prob=maf1)
    
    ## case
    case_pop0 <- which(caseX[,K+1]==0)
    case_pop1 <- which(caseX[,K+1]==1)
    
    pp0 <- c((1-maf0)^2,2*maf0*(1-maf0)*exp(beta),maf0^2*exp(2*beta))
    pp1 <- c((1-maf1)^2,2*maf1*(1-maf1)*exp(beta),maf1^2*exp(2*beta))
    caseX[case_pop0,i] <- sample(0:2,size=length(case_pop0),replace=TRUE,prob=pp0)
    caseX[case_pop1,i] <- sample(0:2,size=length(case_pop1),replace=TRUE,prob=pp1)
  }
  X <- rbind(caseX,conX)
  colnames(X) <- c(paste0("X",1:K),"Population");rownames(X) <- NULL
  case <- c(rep(1,ncase),rep(0,ncon))
  return (data.frame(X,case))
}


MCEM=function(m,data,N) {
  # Input:
  # - m is the value of m
  # - data is the simulated case-control data obtained by simUnmatched()
  # - N is number of Monte Carlo replicates
  # Output: params
  
  model=glm(case~.,data=data,family=binomial(link="logit"),maxit=100)
  initial_params=model$coefficients[c(1,3)]
  params=matrix(0,ncol=2)
  params=rbind(params,initial_params)
  
  p=2
  threshold=1E-04
  
  Weight=function(beta) {
    #data_m=as.matrix(data)
    s=params[p,1]+data$Population*params[p,2]+data$X*beta
    prod(exp(data$case*s)/(1+exp(s)))
  }
  
  Y=rep(data$case,times=N)
  Pop=rep(data$Population,times=N)
  betas=log(rf(N,m,m))
  
  O=numeric() # offset
  for (j in 1:N) {
    O=c(O,data$X*betas[j])
  }
  
  while(norm(as.matrix(params[p,]-params[p-1,]))>=threshold) {
    W_t=numeric() # weight
    for (j in 1:N) {
      W_t[j]=Weight(betas[j])
    }
    W=rep(W_t,each=dim(data)[1])
    
    g=glm(Y~offset(O)+Pop,
          weights=W,family=binomial(link="logit"),maxit=100)
    cat("EM iteration",p-1,":",g$coefficients,"\n")
    p=p+1
    params=rbind(params,g$coefficients)
  }
  return(params[p,])
}


lkhdk=function(params_k,data,m,N) {
  # Input:
  # - alpha_k is the output of MCEM
  # - data is the simulated case-control data obtained by simUnmatched()
  # - m is the value of m
  # - N is number of Monte Carlo replicates
  # Output: Monte Carlo estimate of the profile likelihood
  
  betas=log(rf(N,m,m))
  lvec=rep(NA,N)
  for(j in 1:N) {
    #data_m=as.matrix(data)
    s=params_k[1]+data$Population*params_k[2]+data$X*betas[j]
    lvec[j]=prod(exp(data$case*s)/(1+exp(s)))
  }
  return(log(mean(lvec)))
}


profilelkhd=function(data,mvals,N) {
  # Input:  
  # - data is the simulated case-control data obtained by simUnmatched()
  # mvals is a set of values of m
  # - N is number of Monte Carlo replicates
  # Output: profile likelihood of m
  
  ll=rep(NA,length(mvals))
  for(m in mvals) {
    cat("Estimating profile log-likelihood for m =",m,"\n")
    ll[m]=0
    for(k in 1:K) {
      data_k=data %>% select(c(case,paste0("X",k),Population)) %>% rename(X=paste0("X",k))
      #names(data_k)=c("case","X",paste0("PC",1:10))
      params_k=MCEM(m,data_k,N) 
      ll[m]=ll[m]+lkhdk(params_k,data_k,m,N)
    }
  }
  return(ll)
}


save(simSNPCov_popl,MCEM,lkhdk,profilelkhd,file="MCEM_Population.RData")

