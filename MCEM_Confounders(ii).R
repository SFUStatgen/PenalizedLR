# library(dplyr)
# library(SNPknock)
# library(readxl)
# library(factoextra)
# library(FactoMineR)
# library(ggfortify)

simKnockoffGenotypes=function(n,Beta=log(rf(K,m,m))) {
  # Input:
  # - n is the sample size
  # - Beta is the value of log-OR simulated from log-F(m,m) distribution
  # Output: sim_data
  
  # Fit the Hidden Markov Model on Genotype data with fastPHASE
  sampleID=rownames(X)
  Xinp_file=writeXtoInp(X)
  fp_path="/scratch/yya188/fastPHASE"
  #fp_path="/Users/daisyyu/Desktop/Simulation/fastPHASE"
  fp_outPath=runFastPhase(fp_path, Xinp_file, K=12, numit=30)
  r_file=paste(fp_outPath,"_rhat.txt",sep="")
  alpha_file=paste(fp_outPath,"_alphahat.txt",sep="")
  theta_file=paste(fp_outPath,"_thetahat.txt",sep="")
  char_file=paste(fp_outPath,"_origchars",sep="")
  hmm=loadHMM(r_file,alpha_file,theta_file,char_file,compact=T,phase=T)
  
  # Simulate knockoff genotypes
  seed=sample(10000:50000,20,replace=F)
  data=as.data.frame(NULL)
  for (i in 1:20) {
    Xk=knockoffGenotypes(X,hmm$r,hmm$alpha,hmm$theta,seed=seed[i])
    data=rbind(data,Xk)
  }
  
  # Calculate PCs
  pc=prcomp(data)
  # fviz_eig(pc)
  # sample_info=read_excel("/Users/daisyyu/Desktop/Simulation/1000GenomesProject/20130606_sample_info.xlsx")
  # data$Sample=rep(sampleID,20)
  # data=data %>% left_join(sample_info, by="Sample") %>% mutate(Population=ifelse(Population %in% c("CEU","FIN","GBR","TSI","IBS"), "EUR", "EAS")) 
  # pca.plot=autoplot(pc,data=data,colour='Population')
  # pca.plot
  
  # Calculate disease probability for each individual
  # choose a subset of SNPs to be casually associated with the disease
  cov_ind=sample(1:84,K,replace=F)
  data=cbind(data[,cov_ind],pc$x[,1:10]) 
  data$case=NULL
  
  for (i in 1:dim(data)[1]) {
    s=sum(data[i,1:K]*Beta)
    prob=exp(s)/(1+exp(s))
    data$case[i]=rbinom(1,1,prob)
  }
  
  # Sample case/control status
  # case-control ratio: 1:1
  no_case=n/2
  no_con=n/2
  data_con=data %>% filter(case==0)
  data_case=data %>% filter(case==1)
  sim_data=rbind(data_con[sample(nrow(data_con),no_con,replace=F),],
                 data_case[sample(nrow(data_case),no_case,replace=F),])
  colnames(sim_data)[1:K]=c(paste0("X",1:K));rownames(sim_data)=NULL
  return (sim_data)
}


simUnmatched=function(n,Beta,p,scale=FALSE) {
  # Input:
  # - n is total sample size
  # - beta is value of parameter of interest
  # - p is number of nuisance covariates
  # Output: sim_data
  
  ncase=n/2; ncon=n/2  # assuming 1:1 con:case ratio
  Beta=c(Beta,rep(1,p))
  ncov=p+length(Beta)
  conX=caseX = NULL
  for(i in 1:ncov) {
    conX=cbind(conX,rnorm(ncon,mean=0,sd=1))
    caseX=cbind(caseX,rnorm(ncase,mean=Beta[i],sd=1))
  }
  X=rbind(caseX,conX)
  if(scale) X = round(scale(X))
  colnames(X)=paste0("x",1:ncov);rownames(X) = NULL
  case=c(rep(1,ncase),rep(0,ncon))
  return (data.frame(X,case))
}


MCEM=function(m,data,N) {
  # Input:
  # - m is the value of m
  # - data is the simulated case-control data obtained by simUnmatched()
  # - N is number of Monte Carlo replicates
  # Output: params
  
  model_PC_only=glm(case~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                    data=data,family=binomial(link="logit"),maxit=100)
  params=model_PC_only$coefficients[2:11]
  model_full=glm(case~.,data=data,family=binomial(link="logit"),maxit=100)
  AlphaStar=numeric()
  AlphaStar[1]=0;AlphaStar[2]=model_full$coefficients[1]

  p=2
  threshold=1E-04
  
  Weight=function(beta) {
    data_m=as.matrix(data)
    s=AlphaStar[p]+data_m[,3:12]%*%params+data_m[,2]*beta
    prod(exp(data_m[,1]*s)/(1+exp(s)))
  }
  
  Y=rep(data$case,times=N)
  betas=log(rf(N,m,m))
  
  O=numeric() # offset
  for (j in 1:N) {
    data_m=as.matrix(data)
    O=c(O,data_m[,3:12]%*%params+data_m[,2]*betas[j])
  }
  
  while(abs(AlphaStar[p]-AlphaStar[p-1])>=threshold) {
    W_t=numeric() # weight
    for (j in 1:N) {
      W_t[j]=Weight(betas[j])
    }
    W=rep(W_t,each=dim(data)[1])
    
    g=glm(Y~offset(O),weights=W,family=binomial(link="logit"),maxit=100)
    cat("EM iteration",p-1,":",g$coefficients,"\n")
    p=p+1
    AlphaStar[p]=g$coefficients
  }
  return(c(last(AlphaStar),params))
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
    data_m=as.matrix(data)
    s=params_k[1]+data_m[,c(3:12)]%*%params_k[2:11]+data_m[,2]*betas[j]
    lvec[j]=prod(exp(data_m[,1]*s)/(1+exp(s)))
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
      data_k=data %>% select(c(case,paste0("X",k),PC1:PC10)) %>% rename(X=paste0("X",k))
      #names(data_k)=c("case","X",paste0("PC",1:10))
      params_k=MCEM(m,data_k,N)
      ll[m]=ll[m]+lkhdk(params_k,data_k,m,N)
    }
  }
  return(ll)
}

save(simKnockoffGenotypes,simUnmatched,MCEM,lkhdk,profilelkhd,file="MCEM_Single-SNP.RData")

