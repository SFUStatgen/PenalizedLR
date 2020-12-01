datai<-1
estimate.m<-numeric(100)
fail.count<-0
for (datai in 1:100){
load(paste0("K=50/data_K=50_",datai,".RData"))
y<-sim_data$case
XM<-as.matrix(sim_data[,-c(ncol(sim_data),(ncol(sim_data)-1))])

n<-nrow(sim_data)
p<-ncol(XM)-1


optimLA<-function(alpha,m,X){
  ftracer=0
  tracer<-matrix(0, nrow=1,ncol=2)
  i=1
  for (i in 1:100){
  dlogPenalisedL<-function(beta){
    sum(X*y-(X*exp(alpha+beta*as.numeric(X))/(1+exp(alpha+beta*as.numeric(X)))))-m/2+m*exp(-beta)/(1+exp(-beta))
  }
  beta_max<-uniroot(dlogPenalisedL, c(-20,20))$root
  logLP_betamax<-function(alpha0){
    temp1<-sum(X^2*exp(alpha0+beta_max*X)/(1+exp(alpha0+beta_max*X)))-sum(X^2*(exp(alpha0+beta_max*X)/(1+exp(alpha0+beta_max*X)))^2)
    temp2<-exp(-beta_max)/(1+exp(-beta_max))-(exp(-beta_max)/(1+exp(-beta_max)))^2
    c=temp1+temp2
    LP<-sum(y*(alpha0+beta_max*as.numeric(X))-log(1+exp(alpha0+beta_max*as.numeric(X))))-log(beta(m/2,m/2))-m/2*beta_max-m*log(1+exp(-beta_max))+0.5*log(c)
    LP
  }
  opt.result<-optimize(f=logLP_betamax,lower=-10,upper=10,
                    maximum=TRUE, tol=0.0001)
  if(abs(opt.result$objective-ftracer)>=0.0001*abs(ftracer)){
    alpha=opt.result$maximum
    tracer<-rbind(tracer,c(alpha,opt.result$objective))
    ftracer<-opt.result$objective
  }else{
    break
  }
  }
  tracer
}

ms<-seq(1,10, by=0.1)
logl<-numeric(length(ms))
for (im in 1:length(ms)){
ini.m<-ms[im]
di=1
for (di in 1:(p+1)){
ini.alpha<-(-4)
tracer1<-optimLA(ini.alpha,ini.m,XM[,di])
logl[im]<-logl[im]+tracer1[nrow(tracer1),2]
}
}


if(ms[which.max(logl)]==10){
  print(paste0("dataset ",datai,": ",ms[which.max(logl)]))
  fail.count<-fail.count+1
  if(any(apply(XM,2,mean)==0)){print("all 0 variant exists")
    all0.count<-all0.count+1}
}
estimate.m[datai]<-ms[which.max(logl)]
#plot(ms,logl,type="b")
}

save(estimate.m,file=paste0("original_estimated_m_K=",p+1,".RData"))

mean(estimate.m)
sd(estimate.m)
mean(estimate.m)+1.96*sd(estimate.m)/10
mean(estimate.m)-1.96*sd(estimate.m)/10
