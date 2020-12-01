datai<-1
estimate.m<-numeric(100)
fail.count<-0
for (datai in 1:100){
  load(paste0("K=50/data_K=50_",datai,".RData"))
  y<-sim_data$case
  XM<-as.matrix(sim_data[,-c(ncol(sim_data),(ncol(sim_data)-1))])
  pop<-sim_data$Population
  n<-nrow(sim_data)
  p<-ncol(XM)-1
  
  
  optimLA<-function(alpha,m,X,beta_int){
    ftracer=0
    tracer<-matrix(0, nrow=1,ncol=2)
    i=1
    for (i in 1:50){
      dlogPenalisedL<-function(beta){
        sum(X*y-(X*exp(alpha+(beta+pop*beta_int)*as.numeric(X))/(1+exp(alpha+(beta+beta_int*pop)*as.numeric(X)))))-m/2+m*exp(-beta)/(1+exp(-beta))
      }
      beta_max<-uniroot(dlogPenalisedL, c(-20,20))$root
      logLP_betamax<-function(alpha0beta0){
        alpha0<-alpha0beta0[1]
        beta0<-alpha0beta0[2]
        temp1<-sum(X^2*exp(alpha0+(beta_max+pop*beta0)*X)/(1+exp(alpha0+(beta_max+pop*beta0)*X)))-sum(X^2*(exp(alpha0+(beta_max+pop*beta0)*X)/(1+exp(alpha0+(beta_max+pop*beta0)*X)))^2)
        temp2<-exp(-beta_max)/(1+exp(-beta_max))-(exp(-beta_max)/(1+exp(-beta_max)))^2
        c=temp1+temp2
        LP<-sum(y*(alpha0+(beta_max+pop*beta0)*as.numeric(X))-log(1+exp(alpha0+(beta_max+pop*beta0)*as.numeric(X))))-log(beta(m/2,m/2))-m/2*beta_max-m*log(1+exp(-beta_max))+0.5*log(c)
        LP
      }
      opt.result<-optim(par=c(alpha,beta_int),
                        fn=logLP_betamax,
                        method="Nelder-Mead",
                        control=list(fnscale=-1, reltol=0.001))
      if(abs(opt.result$value-ftracer)>=0.001*abs(ftracer)){
        alpha=opt.result$par[1]
        beta_int=opt.result$par[2]
        tracer<-rbind(tracer,c(alpha,opt.result$objective))
        ftracer<-opt.result$value
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
      tracer1<-optimLA(ini.alpha,ini.m,XM[,di],1)
      logl[im]<-logl[im]+tracer1[nrow(tracer1),2]
    }
  }
  
  print(paste0("dataset ", datai, ": done"))
  
  if(ms[which.max(logl)]==10){
    print(paste0("dataset ",datai,": ",ms[which.max(logl)]))
    fail.count<-fail.count+1
    if(any(apply(XM,2,mean)==0)){print("all 0 variant exists")
      all0.count<-all0.count+1}
  }
  estimate.m[datai]<-ms[which.max(logl)]
  #plot(ms,logl,type="b")
}

save(estimate.m,file=paste0("estimated_m_K=",p+1,".RData"))

mean(estimate.m)
sd(estimate.m)
mean(estimate.m)+1.96*sd(estimate.m)/10
mean(estimate.m)-1.96*sd(estimate.m)/10
