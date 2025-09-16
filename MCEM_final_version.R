MCEM <- function(m,data_k,C_index,N) {
  # Input:
  # - m is the value of m
  # - data_k is the data consisting of the phenotype, SNV_k and confounders
  # - C_index is the index of confounders in data
  # - N is the number of Monte Carlo replicates
  # Output: params_k
  
  model <- glm(Phenotype~.,data=data_k,family=binomial(link="logit"),maxit=100)
  initial_params <- model$coefficients[-2] # remove coefficient of X
  c <- length(C_index) # number of confounders
  params <- matrix(0,ncol=c+1)
  params <- rbind(params,initial_params)
  
  p <- 2
  threshold <- 1E-04
  
  Weight <- function(beta) {
    s <- params[p,1]+data_k$X*beta+as.matrix(data_k[,-c(1,2)])%*%(params[p,-1])
    lkhd <- exp(data_k$Phenotype*s)/(1+exp(s))
    sum(log(lkhd))
  }
  
  Y <- rep(data_k$Phenotype,times=N)
  Cov <- matrix(data=NA,ncol=c,nrow=N*dim(data_k)[1])
  if (c != 0) {
    for (i in 1:c) {
      Cov[,i] <- rep(data_k[,2+i])
    }
  }
  betas <- log(rf(N,m,m))
  
  O <- numeric() # offset
  for (j in 1:N) {
    O <- c(O,data_k$X*betas[j])
  }
  
  while(norm(as.matrix(params[p,]-params[p-1,]))>=threshold) {
    W_t <- numeric() # weight
    for (j in 1:N) {
      W_t[j] <- Weight(betas[j])
    }
    W_t <- exp(W_t-max(W_t))
    W <- rep(W_t,each=dim(data_k)[1])
    
    if (c != 0) {g <- glm(Y~offset(O)+Cov,weights=W,family="binomial",maxit=100)}
    if (c == 0) {g <- glm(Y~offset(O),weights=W,family=binomial(link="logit"),maxit=100)}
    cat("EM iteration",p-1,":",g$coefficients,"\n")
    p <- p+1
    params <- rbind(params,g$coefficients)
  }
  return(params[p,])
}


lkhdk <- function(m,data_k,params_k,N) {
  # Input:
  # - m is the value of m
  # - data_k is the data consisting of the phenotype, SNV_k and confounders
  # - params_k is the output of MCEM()
  # - N is the number of Monte Carlo replicates
  # Output: Monte Carlo estimate of the profile likelihood
  
  betas <- log(rf(N,m,m))
  lvec <- rep(NA,N)
  for(j in 1:N) {
    s <- params_k[1]+data_k$X*betas[j]+as.matrix(data_k[,-c(1,2)])%*%(params_k[-1])
    lvec[j] <- prod(exp(data_k$Phenotype*s)/(1+exp(s)))
  }
  return(log(mean(lvec)))
}


profilelkhd <- function(m,data,weight,Y_index,X_index,C_index,N) {
  # Input: 
  # - m is the value of m
  # - data is the GWAS data
  # - weight is a list of weights for SNVs
  # - Y_index is the index of phenotype in data
  # - X_index is the index of SNVs in data
  # - C_index is the index of confounders in data
  # - N is the number of Monte Carlo replicates
  # Output: profile marginal likelihood of m
  
  ll <- 0
  K <- length(X_index)
  for(k in 1:K) {
    cat("k =",k,"\n")
    #if (length(C_index) == 0) {data_k <- data[,c(Y_index,X_index[k])]}
    #else {data_k <- data[,c(Y_index,X_index[k],C_index)]}
    data_k <- data[,c(Y_index,X_index[k],C_index)]
    colnames(data_k)[1:2]=c("Phenotype","X")
    params_k <- MCEM(m,data_k,C_index,N) 
    ll <- ll+weight[k]*lkhdk(m,data_k,params_k,N)
  }
  return(ll)
} 


optimLA <- function(alpha,m,X,y,beta_int=NULL,confounding_factor=NULL){
  ftracer=0
  tracer<-matrix(0, nrow=1,ncol=3)
  i=1
  if (is.null(confounding_factor)){
    for (i in 1:50){
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
        tracer<-rbind(tracer,c(alpha,NA,opt.result$objective))
        ftracer<-opt.result$objective
      }else{
        break
      }
    }
  }else{
    for (i in 1:50){
      dlogPenalisedL<-function(beta){
        sum(X*y-(X*exp(alpha+confounding_factor*beta_int+beta*as.numeric(X))/(1+exp(alpha+beta_int*confounding_factor+beta*as.numeric(X)))))-m/2+m*exp(-beta)/(1+exp(-beta))
      }
      beta_max<-uniroot(dlogPenalisedL, c(-20,20))$root
      logLP_betamax<-function(alpha0beta0){
        alpha0<-alpha0beta0[1]
        beta0<-alpha0beta0[2]
        temp1<-sum(X^2*exp(alpha0+confounding_factor*beta0+beta_max*X)/(1+exp(alpha0+confounding_factor*beta0+beta_max*X)))-sum(X^2*(exp(alpha0+confounding_factor*beta0+beta_max*X)/(1+exp(alpha0+confounding_factor*beta0+beta_max*X)))^2)
        temp2<-exp(-beta_max)/(1+exp(-beta_max))-(exp(-beta_max)/(1+exp(-beta_max)))^2
        c=temp1+temp2
        LP<-sum(y*(alpha0+confounding_factor*beta0+beta_max*as.numeric(X))-log(1+exp(alpha0+confounding_factor*beta0+beta_max*as.numeric(X))))-log(beta(m/2,m/2))-m/2*beta_max-m*log(1+exp(-beta_max))+0.5*log(c)
        LP
      }
      opt.result<-optim(par=c(alpha,beta_int),
                        fn=logLP_betamax,
                        method="L-BFGS-B",
                        control=list(fnscale=-1))
      if(abs(opt.result$value-ftracer)>=0.0001*abs(ftracer)){
        alpha=opt.result$par[1]
        beta_int=opt.result$par[2]
        tracer<-rbind(tracer,c(alpha,beta_int,opt.result$value))
        ftracer<-opt.result$value
      }else{
        break
      }
    }
  }
  tracer
}


LAapproxL <- function(m,XM,y,ini_alpha,confounding_factor,weight_col){
  n<-nrow(XM)
  p<-ncol(XM)-1
  logl<-0
  if(is.null(ini_alpha)){ini_alpha=-4}
  ini.m<-m
  di=1
  for (di in 1:(p+1)){
    if(is.null(confounding_factor)){
      tracer1<-optimLA(ini_alpha,ini.m,XM[,di],y)}else{
        tracer1<-optimLA(ini_alpha,ini.m,XM[,di],y,0,confounding_factor)
      }
    if(is.null(weight_col)){
      logl<-logl+tracer1[nrow(tracer1),3]}else{
        logl<-logl+weight_col[di]*tracer1[nrow(tracer1),3]
      }
  }
  logl
}


get_m <- function(mvals,data,weight=NULL,Y_index,X_index,C_index=NULL,N,method="MCEM",ini_alpha=NULL){
  require(parallel)
  # Input:
  # - mvals is a list of m values
  # - data is the GWAS data
  # - weight is a list of weights for SNVs
  # - Y_index is the index of phenotype in data
  # - X_index is the index of SNVs in data
  # - C_index is the index of confounders in data
  # - N is the number of Monte Carlo replicates
  # - method is either "MCEM" or "LA"
  # Output: an estimate of m
  if (method == "MCEM") {
    pp <- unlist(mclapply(X=mvals,FUN=profilelkhd,data=data,weight=weight,
                          Y_index=Y_index,X_index=X_index,C_index=C_index,N=N,
                          mc.preschedule=T,mc.cores=4))
    f <- function(x) {stats:::predict.smooth.spline(ss,x)$y}
    ss <- smooth.spline(pp)
    mhat <- optimize(f,lower=mvals[1],upper=mvals[length(mvals)],maximum=T)$maximum
    return(list(lkhd=pp,m=mhat))
  }
  
  if (method == "LA"){
    logl_allm <- unlist(mclapply(X=mvals,FUN=LAapproxL,XM=as.matrix(data[,X_index]),
                                 y=data[,Y_index],ini_alpha=ini_alpha,
                                 confounding_factor=C_index,weight_col=weight,mc.cores=4))
    f <- function(x) {stats:::predict.smooth.spline(ss,x)$y}
    ss <- smooth.spline(logl_allm)
    mhat <- optimize(f,lower=mvals[1],upper=mvals[length(mvals)],maximum=T)$maximum
    return(list(lkhd=logl_allm,m=mhat))
  }
}


logF <- function(form,data,m,control=glm.control()) {
  # Input:
  # - form is an R formula
  # - data is the GWAS data
  # - m is the degree of freedom for the log-F prior
  # - glm.control is the algorithm control arguments to be passed to glm()
  # Output:
  
  #---------------
  # Step 1: Extract (i) the response and (ii) the design matrix
  # from the input formula and data frame so that we can augment them.
  mf <- model.frame(form,data)
  D <- model.response(mf)
  X <- model.matrix(form,data)
  xvars <- colnames(X)
  # Step 2 (augmentation): one pseudo-observation for each covariate,
  #  where the response is m/2 successes and m/2 failures (even if
  #  m is an odd number) and the covariates are all zeros except for
  #  a one indicating the index of the covariate.
  #  Following the recommendation of Greenland and Mansournia (2015; p. 3139)
  #  we do not penalize the intercept.
  n <- rep(1,length(D))
  zeros <- rep(0,ncol(X))
  pseudoD <- m/2; pseudoN <- m
  # for(i in 2:ncol(X)) {
  #   D <- c(D,pseudoD); n <- c(n,pseudoN)
  #   pseudoX <- zeros; pseudoX[i]  <- 1; X <- rbind(X,pseudoX)
  # }
  D <- c(D,pseudoD); n <- c(n,pseudoN)
  pseudoX <- zeros; pseudoX[2]  <- 1; X <- rbind(X,pseudoX) # assume the second variable is SNV
  # Step 3: Set up a response matrix with columns for number of successes
  # and number of failures.
  Y <- cbind(D,n-D)
  # Step 4: Set up X's as a data.frame with null rownames and correct colnames.
  rownames(X) <- NULL
  X <- data.frame(X)
  # Seems that formulas ignore brackets, so  (Intercept) will look
  # for the variable Intercept. Remove brackets.
  xvars <- c("Intercept",xvars[-1])
  names(X) <- xvars
  # Step 5: set up a formula and call glm()
  form <- formula(paste("Y~ -1 + ",paste0(xvars,collapse="+")))
  out <- glm(form,data=X,family=binomial(),control=control)
  return(out)
}

logF_lkhd <- function(form,form_null,data,ff,m){
  # Input:
  # - form is an R formula
  # - form_null is an R formula without the covariate of interest
  # - data is the data set
  # - ff is the fitted LogF model generated by logF()
  # - m is the value of m
  # Output: LRT p-value
  
  # Extract the response and the design matrix from the input formula and data frame
  mf <- model.frame(form,data)
  Y <- model.response(mf)
  X <- model.matrix(form,data)
  
  # Lkhd for null model
  g <- glm(form_null,data,family="binomial")
  coef <- g$coefficients
  dem <- 1+exp(X[,-2]%*%coef)
  num <- exp(Y*(X[,-2]%*%coef))
  beta <- 0 # under null: beta=0
  penalty <- exp(beta*m/2)/(1+exp(beta))^m*1/beta(m/2,m/2)
  lkhd_null <- log(prod(num/dem)*penalty)
  
  # Lkhd for full model
  coef <- ff$coefficients
  dem <- 1+exp(X%*%coef)
  num <- exp(Y*(X%*%coef))
  beta <- as.numeric(coef[2]) # under null: beta=0
  penalty <- exp(beta*m/2)/(1+exp(beta))^m*1/beta(m/2,m/2)
  lkhd_full <- log(prod(num/dem)*penalty)
  return (pchisq(-2*(lkhd_null-lkhd_full),df=1,lower.tail = F))
}


