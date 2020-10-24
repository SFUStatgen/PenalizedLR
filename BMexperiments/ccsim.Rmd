# Case control simulation

ncase <- ncon <- 100
K <- 50
m <- 4
N <- 10000
beta0 <- (-4)
MAF <- 0.5

ccsim <- function(MAF) {
  X <- cbind(1,matrix(rbinom(N*K,size=2,p=MAF),ncol=K))
  beta <- c(beta0,log(rf(K,m,m)))
  linpred <- X %*% beta
  p <- exp(linpred)/(1+exp(linpred))
  cc <- rbinom(N,size=1,prob=p)
  caseind <- sample((1:N)[cc==1],size=ncase,replace=FALSE)
  conind <- sample((1:N)[cc==0],size=ncon,replace=FALSE)
  list(cc=c(rep(1,ncase),rep(0,ncase)),X = X[c(caseind,conind),],
       beta = beta)
}

dat1 <- ccsim(MAF=.5)


