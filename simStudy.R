library(logistf)
library(matrixStats)
library(arm)

## ------------- Utility functions -------------
makeSummaryMat <- function(NREPS) {
  mat <- matrix(NA,nrow=NREPS,ncol=3)
  colnames(mat) <- c("betahat","cover","test.rej")
  return(mat)
}
fitSummary_wald <- function(ff,beta) {
  betahat = ff$coef[2]
  conf.int = cbind(betahat-1.96*summary(ff)$coefficients[2,2],
                   betahat+1.96*summary(ff)$coefficients[2,2])
  cover = is.in(beta,conf.int)
  test.rej <- !is.in(0,conf.int) # if 0 in CI, then don't reject
  return(c(betahat,cover,test.rej))
}
fitSummary_profile <- function(ff,beta) {
  betahat = ff$coef[2]
  conf.int = confint(ff,level=0.95)[paste0("X",k),]
  cover = is.in(beta,conf.int)
  test.rej <- !is.in(0,conf.int) # if 0 in CI, then don't reject
  return(c(betahat,cover,test.rej))
}
fitSummary_firth <- function(ff,beta) {
  betahat = ff$coef[2]
  conf.int = c(ff$ci.lower[2],ff$ci.upper[2])
    #confint(ff,level=0.95)[paste0("X",k),]
  cover = is.in(beta,conf.int)
  test.rej <- !is.in(0,conf.int) # if 0 in CI, then don't reject
  return(c(betahat,cover,test.rej))
}
# find_ci <- function(m_candidate) {
#   ci <- matrix(data=NA,nrow=length(m_candidate),ncol=2)
#   for (i in 1:length(m_candidate)) {
#     ff <- try({ logF(formula,res$data,m_candidate[i],logF_control) })
#     ci[i,] <- confint(ff,level=0.95)[paste0("X",k),]
#   }
#   return (cbind(min(ci[,1]),max(ci[,2])))
# }
is.in <- function(beta,intl) {
  if(any(is.na(intl))) { return(NA) }
  if(beta >= intl[1] & beta <= intl[2]) { return(1) }
  return(0)
}
Sign <- function(x) {
  if (x >= 0) {return(1)}
  if (x < 0) {return(-1)}
}
simSummary <- function(fitSummary,maf,beta) {
  bb <- fitSummary[,"betahat"]
  mean <- mean(bb)
  bias <- mean(Sign(beta)*(bb-beta),na.rm=TRUE)
  SD <- sd(bb,na.rm=TRUE)
  MSE <- bias^2 + SD^2
  cover <- mean(fitSummary[,"cover"],na.rm=TRUE)
  test.rej <- mean(fitSummary[,"test.rej"],na.rm=TRUE)
  return(c(maf,beta,mean,bias,SD,MSE,cover,test.rej))
}
addParams <- function(sum,resNames=c("MAF","Beta","Mean","Bias","SD","MSE","CP","Power")) {
  colnames(sum) <- resNames
  return(sum)
}


K <- 1000
NREPS <- 100
source('~/Desktop/SFU Ph.D./Project 1/New sim/MCEM_final_version.R')
#path <- "~/Desktop/SFU Ph.D./Project 1/New sim/data_500" 
path <- "~/Desktop/data"
logF_control <- glm.control(epsilon=1e-05,maxit=100,trace=FALSE)
Firth_control <- logistf::logistf.control(maxit=100,lconv=1e-05,gconv=1e-05,xconv=1e-05)
MLEMat <- logF_MCEMMat <- logF_LAMat <- FirthMat <- CauchyMat <- matrix(NA,nrow=K,ncol=8)


for (k in 1:K) {
  #print(k)
  #formula <- formula(paste0("case~","X",k,collapse="+"))
  formula <- formula(paste(paste0("case~","X",k,collapse="+")," + Population"))
  fitSumMLE <- fitSumlogF_MCEM <- fitSumlogF_LA <- fitSumFirth <- fitSumCauchy <- makeSummaryMat(NREPS)
  setwd(path)
  for (i in 1:NREPS) {
    load(paste0("res_",i,".RData")) # Load dataset
    #maf = res$MAF[k]
    #lambda = (15+1)*maf/(15*maf+1)
    lambda = (res$MAF[k]*(1-res$MAF[k]))^(-1/2)
    #lambda = 1

    ##-------------------------------------------------------
    ## MLE:
    ff = try({ glm(formula,res$data,family=binomial) })
    if(class(ff) != "try-error" && ff$iter < 100) {
      fitSumMLE[i,] <- fitSummary_wald(ff,res$beta[k])
    }
    ##-------------------------------------------------------

    ##-------------------------------------------------------
    ## LogF-MCEM:
    ff = try({ logF(formula,res$data,res$MCEM$m*lambda,logF_control) })
    if(class(ff) != "try-error" && ff$iter < 100) {
      fitSumlogF_MCEM[i,] <- fitSummary_profile(ff,res$beta[k])
    }
    ##-------------------------------------------------------
    
    ##-------------------------------------------------------
    ## LogF-LA:
    ff = try({ logF(formula,res$data,res$LA$m,logF_control) })
    if(class(ff) != "try-error" && ff$iter < 100) {
      fitSumlogF_LA[i,] <- fitSummary_profile(ff,res$beta[k])
    }
    ##-------------------------------------------------------

    ## -------------------------------------------------------
    ## Firth:
    ff = try({ logistf::logistf(formula,res$data,firth=T,pl=T,control=Firth_control) })
    if(class(ff) != "try-error" && ff$iter < 100) {
      fitSumFirth[i,] <- fitSummary_firth(ff,res$beta[k])
    }
    # ##-------------------------------------------------------
    # 
    # ##-------------------------------------------------------
    # ## Cauchy:
    ff = try({ bayesglm(formula,res$data,family=binomial) })
    if(class(ff) != "try-error" && ff$iter < 100) {
      fitSumCauchy[i,] <- fitSummary_wald(ff,res$beta[k])
    }
    ##-------------------------------------------------------
  }
  MLEMat[k,] <- simSummary(fitSumMLE,res$MAF[k],res$beta[k])
  logF_MCEMMat[k,] <- simSummary(fitSumlogF_MCEM ,res$MAF[k],res$beta[k])
  logF_LAMat[k,] <- simSummary(fitSumlogF_LA,res$MAF[k],res$beta[k])
  FirthMat[k,] <- simSummary(fitSumFirth,res$MAF[k],res$beta[k])
  CauchyMat[k,] <- simSummary(fitSumCauchy,res$MAF[k],res$beta[k])
}
MLEMat <- addParams(MLEMat)
logF_MCEMMat <- addParams(logF_MCEMMat)
logF_LAMat <- addParams(logF_LAMat)
FirthMat <- addParams(FirthMat)
CauchyMat <- addParams(CauchyMat)
sum <- list(MLE=MLEMat,logF_MCEM=logF_MCEMMat,logF_LA=logF_LAMat,Firth=FirthMat,Cauchy=CauchyMat)
setwd("/Users/daisyyu/Desktop")
save(sum,file="Summary.RData")


