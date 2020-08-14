library(dplyr)
library(SNPknock)

load("Input_dataset_1000GP.RData")
load("MCEM_Single-SNP.RData")

args <- commandArgs()
n <- as.numeric(args[6])
print(n)

m=4
K=10
set.seed(2020*n+K)
#sim_data=simUnmatched(n=200,Beta=log(rf(K,m,m)),p=0)     # simulate continuous covariates
sim_data=simKnockoffGenotypes(n=200,Beta=log(rf(K,m,m)))  # simulate SNP covariates
pp=profilelkhd(data=sim_data,mvals=c(1:20),N=1000)

save(pp,file=paste0("K=10_rep",n,".RData"))
save(sim_data,file=paste0("data_K=10_",n,".RData"))










