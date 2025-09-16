# Install required packages
library(dplyr)
# Load required functions
source("/scratch/yya188/MCEM.R")

args <- commandArgs()
n <- as.numeric(args[6])
print(n)

K <- 100
path <- "/scratch/yya188/data"
#files2use <- list.files(path,pattern="RData")
setwd(path)

start_index <- 5*(n-1)+1
seq <- c(start_index,start_index+1,start_index+2,start_index+3,start_index+4)
files2use <- c(paste0("res_",seq,".RData"))
for (i in seq) {
  load(files2use[which(seq==i)])
  set.seed(1234*i+K)
  pp <- profilelkhd(data=res$data,mvals=c(1:10),N=100,weight=res$weight)
  f <- function(x) {stats:::predict.smooth.spline(ss,x)$y}
  ss <- smooth.spline(pp)
  m <- optimize(f,lower=1,upper=10,maximum=T)$maximum
  res$pp <- pp
  res$m <- m
  save(res,file=paste0("res_",i,".RData"))
}










