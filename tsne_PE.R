setwd("~/Documents/Patients/Cytokine Phase 1")

rm(list=ls())
library(SNFtool)
library(readr)
library(lsa)
library(softImpute)
library(Rtsne)
library(rgl)
library(dbscan)
library(bnstruct)

affinityMatrix <- function (diff, K = 20, sigma = 0.5) {
  N <- nrow(diff)
  diff <- (diff + t(diff))/2
  diag(diff) <- 0
  sortedColumns <- as.matrix(t(apply(diff, 2, sort)))
  finiteMean <- function(x) {
    return(mean(x[is.finite(x)]))
  }
  means <- apply(sortedColumns[, 1:K + 1], 1, finiteMean) + 
    .Machine$double.eps
  avg <- function(x, y) {
    return((x + y)/2)
  }
  Sig <- outer(means, means, avg)*2/3 + diff/3 + .Machine$double.eps
  Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
  densities <- dnorm(diff, 0, sigma*Sig, log = FALSE)
  # densities <- exp(-(diff)^2/(sigma*Sig))
  W <- (densities + t(densities))/2
  return(W)
}


K=20
C=4
alpha=0.5
T=25
less=numeric()
more=numeric()

col=c("black","red","blue","green","pink","purple","orange")

P1=c(4,6, 6,4,2,6,5,3,6,6,2,2,1,3,1,3, 5,2,1,6,6, 6,6, 6,1,6,1, 3,1)
P2=c(3,0,-1,-1,0,6,6,5,6,4,0,3,1,3,3,3,-1,1,0,6,6,-1,6,-1,1,6,0,-1,0)

for(i in 1:length(P1)){
  if(P1[i] <= 3){
    less=cbind(less, i)
  }
  else{
    more=cbind(more,i)
  }
}

Data=read.csv("Cy1Data1.csv")
ncol = ncol(Data)
Data = as.matrix(Data[,3:ncol])
Data = knn.impute(Data, k=10, cat.var = 1:ncol(Data), to.impute = 1:nrow(Data), using = 1:nrow(Data))
Data_more = as.matrix(Data[more,])
# fit = softImpute(Data,rank.max=28,lambda=0)
# Data = complete(Data,fit)
SN = standardNormalization(t(standardNormalization(Data)))
dist = (dist2(as.matrix(SN),as.matrix(SN)))^(1/2)
W = affinityMatrix(dist, K, alpha)
group_more = spectralClustering(W, C)

tsne2d=Rtsne(SN, dim=2, initial_dims=29, perplexity=6, theta=0, check_duplicates=FALSE, max_iter=5000)
plot(tsne2d$Y, xlab="", ylab="", t='n', xaxt='n', yaxt='n')
text(tsne2d$Y, labels=1:64, col=col[group_more], cex=1.5)
title('P1 cytokines PE', cex=1.5)

for(i in 1:C){
  assign(paste("more_PE",i,sep=""), which(group_more == i, arr.ind=T))
}

K=10
C_snf=2
more_PE=numeric()

for(j in 1:C){
  Data_new=Data[,get(paste("more_PE",j,sep=""))]
  SN=standardNormalization(Data_new)
  dist=(dist2(as.matrix(SN),as.matrix(SN)))^(1/2)
  W=affinityMatrix(dist, K, alpha)

  groupcy1F=spectralClustering(W, C_snf)

  cygroup1P1=numeric()
  cygroup2P1=numeric()

  group=groupcy1F
  for(g in 1:length(group)){
    if(group[g]==1){
      cygroup1P1=cbind(cygroup1P1,c(g,P1[g],P2[g]))}
    else if(group[g]==2){
      cygroup2P1=cbind(cygroup2P1,c(g,P1[g],P2[g]))}
  }

  t=t.test(cygroup1P1[2,], cygroup2P1[2,])
  print(t[c(3,5)])
  #  print(ncol(cygroup1P1))
  #  print(ncol(cygroup2P1))
  #  print(cygroup1P1[1,])
  #  print(cygroup2P1[1,])
  print(paste('# of cytokines in cluster', j, ':', length(get(paste("more_PE",j,sep="")))))
  if(t[3] < 0.05){
    more_PE=c(more_PE,get(paste("more_PE",j,sep="")))
  }
}

if(length(more_PE) != 0){
  Data_new=Data[,sort(c(more_PE))]
  SN=standardNormalization(Data_new)
  dist=(dist2(as.matrix(SN),as.matrix(SN)))^(1/2)
  W=affinityMatrix(dist, K, alpha)

  groupcy1F=spectralClustering(W, C_snf)

  cygroup1P1=numeric()
  cygroup2P1=numeric()

  group=groupcy1F
  for(g in 1:length(group)){
    if(group[g]==1){
      cygroup1P1=cbind(cygroup1P1,c(g,P1[g],P2[g]))}
    else if(group[g]==2){
      cygroup2P1=cbind(cygroup2P1,c(g,P1[g],P2[g]))}
  }

  print(paste('# of cytokines:', length(more_PE)))
  t=t.test(cygroup1P1[2,], cygroup2P1[2,])
  print(t[c(3,5)])

  CyP1_PE=colnames(Data)[more_PE]
  print(CyP1_PE)
}

K=20
C=4
print('------------------------------------------')
print('------------------------------------------')
print('portal cytokines')
Data=read.csv("Cy1Data2.csv")
ncol=ncol(Data)
Data=as.matrix(Data[,3:ncol])
Data = knn.impute(Data, k=10, cat.var = 1:ncol(Data), to.impute = 1:nrow(Data), using = 1:nrow(Data))
Data_more=as.matrix(Data[more,])
# fit=softImpute(Data, rank.max=28, lambda=0)
# Data = complete(Data, fit)
SN = standardNormalization(t(standardNormalization(Data)))
dist=(dist2(as.matrix((SN)),as.matrix((SN))))^(1/2)
W=affinityMatrix(dist, K, alpha)
group_more=spectralClustering(W, C)

tsne2d=Rtsne((SN), dim=2, initial_dims=29, perplexity=6, theta=0, check_duplicates=FALSE, max_iter=5000)
plot(tsne2d$Y, xlab="", ylab="", t='n', xaxt='n', yaxt='n')
text(tsne2d$Y, labels=1:64, col=col[group_more], cex=1.5)
title('P1 cytokines PO', cex=1.5)

for(i in 1:C){
  assign(paste("more_PO",i,sep=""), which(group_more == i, arr.ind=T))
}

K=10
more_PO=numeric()

for(j in 1:C){
  Data_new=Data[,get(paste("more_PO",j,sep=""))]
  SN=standardNormalization(Data_new)
  dist=(dist2(as.matrix(SN),as.matrix(SN)))^(1/2)
  W=affinityMatrix(dist, K, alpha)

  groupcy1F=spectralClustering(W, C_snf)

  cygroup1P1=numeric()
  cygroup2P1=numeric()

  group=groupcy1F
  for(g in 1:length(group)){
    if(group[g]==1){
      cygroup1P1=cbind(cygroup1P1,c(g,P1[g],P2[g]))}
    else if(group[g]==2){
      cygroup2P1=cbind(cygroup2P1,c(g,P1[g],P2[g]))}
  }

  t=t.test(cygroup1P1[2,], cygroup2P1[2,])
  print(t[c(3,5)])
  print(ncol(cygroup1P1))
  print(ncol(cygroup2P1))
  # print(cygroup1P1[1,])
  # print(cygroup2P1[1,])
  print(paste('# of cytokines in cluster', j, ':', length(get(paste("more_PO",j,sep="")))))
  if(t[3] < 0.05){
    more_PO=c(more_PO,get(paste("more_PO",j,sep="")))
  }
}
#print(more)

Data_new=Data[,sort(c(more_PO))]
SN=standardNormalization(Data_new)
dist=(dist2(as.matrix(SN),as.matrix(SN)))^(1/2)
W=affinityMatrix(dist, K, alpha)

groupcy1F=spectralClustering(W, C_snf)

cygroup1P1=numeric()
cygroup2P1=numeric()

group=groupcy1F
for(g in 1:length(group)){
  if(group[g]==1){
    cygroup1P1=cbind(cygroup1P1,c(g,P1[g],P2[g]))}
  else if(group[g]==2){
    cygroup2P1=cbind(cygroup2P1,c(g,P1[g],P2[g]))}
}

print(paste('# of cytokines:', length(more_PO)))
t=t.test(cygroup1P1[2,], cygroup2P1[2,])
print(t[c(3,5)])

CyP1_PO=colnames(Data)[more_PO]
print(CyP1_PO)