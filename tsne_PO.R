## Phase 1 portal cytokines
## separate data into two groups of patients, based on their severity (ishak scores)
## do tsne on more severe patients cytokines, to see which cytokines cluster together
## do spectral clustering on severe patients' cytokines, to see whether the clustering agrees with that of tsne
## using the clusteres found by SC, do SC on patients, and see if any clusters can distinguish less severe from more severe


setwd("~/Documents/Patients/Cytokine Phase 1")

#rm(list=ls())
library(SNFtool)
library(readr)
library(lsa)
library(softImpute)
library(Rtsne)
library(rgl)
library(dbscan)

K=20
C=3
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

Data=read.csv("Cy1Data2.csv")
ncol=ncol(Data)
Data=as.matrix(Data[,3:ncol])
fit=softImpute(Data, rank.max=28, lambda=0)
Data = complete(Data, fit)
SN = standardNormalization(Data)
dist=dist2(as.matrix(t(SN)),as.matrix(t(SN)))
W=affinityMatrix(dist, K, alpha)
group_more=spectralClustering(W, C)

tsne2d=Rtsne(t(SN), dim=2, initial_dims=29, perplexity=6, theta=0, check_duplicates=FALSE, max_iter=5000)
plot(tsne2d$Y, xlab="", ylab="", t='n', xaxt='n', yaxt='n')
text(tsne2d$Y, labels=1:64, col=col[group_more], cex=1.5)
title('P1 cytokines PO', cex=1.5)

for(i in 1:C){
  assign(paste("more",i,sep=""), which(group_more == i, arr.ind=T))
}

K=10
C_snf=2
more_PO=numeric()

for(j in 1:C){
  Data_new=Data[,get(paste("more",j,sep=""))]
  SN=standardNormalization(Data_new)
  dist=dist2(as.matrix(SN),as.matrix(SN))
  W=affinityMatrix(dist, K, alpha)
  
  groupcy1F=spectralClustering(W, C_snf)
  
  cygroup1P1=numeric()
  cygroup2P1=numeric()
  cygroup3P1=numeric()
  
  group=groupcy1F
  for(g in 1:length(group)){
    if(group[g]==1){
      cygroup1P1=cbind(cygroup1P1,c(g,P1[g],P2[g]))}
    else if(group[g]==2){
      cygroup2P1=cbind(cygroup2P1,c(g,P1[g],P2[g]))}
    else if(group[g]==3){
      cygroup3P1=cbind(cygroup3P1,c(g,P1[g],P2[g]))}
  }
  
  t=t.test(cygroup1P1[2,], cygroup2P1[2,])
  print(t[c(3,5)])
  # print(ncol(cygroup1P1))
  # print(ncol(cygroup2P1))
  # print(cygroup1P1[1,])
  # print(cygroup2P1[1,])
  # print(length(get(paste("more",j,sep=""))))
  if(t[3] < 0.05){
    more_PO=c(more_PO,get(paste("more",j,sep="")))
  }
}
#print(more)

Data_new=Data[,sort(c(more_PO))]
SN=standardNormalization(Data_new)
dist=dist2(as.matrix(SN),as.matrix(SN))
W=affinityMatrix(dist, K, alpha)

groupcy1F=spectralClustering(W, C_snf)

cygroup1P1=numeric()
cygroup2P1=numeric()
cygroup3P1=numeric()

group=groupcy1F
for(g in 1:length(group)){
  if(group[g]==1){
    cygroup1P1=cbind(cygroup1P1,c(g,P1[g],P2[g]))}
  else if(group[g]==2){
    cygroup2P1=cbind(cygroup2P1,c(g,P1[g],P2[g]))}
  else if(group[g]==3){
    cygroup3P1=cbind(cygroup3P1,c(g,P1[g],P2[g]))}
}

print(length(more_PO))
t=t.test(cygroup1P1[2,], cygroup2P1[2,])
print(t[c(3,5)])

CyP1_PO=colnames(Data)[more_PO]
print(CyP1_PO)