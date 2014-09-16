## Program to calculate the type of copula for each pair of variables and then estimate
## parameters, given data
library(MASS)
library(VineCopula)
library(copula)
library(network)
library(qgraph)

setwd('C:/study/my projects/depth')
source('functions.R')

## Simulated data
# n = 1000; p = 2
# X = BiCopSim(n, family=23, par=-3)

# real data: Breast Tissue data, UCI repository: https://archive.ics.uci.edu/ml/datasets/Breast+Tissue
BreastTissue = read.csv("C:/Study/My projects/depth/BreastTissue.csv", header=T)
X = BreastTissue[,3:11]
p=9; npairs=p*(p-1)/2
U = pobs(X); pairs(U)

graph.ml = CopNetwork(X, is.U=F, method='mle', quiet=F)
plot.CopNetwork(graph.ml, varnames=colnames(X),"Method- maximum likelihood", leg=F)

graph.md = CopNetwork(X, is.U=F, method='mde', quiet=F)
plot.CopNetwork(graph.md, varnames=colnames(X),"Method- maximum likelihood depth", leg=T)

cbind(graph.ml,graph.md[,3:5])


# Cardiotocography data: https://archive.ics.uci.edu/ml/datasets/Cardiotocography#
CTG = read.csv("C:/Study/My projects/depth/CTG.csv", header=T)
X = CTG[,-c(6,22)]
U = pobs(X)

#graph.ml.ct = CopNetwork(X, is.U=F, method='mle', quiet=F)
graph.ml.ct = read.table("ctgraph.txt", header=T, quote="\"")
plot.CopNetwork(graph.ml.ct, varnames=colnames(X), "Method- maximum likelihood")

## summarize full graph
summary.CopNetwork(graph.ml.ct)

## class-wise analysis
g1 = X[which(CTG$NSP==1),]
g1a.mlct = CopNetwork(g1[,1:6], is.U=F, method='mle', quiet=F)
plot.CopNetwork(g1a.mlct, varnames=colnames(X[1:6]), "Fetal class - normal, group A variables",leg=F)
g1b.mlct = CopNetwork(g1[,7:10], is.U=F, method='mle', quiet=F)
plot.CopNetwork(g1b.mlct, varnames=colnames(X[7:10]), "Fetal class - normal, group B variables",leg=F)
g1c.mlct = CopNetwork(g1[,11:20], is.U=F, method='mle', quiet=F)
g1c.mlct[,1:2] = g1c.mlct[,1:2]+10
summary.CopNetwork(g1c.mlct)

g2 = X[which(CTG$NSP==2),]
g2a.mlct = CopNetwork(g2[,1:6], is.U=F, method='mle', quiet=F)
plot.CopNetwork(g2a.mlct, varnames=colnames(X[1:6]), "Fetal class - suspect, group A variables",leg=F)
g2b.mlct = CopNetwork(g2[,7:10], is.U=F, method='mle', quiet=F)
plot.CopNetwork(g2b.mlct, varnames=colnames(X[7:10]), "Fetal class - suspect, group B variables",leg=F)
g2c.mlct = CopNetwork(g2[,11:20], is.U=F, method='mle', quiet=F)
g2c.mlct[,1:2] = g2c.mlct[,1:2]+10
summary.CopNetwork(g2c.mlct)

g3 = X[which(CTG$NSP==3),]
g3a.mlct = CopNetwork(g3[,1:6], is.U=F, method='mle', quiet=F)
plot.CopNetwork(g3a.mlct, varnames=colnames(X[1:6]), "Fetal class - pathologic, group A variables",leg=T)
g3b.mlct = CopNetwork(g3[,7:10], is.U=F, method='mle', quiet=F)
plot.CopNetwork(g3b.mlct, varnames=colnames(X[7:10]), "Fetal class - pathologic, group B variables",leg=T)
g3c.mlct = CopNetwork(g3[,11:20], is.U=F, method='mle', quiet=F)
g3c.mlct[,1:2] = g3c.mlct[,1:2]+10
summary.CopNetwork(g3c.mlct)
