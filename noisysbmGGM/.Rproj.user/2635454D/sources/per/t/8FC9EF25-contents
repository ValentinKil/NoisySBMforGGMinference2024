#We analysis the dataset Sachs
rm(list=ls())
library(parallel)
library(huge)  #si huge
library(MASS)  #si mvtnorm
library(SILGGM)
library(igraph)
library(noisysbmGGM)


# Define a custom layout function with user-specified coordinates
custom_layout <- function(graph, coordinates) {
  layout <- matrix(0, nrow = vcount(graph), ncol = 2)
  
  # Update the layout matrix with user-specified coordinates
  layout[, 1] <- coordinates[, 1]
  layout[, 2] <- coordinates[, 2]
  
  return(layout)
}

# Specify the user-defined coordinates for the vertices
user_coordinates <- matrix(c(0, 1,
                             1, 1,
                             -1.5, 0,
                             -1, -0.5,
                             -1.5, -1,
                             1.5,0,
                             1.5,-1,
                             0.5,0,
                             -0.5,0,
                             0,0.5,
                             0, -0.5), ncol = 2, byrow = TRUE)


Names=c("Raf","Mek1/2","PLCg","PIP2","PIP3","Erk1/2","Akt","PKA","PKC","p38","JNK")

A.true=matrix(0,nrow=11,ncol=11)
A.true[1,]=c(0,1,0,0,0,0,0,1,1,0,0)
A.true[2,]=c(0,0,0,0,0,1,0,0,0,0,0)
A.true[3,]=c(0,0,0,1,1,0,0,0,1,0,0)
A.true[4,]=c(0,0,0,0,1,0,0,0,1,0,0)
A.true[5,]=c(0,0,0,0,0,0,1,0,0,0,0)
A.true[6,]=c(0,0,0,0,0,0,1,0,0,0,0)
A.true[7,]=c(0,0,0,0,0,0,0,1,0,0,0)
A.true[8,]=c(0,0,0,0,0,1,0,0,0,0,0)
A.true[9,]=c(0,0,0,0,0,0,0,1,0,0,0)
A.true[10,]=c(0,0,0,0,0,0,0,1,1,0,0)
A.true[11,]=c(0,0,0,0,0,0,0,1,1,0,0)


A.true=A.true+t(A.true)
G.true=graph_from_adjacency_matrix(A.true, mode="undirected")



#--------------full dataset---------------
load("Sachs/Sachs.Rdata")
p=11
data=data.matrix(Sachs)
X=data

alpha=0.05
nbCores=6
nbOfZ=50
Qup=11

set.seed(1)
outlist <- SILGGM(X, method = "GFC_L", global=TRUE, true_graph = A.true, alpha=alpha)   #GFC_L = names of LiuL in SILGGM package
A=outlist$global_decision[[1]]
G = graph_from_adjacency_matrix(A, mode="undirected")
plot.igraph(G,main="LiuL-classical ",layout = custom_layout(G, user_coordinates), vertex.label=Names,vertex.shape="none" ,vertex.color=NA, edge.arrow.mode=0, edge.width=2)

dataMatrix =outlist$T_stat
result=main_noisySBM(dataMatrix, Qup=Qup,NIG=TRUE,threshold=0.5, Nbrepet=2, nbCores=nbCores, nbOfZ=nbOfZ, sigma0=1, alpha=alpha)
Ahat= result$A
Zhat= result$Z
Ghat=graph_from_adjacency_matrix(Ahat, mode="undirected")
plot.igraph(Ghat,main='LiuL-NSBM',layout = custom_layout(Ghat, user_coordinates), vertex.label=Names,vertex.shape="none" ,vertex.color=NA, edge.arrow.mode=0)


# ------------ n=20 ---------------
n=20
nbOfZ=50
Nsim=200
Aliu.vec = matrix(NA, nrow=Nsim, ncol=p*(p-1)/2)
Ahat.NIG.vec = matrix(NA, nrow=Nsim, ncol=p*(p-1)/2)

for (nsim in 1:Nsim){
  X=data[sample(1:902,n,replace=FALSE),]
  outlist_Liu_L <- SILGGM(X, method = "GFC_L", alpha=alpha)
  Aliu=outlist_Liu_L$global_decision[[1]]
  Aliu.vec[nsim,] = Aliu[lower.tri(diag(p))]
  
  dataMatrix =outlist_Liu_L$T_stat
  resultNIG=main_noisySBM(dataMatrix, Qup=Qup,NIG=TRUE,threshold=0.5, Nbrepet=2, nbCores=nbCores, nbOfZ=nbOfZ, sigma0=1, alpha=alpha)
  Ahat.NIG= resultNIG$A
  Ahat.NIG.vec[nsim,] = Ahat.NIG[lower.tri(diag(p))]
}

ind.all=listNodePairs(p)
Names[ind.all[,1]]
Names=c("Raf","Mek1/2","PLCg","PIP2","PIP3","Erk1/2","Akt","PKA","PKC","p38","JNK")
detected.edges.n20 = cbind(Names[ind.all[,1]], Names[ind.all[,2]],apply(Aliu.vec, 2, sum),  apply(Ahat.NIG.vec, 2, sum))
save(list=ls(), file="studySachsn20paper.Rdata")


load("studySachsn20paper.Rdata")
selectedEdges = c(1,20,21,28,41,42,46,53,54,55)   # edges of the benchmarck graph
detected.edges.n20[selectedEdges,]
[,1]     [,2]     [,3]  [,4] 
[1,] "Raf"    "Mek1/2" "200" "200"
[2,] "PLCg"   "PIP2"   "28"  "44" 
[3,] "PLCg"   "PIP3"   "121" "143"
[4,] "PIP2"   "PIP3"   "167" "172"
[5,] "Erk1/2" "Akt"    "199" "199"
[6,] "Erk1/2" "PKA"    "42"  "54" 
[7,] "Akt"    "PKA"    "92"  "121"
[8,] "PKC"    "p38"    "146" "158"
[9,] "PKC"    "JNK"    "107" "121"
[10,] "p38"    "JNK"    "102" "115"

