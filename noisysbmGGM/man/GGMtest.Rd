% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data.R
\docType{data}
\name{GGMtest}
\alias{GGMtest}
\title{GGM for test}
\format{
\describe{
\item{\code{dataMatrix}}{A n-sample of a p Gaussian Vector associated to a GGM G}
\item{\code{Z.true}}{True latent clustering}
\item{\code{A.true}}{True latent adjacency matrix of the graph G}
}
}
\usage{
GGMtest
}
\description{
Example of a GGM
}
\examples{
main_noisySBM_GGM(GGMtest$dataMatrix,Meth="Ren",NIG=TRUE,Qup=10,nbOfZ=1,nbCores=1)

#Note : These data were created using the following instructions
n=30
p=10
u=0.1
v=0.3
theta=list(pi=c(1/3,2/3),w=0.25*cbind(c(1/6,1/120),c(1/120,1/6)))
Q=2
Z <- sample(1:Q, p, replace=TRUE, prob=theta$pi)
A <- matrix(0, p, p)
for (i in 1:(p-1)){
    A[i,(i+1):p] <- stats::rbinom(p-i, 1, theta$w[Z[i],Z[(i+1):p]])}
 A.true <- A + t(A)
Omega <- A.true*v
diag(Omega) = abs(min(eigen(Omega)$values)) + 0.1 + u
Sigma = stats::cov2cor(solve(Omega))
X = MASS::mvrnorm(n, rep(0, p), Sigma)
GGMtest=list(dataMatrix=X,Z.true=Z,A.true=A.true)
}
\keyword{datasets}
