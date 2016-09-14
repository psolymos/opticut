\name{opticut-package}
\alias{opticut-package}
\docType{package}
\title{
\packageTitle{opticut}
}
\description{
\packageDescription{opticut}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{opticut}
\packageIndices{opticut}

The main user interface is the \code{\link{opticut}} function
to find the optimal binary partitioning.

\code{\link{optilevel}} finds the optimal number of factor levels.
}
\author{
\packageAuthor{opticut}

Maintainer: \packageMaintainer{opticut}
}
\references{
Kemencei, Z., Farkas, R., Pall-Gergely, B., Vilisics, F., Nagy, A., Hornung,
E. & Solymos, P. (2014): Microhabitat associations of land snails in
forested dolinas: implications for coarse filter conservation.
Community Ecology 15:180--186.
}
\examples{
## simple example from Legendre 2013
## Indicator Species: Computation, in
## Encyclopedia of Biodiversity, Volume 4
## http://dx.doi.org/10.1016/B978-0-12-384719-5.00430-5
gr <- as.factor(paste0("X", rep(1:5, each=5)))
spp <- cbind(Species1=rep(c(4,6,5,3,2), each=5),
    Species2=c(rep(c(8,4,6), each=5), 4,4,2, rep(0,7)),
    Species3=rep(c(18,2,0,0,0), each=5))
rownames(spp) <- gr
spp
summary(oc <- opticut(spp ~ 1, strata=gr, comb="all"))
plot(oc)

ol <- optilevels(spp[,"Species3"], as.factor(gr))
ol[c("delta", "coef", "rank", "levels")]
}
\keyword{ package }
