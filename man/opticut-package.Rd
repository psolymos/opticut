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
Make sure to evaluate \code{\link{uncertainty}}.

Go beyond binary partitions:
\code{\link{multicut}} fits multinomial-response models;
\code{\link{optilevels}} finds the optimal number of factor levels.
}
\author{
\packageAuthor{opticut}

Maintainer: \packageMaintainer{opticut}
}
\references{
Kemencei, Z., Farkas, R., Pall-Gergely, B., Vilisics, F., Nagy, A., Hornung,
E. & Solymos, P., 2014. Microhabitat associations of land snails in
forested dolinas: implications for coarse filter conservation.
Community Ecology 15:180--186.
<doi:10.1556/ComEc.15.2014.2.6>
}
\examples{
## community data
y <- cbind(
    Sp1=c(4,6,3,5, 5,6,3,4, 4,1,3,2),
    Sp2=c(0,0,0,0, 1,0,0,1, 4,2,3,4),
    Sp3=c(0,0,3,0, 2,3,0,5, 5,6,3,4))

## stratification
g <-    c(1,1,1,1, 2,2,2,2, 3,3,3,3)

## find optimal partitions for each species
oc <- opticut(formula = y ~ 1, strata = g, dist = "poisson")
summary(oc)

## visualize the results
plot(oc, cut = -Inf)

## quantify uncertainty
uc <- uncertainty(oc, type = "asymp", B = 999)
summary(uc)

## go beyond binary partitions

mc <- multicut(formula = y ~ 1, strata = g, dist = "poisson")
summary(mc)

ol <- optilevels(y[,"Sp2"], as.factor(g))
ol[c("delta", "coef", "rank", "levels")]
}
\keyword{ package }
