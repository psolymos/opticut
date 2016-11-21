library(opticut)

data(dolina)
ddata <- dolina$samp[dolina$samp$method == "T",]
Y <- dolina$xtab[dolina$samp$method == "T",]
Y <- ifelse(Y[,colSums(Y > 0) >= 20] > 0, 1, 0)
dol <- opticut(Y, strata=ddata$mhab, dist="binomial")

## parallel computing for uncertainty
library(parallel)
cl <- makeCluster(2)
object <- uncertainty(dol, type="multi", B=25, cl=cl)
stopCluster(cl)


strata.opticut <- opticut:::strata.opticut
strata.uncertainty <-
function (object, ...) {
    strata.opticut(object, ...)
}

bsmooth <- function (object, ...)
    UseMethod("bsmooth")

bsmooth.uncertainty1 <- function (object, ...) {
    bp <- bestpart(object)
    t(t(bp) * object$mu1 + (1 - t(bp)) * object$mu0)
}
bsmooth.uncertainty <- function (object, ...) {
    lapply(object$uncertainty, bsmooth)
}

yy <- mefa4::groupMeans(Y, 1, strata(object))
boxplot(t(bsmooth(object$uncertainty[["pvic"]])), range=0)
points(yy[,"pvic"], col=2, pch=4, cex=2)

xx <- sapply(bs, rowMeans)
xx <- mefa4::groupMeans(xx, 1, strata(object))
heatmap(t(xx), scale="none", col=occolors()(25),
    distfun=function(x) dist(x, "manhattan"))

plot(yy, xx, xlim=c(0,1), ylim=c(0,1));abline(0,1)

