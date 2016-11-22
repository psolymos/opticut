#devtools::install_github("psolymos/opticut")
library(opticut)

data(dolina)
ddata <- dolina$samp[dolina$samp$method == "T",]
Y <- dolina$xtab[dolina$samp$method == "T",]
Y <- ifelse(Y[,colSums(Y > 0) >= 20] > 0, 1, 0)
dol <- opticut(Y, strata=ddata$mhab, dist="binomial")
library(parallel)
cl <- makeCluster(2)
object <- uncertainty(dol, type="multi", B=25, cl=cl)
stopCluster(cl)




yy <- mefa4::groupMeans(Y, 1, strata(object))
boxplot(t(bsmooth(object$uncertainty[["pvic"]])), range=0)
points(yy[,"pvic"], col=2, pch=4, cex=2)

bsmooth(object)

xx <- sapply(bs, rowMeans)
xx <- mefa4::groupMeans(xx, 1, strata(object))
heatmap(t(xx), scale="none", col=occolors()(25),
    distfun=function(x) dist(x, "manhattan"))

plot(yy, xx, xlim=c(0,1), ylim=c(0,1));abline(0,1)

