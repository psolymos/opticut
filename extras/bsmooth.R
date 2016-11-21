data(dolina)
## stratum as ordinal
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
## filter species to speed up things a bit
Y <- ifelse(dolina$xtab[,colSums(dolina$xtab > 0) >= 20] > 0, 1, 0)
## opticut results, note the cloglog link function
#dol <- opticut(Y ~ stratum + lmoist + method, data=dolina$samp,
#    strata=dolina$samp$mhab, dist="binomial:cloglog")
dol <- opticut(Y, strata=dolina$samp$mhab, dist="binomial:cloglog")

## parallel computing for uncertainty
library(parallel)
cl <- makeCluster(2)
object <- uncertainty(dol, type="multi", B=25, cl=cl)
stopCluster(cl)


strata.uncertainty <-
function (object, ...) {
    strata.opticut(object, ...)
}

bsmooth <- function (object, ...)
    UseMethod("bsmooth")

bsmooth.uncertainty <- function (object, ...) {
    st <- as.character(strata(object))
    out <- object$uncertainty
    for (i in seq_len(length(out))) {
        uc1 <- object$uncertainty[[i]]
        bp1 <- bestpart(object$uncertainty[[i]])
        m <- matrix(NA, nobs(object), ncol(bp1))
        for (j in seq_len(nobs(object))) {
            m[j,] <- bp1[st[j],] * uc1$mu1 + (1 - bp1[st[j],]) * uc1$mu0
        }
        out[[i]] <- m
    }
    out
}

bs <- bsmooth(object)
boxplot(rowMeans(bs[[1]]) ~ strata(object))

xx <- sapply(bs, rowMeans)
xx <- mefa4::groupMeans(xx, 1, strata(object))
heatmap(t(xx), scale="none", col=occolors()(25),
    distfun=function(x) dist(x, "manhattan"))

yy <- mefa4::groupMeans(Y, 1, strata(object))
plot(yy, xx);abline(0,1)

