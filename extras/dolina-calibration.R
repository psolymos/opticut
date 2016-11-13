data(dolina)
library(parallel)
library(mefa4)


Y <- ifelse(dolina$xtab[,colSums(dolina$xtab > 0) >= 20] > 0, 1, 0)
x <- dolina$samp
i <- 1:nrow(Y) %in% sample.int(nrow(Y), round(0.9*nrow(Y)))
Yxv <- Y[!i,]
xxv <- x[!i,]
Y <- Y[i,]
x <- x[i,]

dol <- opticut(Y ~ 1, data=x, strata=x$mhab, dist="binomial")

cl <- makeCluster(2)
ucdol <- uncertainty(dol, type="multi", B=25)
stopCluster(cl)

bp <- bestpart(ucdol)

m <- sapply(rownames(bp), function(i) colSums(t(Yxv) * bp[i,]))
df <- data.frame(find_max(m), mhab=xxv$mhab)
tt <- table(pred=df$index, true=df$mhab)
tt <- tt[rownames(bp),rownames(bp)]
tt

## need to apply I and selection values for each iteration combined with detections
## then see how many times one was better supported

