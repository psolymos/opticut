## This is opticut based calibration
##
## 1. simulate some data
if (!require(opticut)) {
    if (!require(devtools))
        install.packages("devtools")
    devtools::install_github("psolymos/opticut")
}
library(opticut)
library(mefa4)

S <- 20 # number of species
N <- 200 # number of possible sampling loactions
n <- 100 # number of sampling locations out of N

set.seed(1234321)
## random coefs for species responses
cfmat <- matrix(rnorm(S*3), S, 3)
## predictor variable
x <- sort(runif(N))
X <- model.matrix(~x + I(x^2))

## measured x
xm <- x + 0 # rnorm(N, 0, 0.1)
br <- c(-Inf, 0.25, 0.5, 0.75, Inf)
#br <- c(-Inf, 0.33, 0.66, Inf)
#br <- c(-Inf, 0.5, Inf)
table(cut(x, br), cut(xm, br))

mu <- apply(cfmat, 1, function(z) X %*% z)

matplot(x, exp(mu), type="l", lty=1)

Y <- rpois(N*S, exp(mu))
dim(Y) <- dim(mu)
colnames(Y) <- paste0("Sp", 1:S)

samp_i <- sample.int(N, n)
uns_i <- seq_len(N)[!(seq_len(N) %in% samp_i)]
z <- cut(xm, br)
levels(z) <- LETTERS[1:nlevels(z)]

Ysamp <- Y[samp_i,]
zsamp <- z[samp_i]
ocop <- ocoptions(collapse="+", sort=FALSE, cut=-Inf)
oc <- opticut(Ysamp ~ 1, strata=zsamp, dist="poisson", comb="rank")
#plot(oc,sort=TRUE)

ocb <- uncertainty(oc, type="multi")

getProps <- function(z) {
    LEV <- names(attr(z, "est"))
    LAB <- strsplit(levels(z$best), attr(z, "collapse"), fixed=TRUE)
    LAB <- LAB[as.integer(z$best)]
    TAB <- table(unlist(LAB))
    out <- numeric(length(LEV))
    names(out) <- LEV
    out[names(TAB)] <- TAB
    out / nrow(z)
}

get_tmat <- function(object, useI=FALSE, lo=0, hi=1) {
    LEV <- names(attr(object, "est"))
    LAB <- strsplit(levels(object$best), attr(object, "collapse"), fixed=TRUE)
    LAB <- LAB[as.integer(object$best)]
    B <- length(LAB)
    tmat <- matrix(lo, B, length(LEV))
    colnames(tmat) <- LEV
    for (i in seq_len(B))
        tmat[i,LAB[[i]]] <- hi
    if (useI)
        tmat <- tmat * object$I
    tmat
}

## new observations (z unknown)
Yuns <- Y[uns_i,]
zuns <- z[uns_i]

## weights
bp <- summary(oc)$bestpart
#bp1 <- (bp*2 -1) * summary(oc)$summary$I
#bp1 <- ifelse(bp > 0, 0, -1) * summary(oc)$summary$I
bp1 <- bp * summary(oc)$summary$I

w <- t(sapply(ocb$uncertainty, getProps))
#w2 <- t(sapply(ocb$uncertainty, function(z) colMeans(get_tmat(z, useI=TRUE, lo=-1, hi=1))))
#w2 <- t(sapply(ocb$uncertainty, function(z) colMeans(get_tmat(z, useI=TRUE, lo=-1, hi=0))))
w2 <- t(sapply(ocb$uncertainty, function(z) colMeans(get_tmat(z, useI=TRUE, lo=0, hi=1))))



#### Weighted averaging, correlation, or distance based ranking
## averaging and correlation works quite well
## distance sucks

## y is a site x species matrix (rows standardized to add up to 1)
## w is a species x partition matrix with weights (columns sum to 1)
wavg_fun <- function(y, w, type="avg", ...) {
    if (is.function(type)) {
        dist_fun <- type
        type <- "dis"
    } else {
        type <- match.arg(type, c("avg", "cor", "dis"))
        if (type == "dis")
            dist_fun <- stats::dist
    }
    p <- apply(y, 1, function(z) z / sum(z)) # transposed
    if (type == "avg")
        out <- sapply(seq_len(ncol(w)), function(k, ...) {
            colSums(p * w[,k], ...) / sum(w[,k], ...)
        })
    if (type == "cor")
        out <- sapply(seq_len(ncol(w)), function(k, ...) {
            sapply(seq_len(ncol(p)), function(i, ...) {
                cor(p[,i], w[,k], ...)
            })
        })
    if (type == "dis")
        out <- -sapply(seq_len(ncol(w)), function(k, ...) {
            sapply(seq_len(ncol(p)), function(i, ...) {
                as.numeric(dist_fun(rbind(p[,i], w[,k]), ...))
            })
        })
    dimnames(out) <- list(rownames(y), colnames(w))
    out
}

TYPE <- "cor"

## pure partitioning (0/1)
res <- wavg_fun(Yuns, bp, type=TYPE)
resv <- find_max(res)
table(pred=resv$index, zuns)
100*sum(diag(table(pred=resv$index, zuns))) / nrow(Yuns)

## I contrast
res1 <- wavg_fun(Yuns, bp1, type=TYPE)
resv1 <- find_max(res1)
table(pred=resv1$index, zuns)
100*sum(diag(table(pred=resv1$index, zuns))) / nrow(Yuns)

## bootstrap based support w
res2 <- wavg_fun(Yuns, w, type=TYPE)
resv2 <- find_max(res2)
table(pred=resv2$index, zuns)
100*sum(diag(table(pred=resv2$index, zuns))) / nrow(Yuns)

## w * I
res3 <- wavg_fun(Yuns, w2, type=TYPE)
resv3 <- find_max(res3)
table(pred=resv3$index, zuns)
100*sum(diag(table(pred=resv3$index, zuns))) / nrow(Yuns)



## cross validation

B <- 5
n <- nrow(Ysamp)
#oc <- opticut(Y ~ 1, strata=z, dist="poisson", comb="rank")
formula <- Y ~ 1

bv <- sample(seq_len(B), n, replace=TRUE)

Ilist <- list()
tlist <- list()
xvavg <- list()
xvcor <- list()

#bi <- 2
for (bi in seq_len(B)) {

Ytr <- Ysamp[bv != bi,]
Ztr <- zsamp[bv != bi]
Yva <- Ysamp[bv == bi,]
Zva <- zsamp[bv == bi]

cat("fold", bi, "of", B, "\n")
#oci <- opticut(Ytr ~ 1, strata=Ztr, dist="poisson", comb="rank")
oci <- opticut(Ysamp ~ 1, strata=zsamp, dist="poisson", comb="rank",
    sset=bv != bi)
Ilist[[bi]] <- summary(oci)$bestpart
tlist[[bi]] <- summary(oci)$summary$I
names(tlist[[bi]]) <- rownames(summary(oci)$summary)

ri_avg <- wavg_fun(Yva, Ilist[[bi]], type="avg")
ri_cor <- wavg_fun(Yva, Ilist[[bi]], type="cor")
pi_avg <- factor(colnames(ri_avg)[apply(ri_avg, 1, which.max)], colnames(ri_avg))
pi_cor <- factor(colnames(ri_cor)[apply(ri_cor, 1, which.max)], colnames(ri_cor))

mc_avg <- table(True=Zva, Pred=pi_avg)
attr(mc_avg, "UA") <- diag(mc_avg) / colSums(mc_avg)
attr(mc_avg, "PA") <- diag(mc_avg) / rowSums(mc_avg)
attr(mc_avg, "OA") <- sum(diag(mc_avg)) / sum(mc_avg)
mc_cor <- table(True=Zva, Pred=pi_cor)
attr(mc_cor, "UA") <- diag(mc_cor) / colSums(mc_cor)
attr(mc_cor, "PA") <- diag(mc_cor) / rowSums(mc_cor)
attr(mc_cor, "OA") <- sum(diag(mc_cor)) / sum(mc_cor)

xvavg[[bi]] <- mc_avg
xvcor[[bi]] <- mc_cor

#unclass(mc_avg)
}

for (bi in 1:B) {
    if (bi == 1) {
        xx <- xvavg[[bi]]
    } else {
        xx <- xx + xvavg[[bi]]
    }
}
attr(xx, "UA") <- diag(xx) / colSums(xx)
attr(xx, "PA") <- diag(xx) / rowSums(xx)
attr(xx, "OA") <- sum(diag(xx)) / sum(xx)
unclass(xx)

for (bi in 1:B) {
    if (bi == 1) {
        xx <- xvcor[[bi]]
    } else {
        xx <- xx + xvcor[[bi]]
    }
}
attr(xx, "UA") <- diag(xx) / colSums(xx)
attr(xx, "PA") <- diag(xx) / rowSums(xx)
attr(xx, "OA") <- sum(diag(xx)) / sum(xx)
unclass(xx)

## OA = overall accuracy
## UA = user's accuracy (reliability, error of commission)
##    prop of cases predicted represents the true class
## PA = producer's accuracy (error of omission)
##    prop of true class being classified the same under predicted

## we can use bootstrap results to make B prediction
## this could give the uncertainty in prediction

## cross validation: dolina data

data(dolina)
## stratum as ordinal
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
## filter species to speed up things a bit
Y <- dolina$xtab[,colSums(dolina$xtab > 0) >= 20]
Y01 <- ifelse(Y>0, 1, 0)
## opticut results, note the cloglog link function
dol <- opticut(Y01 ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="poisson")
summary(dol)
plot(dol)

B <- 5
n <- nrow(Y01)
bv <- sample(seq_len(B), n, replace=TRUE)

Ilist <- list()
tlist <- list()
xvavg <- list()
xvcor <- list()

#bi <- 2
for (bi in seq_len(B)) {

cat("fold", bi, "of", B, "\n")
oci <- opticut(Y01 ~ stratum + lmoist + method, data=dolina$samp,
    strata=dolina$samp$mhab, dist="poisson", sset=bv != bi)
Ilist[[bi]] <- summary(oci)$bestpart
tlist[[bi]] <- summary(oci)$summary$I
names(tlist[[bi]]) <- rownames(summary(oci)$summary)

ss <- bv == bi & rowSums(Y01) > 0

ri_avg <- wavg_fun(Y01[ss,,drop=FALSE], Ilist[[bi]], type="avg")
ri_cor <- wavg_fun(Y01[ss,,drop=FALSE], Ilist[[bi]], type="cor")
pi_avg <- factor(colnames(ri_avg)[apply(ri_avg, 1, which.max)], colnames(ri_avg))
pi_cor <- factor(colnames(ri_cor)[apply(ri_cor, 1, which.max)], colnames(ri_cor))

mc_avg <- table(True=as.character(dolina$samp$mhab[ss]), Pred=pi_avg)
attr(mc_avg, "UA") <- diag(mc_avg) / colSums(mc_avg)
attr(mc_avg, "PA") <- diag(mc_avg) / rowSums(mc_avg)
attr(mc_avg, "OA") <- sum(diag(mc_avg)) / sum(mc_avg)
mc_cor <- table(True=as.character(dolina$samp$mhab[ss]), Pred=pi_cor)
attr(mc_cor, "UA") <- diag(mc_cor) / colSums(mc_cor)
attr(mc_cor, "PA") <- diag(mc_cor) / rowSums(mc_cor)
attr(mc_cor, "OA") <- sum(diag(mc_cor)) / sum(mc_cor)

#unclass(mc_avg)
#unclass(mc_cor)

xvavg[[bi]] <- mc_avg
xvcor[[bi]] <- mc_cor

#unclass(mc_avg)
}

for (bi in 1:B) {
    if (bi == 1) {
        xx <- xvavg[[bi]]
    } else {
        xx <- xx + xvavg[[bi]]
    }
}
attr(xx, "UA") <- diag(xx) / colSums(xx)
attr(xx, "PA") <- diag(xx) / rowSums(xx)
attr(xx, "OA") <- sum(diag(xx)) / sum(xx)
unclass(xx)

for (bi in 1:B) {
    if (bi == 1) {
        xx <- xvcor[[bi]]
    } else {
        xx <- xx + xvcor[[bi]]
    }
}
attr(xx, "UA") <- diag(xx) / colSums(xx)
attr(xx, "PA") <- diag(xx) / rowSums(xx)
attr(xx, "OA") <- sum(diag(xx)) / sum(xx)
unclass(xx)

## zero samples throw things off
## soring the confusion matrix and making sure that missing classes are represented

