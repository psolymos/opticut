library(opticut)
library(dclone)
library(rjags)
source("~/repos/opticut/extras/calibrate.R")
#source("~/repos/opticut/R/ipredict.R")
#source("~/repos/opticut/R/ipredict.default.R")
#source("~/repos/opticut/R/ipredict.opticut.R")
source("~/repos/opticut/extras/multiclass.R")

#Inverse prediction might be a good metric of overall accuracy and scaling species contribution (i.e. the real indicator power).
#Implement binary-multinomial-continuous comparison.
#See how non-linear (unimodal, multimodal) responses affect inverse prediction
#Somehow show how # species affects accuracy.

## LOO

## define data

if (FALSE) {
Proj <- "dolina"
Dist <- "binomial"
data(dolina)
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
Y <- dolina$xtab[dolina$samp$method=="Q",]
X <- dolina$samp[dolina$samp$method=="Q",]
Y <- Y[,colSums(Y > 0) >= 20]
if (Dist == "binomial")
    Y <- ifelse(Y > 0, 1, 0)
s_col <- "mhab"
}

if (FALSE) {
Proj <- "simul"
Dist <- "binomial"
set.seed(1)
K <- 5
m <- 20
g <- rep(LETTERS[1:K], each=20)
n <- length(g)
cf <- matrix(rnorm(2*m), 2, m)
bp <- array(NA, c(K, m), list(LETTERS[1:K], paste0("Sp", 1:m)))
bpcf <- bp
for (j in 1:m) {
    bp[,j] <- rbinom(K, 1, rbeta(1,10,15))
    bpcf[,j] <- cf[bp[,j]+1,j]
}
mu <- bpcf[match(g, rownames(bp)),]
Y <- rbinom(length(mu), 1, plogis(drop(mu)))
dim(Y) <- dim(mu)
dimnames(Y) <- dimnames(mu)
X <- data.frame(g=g)
s_col <- "g"
}

cl <- makeCluster(4)
nn <- nrow(Y)
mm <- ncol(Y)
## All species + LOO
gnew0 <- character(0)
pm0 <- matrix(NA, 0, nlevels(X[,s_col]))
for (i in 1:nn) {
    if (interactive()) {
        cat("<<< All species --- Run", i, "of", nn, ">>>\n")
    }
    ii <- seq_len(nrow(Y)) != i
    y_trn <- Y[ii,,drop=FALSE]
    y_new <- Y[!ii,,drop=FALSE]
    x_trn <- X[ii,,drop=FALSE]
    #x_new <- X[!ii,,drop=FALSE]
    o <- opticut(y_trn ~ 1, strata=x_trn[,s_col], dist=Dist, cl=cl)
    cal <- calibrate(o, y_new, n.chains=length(cl), cl=cl)
    gnew0 <- c(gnew0, cal$gnew)
    pm0 <- rbind(pm0, cal$pi)
}
## -1 species + LOO
gnew_list <- list()
pm_list <- list()
for (j in 1:mm) {
    gnew <- character(0)
    pm <- matrix(NA, 0, nlevels(X[,s_col]))
    for (i in 1:nn) {
        if (interactive()) {
            cat("<<< -Spp", j, "of", mm, "--- Run", i, "of", nn, ">>>\n")
        }
        ii <- seq_len(nrow(Y)) != i
        jj <- seq_len(ncol(Y)) != j
        y_trn <- Y[ii,jj,drop=FALSE]
        y_new <- Y[!ii,jj,drop=FALSE]
        x_trn <- X[ii,,drop=FALSE]
        #x_new <- X[!ii,,drop=FALSE]
        o <- opticut(y_trn ~ 1, strata=x_trn[,s_col], dist=Dist, cl=cl)
        cal <- calibrate(o, y_new, n.chains=length(cl), cl=cl)
        gnew <- c(gnew, cal$gnew)
        pm <- rbind(pm, cal$pi)
    }
    gnew_list[[j]] <- gnew
    pm_list[[j]] <- pm
}
stopCluster(cl)

save(X, Y, s_col, gnew0, pm0, gnew_list, pm_list,
    file=file.path("~/Dropbox/collaborations/opticut/R",
    paste0("calibr-", Proj, "-", Dist, ".Rdata")))

mcm(X[1:nn,s_col], factor(gnew0, levels(X[,s_col])))
(tt <- table(X[1:nn,s_col], factor(gnew0, levels(X[,s_col]))))
sum(diag(tt))/sum(tt)

MCM0 <- mcm(X[,s_col], factor(gnew0, levels(X[,s_col])))
MCM <- lapply(gnew_list, function(z) mcm(X[,s_col], factor(z, levels(X[,s_col]))))

stat_fun <- function(x) x$accuracy

Stat0 <- stat_fun(MCM0)
Stat <- t(sapply(MCM, stat_fun))

Acc0 <- mean(Stat0)
Acc <- rowMeans(Stat)
names(Acc) <- colnames(Y)

plot(sort(Acc), type="b", pch=19)
abline(h=Acc0, col=2, lty=2)
dAcc <- Acc - Acc0

dAccH <- t(t(Stat) - Stat0)
rownames(dAccH) <- colnames(Y)

matplot(apply(Stat, 2, sort), type="l", lty=1, lwd=2)
for (i in 1:length(Stat0))
    abline(h=Stat0[i], lty=2, col=i)

## weighted aweraging

Proj <- "dolina"
Dist <- "binomial"
data(dolina)
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
Y <- dolina$xtab[dolina$samp$method=="Q",]
X <- dolina$samp[dolina$samp$method=="Q",]
Y <- Y[,colSums(Y > 0) >= 20]
if (Dist == "binomial")
    Y <- ifelse(Y > 0, 1, 0)
s_col <- "mhab"

i <- 111
cl <- NULL
ii <- seq_len(nrow(Y)) != i
y_trn <- Y[ii,,drop=FALSE]
y_new <- Y[!ii,,drop=FALSE]
x_trn <- X[ii,,drop=FALSE]
o <- opticut(y_trn ~ 1, strata=x_trn[,s_col], dist=Dist, cl=cl)
u <- uncertainty(o, type="multi", B=100)

ynew <- y_new
object <- o
wa <- function(object, ...)
    UseMethod("wa")
wa.opticut <- function(object, ynew) {
    s <- summary(object)
    bp <- s$bestpart
    I <- s$summary$I
#    p <- t(apply(ynew, 1, function(z) z / sum(z)))
#    t(apply(p, 1, function(y)
#        colMeans(apply(bp, 2, function(z) y * z * I))))
    t(apply(ynew, 1, function(y)
        colSums(apply(bp, 2, function(z)
            abs((I*z + (1-z)*(1-I)) - y)))))
#    t(apply(ynew, 1, function(y)
#        colSums(apply(bp, 2, function(z)
#            abs((I*z + (1-z)*(-I)) - y)))))
}
object <- u
wa.uncertainty <- function(object, ynew) {
    bbp <- lapply(object$uncertainty, bestpart)
    K <- nrow(bbp[[1L]])
    bp <- array(NA, c(length(bbp), K, object$B+1L))
    for (i in seq_len(length(bbp)))
        bp[i,,] <- bbp[[i]]
    I <- sapply(object$uncertainty, function(z) z$I)
    p <- t(apply(ynew, 1, function(z) z / sum(z)))
    fun <- function(i,j) {
        colMeans(apply(bp[,,j], 2, function(z) p[i,] * z * I[j,]))
    }
    xx <- lapply(seq_len(nrow(ynew)), function(i)
        sapply(seq_len(object$B+1L), function(j) fun(i,j)))
    mm <- lapply(xx, function(z) apply(z, 2, which.max))
    f <- function(x, K) {
        out <- numeric(K)
        for (k in 1:K)
            out[k] <- sum(x==k)
        out
    }
    out <- t(sapply(mm, f, K=K) / (object$B+1L))
    rownames(out) <- rownames(ynew)
    colnames(out) <- rownames(bbp[[1]])
    out
}
wa(o, y_new)
wa(u, y_new)

#cl <- makeCluster(4)
cl <- NULL
nn <- nrow(Y)
mm <- ncol(Y)
## All species + LOO
gnew0 <- character(nn)
pm0 <- matrix(NA, nn, nlevels(X[,s_col]))
gnew0u <- character(nn)
pm0u <- matrix(NA, nn, nlevels(X[,s_col]))
for (i in 1:nn) {
    if (interactive()) {
        cat("<<< All species --- Run", i, "of", nn, ">>>\n")
    }
    ii <- seq_len(nrow(Y)) != i
    y_trn <- Y[ii,,drop=FALSE]
    y_new <- Y[!ii,,drop=FALSE]
    x_trn <- X[ii,,drop=FALSE]
    o <- opticut(y_trn ~ 1, strata=x_trn[,s_col], dist=Dist, cl=cl)
#    u <- uncertainty(o, type="multi", B=100, cl=cl)
    cal <- wa(o, y_new)
#    calu <- wa(u, y_new)
    if (!any(is.na(cal)))
        gnew0[i] <- colnames(cal)[which.max(cal)]
    pm0[i,] <- cal[1,]
#    if (!any(is.na(cal)))
#        gnew0u[i] <- colnames(calu)[which.max(calu)]
#    pm0u[i,] <- calu[1,]
}
#stopCluster(cl)
iii <- gnew0 != ""
mcm(X[iii,s_col], factor(gnew0[iii], levels(X[,s_col])))
mcm(X[iii,s_col], factor(gnew0u[iii], levels(X[,s_col])))

<<<<<<< HEAD

## example
if (FALSE) {

library(opticut)
library(dclone)
library(rjags)
#source("~/repos/opticut/extras/calibrate.R")
source("~/repos/opticut/R/ipredict.R")
source("~/repos/opticut/R/ipredict.default.R")
source("~/repos/opticut/R/ipredict.opticut.R")
source("~/repos/opticut/extras/multiclass.R")

data(dolina)
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
Y <- dolina$xtab[dolina$samp$method=="Q",]
X <- dolina$samp[dolina$samp$method=="Q",]
Y <- ifelse(Y > 0, 1, 0)

set.seed(i)
i <- sample.int(nrow(Y), 100)
ii <- seq_len(nrow(Y)) %in% i
y_trn <- Y[ii,,drop=FALSE]
y_trn <- y_trn[,colSums(y_trn > 0) >= 20,drop=FALSE]
y_new <- Y[!ii,colnames(y_trn),drop=FALSE]
x_trn <- X[ii,,drop=FALSE]
x_new <- X[!ii,,drop=FALSE]

oc1 <- opticut(y_trn ~ 1, strata=x_trn[,"mhab"], dist=Dist)
ip1 <- ipredict(oc1, y_new, n.iter=1000)

oc2 <- opticut(y_trn ~ lmoist, data=x_trn, strata=x_trn[,"mhab"], dist=Dist)
ip2 <- ipredict(oc2, y_new, x_new, n.iter=1000)

mod1 <- lapply(1:ncol(y_trn), function(i) {
    glm(y_trn[,i] ~ mhab, data=x_trn, family=binomial)
})
ip3 <- ipredict(mod1, y_new, K=4, n.iter=1000)

mod2 <- lapply(1:ncol(y_trn), function(i) {
    glm(y_trn[,i] ~ mhab + lmoist, data=x_trn, family=binomial)
})
ip4 <- ipredict(mod2, y_new, xnew=model.matrix(~ lmoist, x_new), K=4, n.iter=1000)

}
=======
colnames(pm0) <- colnames(cal)
gnew0 <- mefa4::find_min(pm0)
mcm(X[iii,s_col], factor(gnew0[,1], levels(X[,s_col])))

(tt <- table(data=X[iii,s_col], est=factor(gnew0[,1], levels(X[,s_col]))))

load("~/Dropbox/collaborations/opticut/R/calibr-dolina-binomial.Rdata")
(tt <- table(data=X[,s_col], class=factor(gnew0, levels(X[,s_col]))))

## InvPred
#    class
#data LI DW TL RO
#  LI 71 11 22  8
#  DW 10 14 14 10
#  TL 30  7  3  8
#  RO  4  6  3 35

## -I
#    est
#data LI DW TL RO
#  LI 26  6  0 80
#  DW  2  9  0 37
#  TL  6  9  0 33
#  RO  1  3  0 44

## 1-I
#    est
#data  LI  DW  TL  RO
#  LI   3   7 102   0
#  DW   5  15  28   0
#  TL   8   5  34   1
#  RO   4  17  16  11
>>>>>>> master
