library(opticut)
library(dclone)
library(rjags)
source("~/repos/opticut/extras/calibrate.R")

#Inverse prediction might be a good metric of overall accuracy and scaling species contribution (i.e. the real indicator power).
#Implement binary-multinomial-continuous comparison.
#See how non-linear (unimodal, multimodal) responses affect inverse prediction
#Somehow show how # species affects accuracy.

## Multiclass classification measures
## x is data class, y is classification

## simple binary confusion matrix
## Might have to check valid levels (x, y, union)
#cmat(c(1,0,0,1,1,1,0,0,0,0), c(1,1,1,0,0,0,0,0,0,0))
#cmat(c(1,0,0,1,1,1,0,0,0,0), c(1,1,1,0,0,0,0,0,0,0), drop=FALSE)

ctable <- function(x, y) {
    if (length(y) != length(x))
        stop("Hey! x and y lengths must match")
    x <- as.factor(x)
    x <- droplevels(x)
    levs <- sort(levels(x))
    x <- factor(x, levels=levs)
    y <- factor(y, levels=levs)
    if (any(is.na(y)))
        stop("NAs, and levels in y that are not in x are not allowed, sorry")
    as.matrix(table(Reference=x, Prediction=y))
}
btable <- function(table) {
    if (ncol(table) != nrow(table))
        stop("dimension mismatch")
    if (!all(colnames(table) == rownames(table)))
        stop("dimnames must match")
    f <- function(i, x) {
        c(tp=sum(x[i,i]), fp=sum(x[-i,i]),
        fn=sum(x[i,-i]), tn=sum(x[-i,-i]))
    }
    out <- sapply(seq_len(ncol(table)), f, x=table)
    colnames(out) <- colnames(table)
    out
}
untable <- function(table) {
    x <- y <- integer(sum(table))
    S <- c(0, cumsum(table))
    i <- row(table)
    j <- col(table)
    for (k in 2:length(S)) {
        if (S[k]-S[k-1] > 0) {
            x[(S[k-1]+1):S[k]] <- i[k-1]
            y[(S[k-1]+1):S[k]] <- j[k-1]
        }
    }
    list(x = factor(rownames(table)[x], levels=rownames(table)),
        y = factor(colnames(table)[y], levels=colnames(table)))
}

multiclass <- function(x, y=NULL, beta=1) {
    if (!is.null(y))
        x <- ctable(x, y)
    cm <- btable(x)
    Stat <- c(
        Acc = mean((cm["tp",] + cm["tn",]) / N),
        Err = mean((cm["fp",] + cm["fn",]) / N),
        Prec_m = sum(cm["tp",]) / sum(cm["tp",] + cm["fp",]),
        Prec_M = mean(cm["tp",] / (cm["tp",] + cm["fp",])),
        Spec_m = sum(cm["tn",]) / sum(cm["fp",] + cm["tn",]),
        Spec_M = mean(cm["tn",] / (cm["fp",] + cm["tn",])),
        Rec_m = sum(cm["tp",]) / sum(cm["tp",] + cm["fn",]),
        Rec_M = mean(cm["tp",] / (cm["tp",] + cm["fn",])))
    Stat <- c(Stat,
        ## F-score is harmonic mean
        ## G-score is geometric mean sqrt(prec * recall)
        F_m = unname((beta^2 + 1) * Stat["Prec_m"] * Stat["Rec_m"] /
            (beta^2 * Stat["Prec_m"] + Stat["Rec_m"])),
        F_M = unname((beta^2 + 1) * Stat["Prec_M"] * Stat["Rec_M"] /
            (beta^2 * Stat["Prec_M"] + Stat["Rec_M"])))
    Mat <- rbind(Accuracy=(cm["tp",] + cm["tn",]) / N,
        Error=(cm["fp",] + cm["fn",]) / N,
        Precision=cm["tp",] / (cm["tp",] + cm["fp",]),
        Specificity=cm["tn",] / (cm["fp",] + cm["tn",]),
        Recall=cm["tp",] / (cm["tp",] + cm["fn",])) # recall = sensitivity
    #Mat <- rbind(Mat, jouden=Mat["Specificity",]+Mat["Recall",]-1)
    out <- list(ctable=x, btable=cm, average=Stat, beta=beta,
        single=Mat)
    class(out) <- "multiclass"
    out
}


N <- 100
K <- 5
x <- sample(LETTERS[1:K], N, replace=TRUE, prob=sqrt(2^(0:(K-1))))
y <- x
for (i in 1:length(x))
    if (runif(1) > 0.1)
        y[i] <- sample(LETTERS[1:K], 1)
(ct <- ctable(x, y))
(cm <- btable(ct))
untable(ct)
stopifnot(all(ctable(untable(ct)$x, untable(ct)$y) == ctable(x, y)))
(x <- multiclass(x, y))

#https://en.wikipedia.org/wiki/Accuracy_paradox

etable <- function(table, type="cohen", w=NULL) {
    type <- match.arg(type,
        c("majority", "random", "weighted", "cohen"))
    N <- sum(table)
    K <- ncol(table)
    p <- rowSums(table) / N
    q <- colSums(table) / N
    if (type == "majority") {
        ## mcc: majority class classifier, No Information Rate (NIR)
        ecm <- table
        ecm[] <- 0
        ecm[,which.max(p)] <- rowSums(table)
    }
    if (type == "random") {
        ## rgc: random-guess classifier
        ecm <- table
        ecm[] <- N/K * p
    }
    if (type == "weighted") {
        ## rwgc: random-weighted-guess classifier
        ecm <- table
        if (is.null(w)) {
            w <- p
        } else {
            w <- w / sum(w)
        }
        ecm[] <- N * p %*% t(p)
    }
    if (type == "cohen") {
        ## Cohen
        ecm <- table
        ecm[] <- N * p %*% t(q)
    }
    ecm
}
rtable <- function(n, table, type=c("r", "rc")) {
    type <- match.arg(type)
    if (type == "rc") {
        out <- array(unlist(r2dtable(n, rowSums(table), colSums(table))),
            c(dim(table), n))
    }
    if (type == "r") {
        r <- rowSums(table)
        K <- length(r)
        f <- function() {
            t(sapply(seq_len(K), function(i)
                rmultinom(1, r[i], table[i,])))
        }
        out <- replicate(n, f())
    }
    #out <- rmultinom(n, sum(table), table) # this does not keep margins
    dimnames(out) <- list(rownames(table), colnames(table), NULL)
    out
}

kappa <- function(predicted, reference) {
    a <- sum(diag(predicted)) / sum(predicted)
    a0 <- sum(diag(reference)) / sum(reference)
    k <- (a - a0) / (1 - a0)
    c(a=a, a0=a0, k=k)
}

eval_fun <- function(table, n=0, ptype="cohen", rtype="r", w=NULL) {
    ref <- etable(table, type=ptype, w=w)
    D <- c(dim(table), n+1)
    rnd <- if (n > 0)
        c(table, rtable(n, ref, type=rtype)) else c(table)
    dim(rnd) <- D
    k <- sapply(seq_len(n+1), function(i) kappa(table, rnd[,,i]))
    k
}

etable(x$ctable)

sum(diag(ecm))/N
sum(p*q)
sum(diag(table))/N
kappa <- (sum(diag(table))/N - sum(diag(ecm))/N) / (1 - sum(diag(ecm))/N)

rmn <- rmultinom(10^4, sum(aa), etable(aa, "c"))
v <- apply(rmn, 2, function(z){
    dim(z) <- dim(aa)
    sum(diag(z))/sum(z)
})


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
k <- 5
m <- 20
g <- rep(LETTERS[1:k], each=20)
n <- length(g)
cf <- matrix(rnorm(2*m), 2, m)
bp <- array(NA, c(k, m), list(LETTERS[1:k], paste0("Sp", 1:m)))
bpcf <- bp
for (j in 1:m) {
    bp[,j] <- rbinom(k, 1, rbeta(1,10,15))
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
    p <- t(apply(ynew, 1, function(z) z / sum(z)))
    t(apply(p, 1, function(y)
        colMeans(apply(bp, 2, function(z) y * z * I))))
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

cl <- makeCluster(4)
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
    u <- uncertainty(o, type="multi", B=100, cl=cl)
    cal <- wa(o, y_new)
    calu <- wa(u, y_new)
    if (!any(is.na(cal)))
        gnew0[i] <- colnames(cal)[which.max(cal)]
    pm0[i,] <- cal[1,]
    if (!any(is.na(cal)))
        gnew0u[i] <- colnames(calu)[which.max(calu)]
    pm0u[i,] <- calu[1,]
}
stopCluster(cl)
iii <- gnew0 != ""
mcm(X[iii,s_col], factor(gnew0[iii], levels(X[,s_col])))
mcm(X[iii,s_col], factor(gnew0u[iii], levels(X[,s_col])))
