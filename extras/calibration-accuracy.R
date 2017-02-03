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
cmat <- function(x, y, drop=TRUE) {
    out <- c(
        tp=sum(x * y),
        fp=sum((1-x) * y),
        fn=sum(x * (1-y)),
        tn=sum((1-x) * (1-y)))
    if (drop)
        out else array(out, c(2L, 2L),
            list(data=c("t", "f"), class=c("p", "n")))
}
cmat(c(1,0,0,1,1,1,0,0,0,0), c(1,1,1,0,0,0,0,0,0,0))
cmat(c(1,0,0,1,1,1,0,0,0,0), c(1,1,1,0,0,0,0,0,0,0), drop=FALSE)

mcm <- function(x, y, beta=1) {
    levs <- if (is.factor(x))
        sort(levels(x)) else sort(unique(x))
    N <- length(x)
    if (length(y) != N)
        stop("x and y lengths must match")
    if (length(setdiff(y, levs)) > 0)
        warning("y has more levels than x")
    cm <- sapply(levs, function(z) {
        xx <- ifelse(x == z, 1, 0)
        yy <- ifelse(y == z, 1, 0)
        cmat(xx, yy)
    })
    Stat <- c(
        Acc = mean((cm["tp",] + cm["tn",]) / N),
        Err = mean((cm["fp",] + cm["fn",]) / N),
        Prc_m = sum(cm["tp",]) / sum(cm["tp",] + cm["fp",]),
        Prc_M = mean(cm["tp",] / (cm["tp",] + cm["fp",])),
        Rec_m = sum(cm["tp",]) / sum(cm["tp",] + cm["fn",]),
        Rec_M = mean(cm["tp",] / (cm["tp",] + cm["fn",])))
    Stat <- c(Stat,
        F_m = (beta^2 + 1) * Stat["Prc_m"] * Stat["Rec_m"] /
            (beta^2 * Stat["Prc_m"] + Stat["Rec_m"]),
        F_M = (beta^2 + 1) * Stat["Prc_M"] * Stat["Rec_M"] /
            (beta^2 * Stat["Prc_M"] + Stat["Rec_M"]))
    list(cmat=cm, stats=Stat,
        accuracy=(cm["tp",] + cm["tn",]) / N,
        error=(cm["fp",] + cm["fn",]) / N,
        precision=cm["tp",] / (cm["tp",] + cm["fp",]),
        recall=cm["tp",] / (cm["tp",] + cm["fn",]),
        beta=beta)
}

x <- sample(LETTERS[1:4], 20, replace=TRUE)
y <- x
for (i in 1:length(x))
    if (runif(1) > 0.1)
        y[i] <- sample(LETTERS[1:4], 1)
mcm(x, y)


## LOO

## define data

if (FALSE) {
data(dolina)
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
Y <- dolina$xtab[dolina$samp$method=="Q",]
X <- dolina$samp[dolina$samp$method=="Q",]
Y <- Y[,colSums(Y > 0) >= 20]
s_col <- "mhab"
}

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
Dist <- "binomial"

cl <- makeCluster(4)
nn <- nrow(Y)
gnew <- character(0)
pm <- matrix(NA, 0, nlevels(X[,s_col]))
for (i in 1:nn) {
    if (interactive()) {
        cat("\n<<< Run", i, "of", nn, ">>>\n\n")
    }
    ii <- seq_len(nrow(Y)) != i
    y_trn <- Y[ii,,drop=FALSE]
    y_new <- Y[!ii,,drop=FALSE]
    x_trn <- X[ii,,drop=FALSE]
    #x_new <- X[!ii,,drop=FALSE]
    o <- opticut(y_trn ~ 1, strata=x_trn[,s_col], dist=Dist, cl=cl)
    cal <- calibrate(o, y_new, n.chains=length(cl), cl=cl)
    gnew <- c(gnew, cal$gnew)
    pm <- rbind(pm, cal$pi)
}
stopCluster(cl)

mcm(X[1:nn,s_col], factor(gnew, levels(X[,s_col])))
(tt <- table(X[1:nn,s_col], factor(gnew, levels(X[,s_col]))))
sum(diag(tt))/sum(tt)

gnew0 <- gnew
pm0 <- pm

cl <- makeCluster(4)
nn <- nrow(Y)
mm <- ncol(Y)
gnew_list <- list()
pm_list <- list()
for (j in 1:mm) {
    gnew <- character(0)
    pm <- matrix(NA, 0, nlevels(X[,s_col]))
    for (i in 1:nn) {
        if (interactive()) {
            cat("<<< Spp", j, "of", mm, "--- Run", i, "of", nn, ">>>\n")
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

