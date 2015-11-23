## -------------- OptiCut -----------

## higher than kmax is complement,
## e.g. 100 is same as 011 for our purposes
## this returns a 'contrast' matrix corresponding to
## all possible binary partitions of the factor levels n
kComb <-
function(k)
{
    k <- as.integer(k)
    if (k < 2)
        stop("k must be at least 2")
    kmax <- floor(k/2)
    s <- seq_len(k)
    clist <- lapply(seq_len(kmax), function(kk) combn(k, kk))
    ## if kmax is even, take care of cases like
    ## 1100 and 0011
    if (kmax == k/2) {
        COL <- seq_len(ncol(clist[[kmax]])/2)
        clist[[kmax]] <- clist[[kmax]][,COL, drop=FALSE]
    }
    m <- sapply(clist, ncol)
    out <- matrix(0L, k, sum(m))
    z <- 1
    for (i in seq_len(length(clist))) {
        for (j in seq_len(m[i])) {
            out[s %in% clist[[i]][,j],z] <- 1L
            z <- z + 1
        }
    }
    out
}

## this takes a classification vector
## with at least 2 levels
## and returns a model matrix with binary partitions
allComb <-
function(x, collapse = " ")
{
    f <- droplevels(as.factor(x))
    LEVELS <- gsub("\\s", "", levels(f))
    i <- as.integer(f)
    n <- max(i)
    s <- seq_len(n)
    ac <- kComb(n)
    LABELS <- apply(ac, 2, function(z)
        paste(LEVELS[as.logical(z)], collapse=collapse))
    out <- apply(ac, 2, function(z) z[match(i, s)])
    rownames(out) <- f
    colnames(out) <- LABELS
    out
}

## this checks a design matrix for complementary rows
## e.g. 1100 vs 0011
checkComb <- function(x) {
    n <- NCOL(x)
    if (n < 2)
        return(TRUE)
    mat <- matrix(FALSE, n, n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            ## upper.tri
            mat[i,j] <- all(x[,i] == 1-x[,j]) # comp(lementary)
            ## lower.tri
            mat[j,i] <- all(x[,i] == x[,j]) # same
        }
    }
    out <- !any(mat)
    attr(out, "comp") <- cbind(
        i=row(mat)[upper.tri(row(mat))][which(mat[upper.tri(mat)])],
        j=col(mat)[upper.tri(col(mat))][which(mat[upper.tri(mat)])])
    attr(out, "same") <- cbind(
        i=row(mat)[lower.tri(row(mat))][which(mat[lower.tri(mat)])],
        j=col(mat)[lower.tri(col(mat))][which(mat[lower.tri(mat)])])
    out
}


## x is a named vector of ranks, referring to factor levels
## in some classification vector, 1=highest abundance.
oComb <-
function(x, collapse = " ")
{
    if (length(x) < 2L)
        stop("length of x must be >1")
    if (is.null(names(x)))
        names(x) <- seq_len(length(x))
    o <- x[order(x, decreasing = FALSE)]
    out <- diag(1L, length(o))
    out[upper.tri(out)] <- 1L
    out <- out[,-ncol(out)]
    rownames(out) <- names(o)
    colnames(out) <- seq_len(ncol(out))
    for (i in seq_len(ncol(out))) {
        colnames(out)[i] <- paste(names(x)[names(x) %in%
            rownames(out)[out[,i] > 0]],
            collapse = collapse)
    }
    #out <- 1L - out
    attr(out, "rank") <- o
    out
}

rankComb <-
function(Y, X, Z, dist="gaussian", collapse = " ", ...)
{
    if (!is.factor(Z))
        stop("Z must be a factor")
    Z0 <- model.matrix(~Z)
    m <- .opticut1(Y, X, Z1=Z0[,-1L,drop=FALSE],
        linkinv=TRUE, dist=dist, ...)
    lc <- c(m$coef[1], m$coef[1] + m$coef[2:ncol(Z0)])
    names(lc) <- levels(Z)
    x <- rank(-lc)
    oc <- oComb(x, collapse = collapse)
    out <- oc[match(Z, rownames(oc)),]
    attr(out, "est") <- m$linkinv(lc)
    out
}

## Z1 is:
## * NULL (this is used to fit model for H0, i.e. no partition)
## * a single column from Z matrix
## * or a matrix itself (design matrix w/o intercept)
.opticut1 <-
function(Y, X, Z1=NULL,
dist="gaussian", linkinv, full_model=FALSE, ...)
{
    if (missing(linkinv))
        linkinv <- is.null(Z1)
    if (is.null(Z1)) {
        XX <- as.data.frame(X)
    } else {
        if (is.null(dim(X)))
            X <- data.matrix(X)
        XX <- as.data.frame(cbind(X[,1,drop=FALSE], Z1, X[,-1,drop=FALSE]))
    }
    if (!is.function(dist)) {
        dist <- match.arg(dist,
            c("gaussian","poisson","binomial","negbin",
            "beta","zip","zinb","ordered", "rspf"))
        if (dist == "gaussian") {
            link <- list(...)$link
            if (is.null(link))
                link <- "identity"
            mod <- stats::glm(Y ~ .-1, data=XX, family=gaussian(link), ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- family(mod)$linkinv
        }
        if (dist == "poisson") {
            link <- list(...)$link
            if (is.null(link))
                link <- "log"
            mod <- stats::glm(Y ~ .-1, data=XX, family=poisson(link), ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- family(mod)$linkinv
        }
        if (dist == "binomial") {
            link <- list(...)$link
            if (is.null(link))
                link <- "logit"
            mod <- stats::glm(Y ~ .-1, data=XX, family=binomial(link), ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- family(mod)$linkinv
        }
        if (dist == "negbin") {
            mod <- MASS::glm.nb(Y ~ .-1, data=XX, ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- family(mod)$linkinv
        }
        if (dist == "beta") {
            mod <- betareg::betareg(Y ~ .-1, data=XX, ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- mod$link$mean$linkinv
        }
        if (dist == "zip") {
            mod <- pscl::zeroinfl(Y ~ .-1 | 1, data=XX, dist="poisson", ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- mod$linkinv
        }
        if (dist == "zinb") {
            mod <- pscl::zeroinfl(Y ~ .-1 | 1, data=XX, dist="negbin", ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- mod$linkinv
        }
        if (dist == "ordered") {
            if (!is.null(list(...)$method))
                if (list(...)$method != "logistic")
                    stop("only logisic model allowed for ordered")
            Y <- as.ordered(Y)
            if (nlevels(Y) > 2) { # ordinal
                ## need to keep the intercept
                if (ncol(XX)==1) {
                    if (!is.null(list(...)$data))
                        stop("data argument should not be provided as part of ...")
                    mod <- MASS::polr(Y ~ 1, method="logistic", ...)
                } else {
                    mod <- MASS::polr(Y ~ ., data=data.frame(XX[,-1,drop=FALSE]),
                        method="logistic", ...)
                }
                cf <- c(0, coef(mod))
            } else {
                mod <- stats::glm(Y ~ .-1, data=XX, family=binomial("logit"), ...)
                cf <- coef(mod)
            }
            ll <- as.numeric(logLik(mod))
            linv <- binomial("logit")$linkinv
        }
        if (dist == "rspf") {
                m <- list(...)$m
            if (is.null(m))
                stop("'m' must be provided, see ?rspf")
            link <- list(...)$link
            if (is.null(link))
                link <- "logit"
            if (link == "log") {
                stop("rsf (rspf with log link) not implemented")
            } else {
                mod <- ResourceSelection::rspf(Y ~ ., data=XX[,-1,drop=FALSE],
                    link=link, ...)
            }
            cf <- mod$coefficients
            ll <- as.numeric(mod$loglik)
            linv <- binomial(link)$linkinv
        }
        if (!linkinv)
            linv <- NULL
        out <- list(coef=cf, logLik=ll, linkinv=linv)
    } else {
        if (full_model)
            stop("custom distribution function: cannot return full model")
        out <- dist(Y, XX, linkinv, ...)
    }
    if (full_model)
        mod else out
}

## todo: if Z inherits from class optilevel,
## use that as best binary partition (check no. of levels)

## Y is abundance vector
## X is model matrix for nuisance variables
## Z is design matrix for binary splits or a factor (using rankComb)
opticut1 <-
function(Y, X, Z, dist="gaussian", ...)
{
    X <- data.matrix(X)
    if (is.null(rownames(X)))
        rownames(X) <- seq_len(nrow(X))
    if (is.factor(Z)) {
        Z <- rankComb(Y, X, Z, dist=dist, ...)
        Est <- attr(Z, "est")
        Comb <- "rank"
    } else {
        Est <- NA
        Comb <- "all"
    }
    Z <- data.matrix(Z)
    if (is.null(colnames(Z)))
        colnames(Z) <- paste0("split.", seq_len(ncol(Z)))
    if (!checkComb(Z))
        stop("complementary design variables found")
    if (length(unique(c(length(Y), nrow(X), nrow(Z)))) > 1)
        stop("dimension mismatch")
    N <- ncol(Z)
    res0 <- .opticut1(Y, X, Z1=NULL, dist=dist, ...)
    cf <- matrix(0, N, length(res0$coef)+1)
    rownames(cf) <- colnames(Z)
    ll <- numeric(N)
    names(ll) <- colnames(Z)
    for (i in seq_len(N)) {
        res <- .opticut1(Y, X, Z1=Z[,i], dist=dist, ...)
        cf[i,] <- res$coef
        ll[i] <- res$logLik
    }
    dll <- ll - max(ll)
    w <- exp(dll) / sum(exp(dll))
    cfnull <- res0$linkinv(res0$coef[1L])
    cf0 <- res0$linkinv(cf[,1L])
    cf1 <- res0$linkinv(cf[,1L] + cf[,2L])
    h <- sign(cf[,2L])
    I <- 1 - (pmin(cf0, cf1) / pmax(cf0, cf1))
    out <- data.frame(assoc=h, I=I,
        null=cfnull,
        mu0=cf0, mu1=cf1,
        logL=ll, logLR=ll-res0$logLik, w=w)
    rownames(out) <- colnames(Z)
    attr(out, "logL_null") <- res0$logLik
    attr(out, "H") <- sum(w^2)
    attr(out, "dist") <- if (is.function(dist))
        deparse(substitute(dist)) else dist
    attr(out, "comb") <- Comb
    attr(out, "est") <- Est
    class(out) <- c("opticut1", "data.frame")
    out
}

## this is the main user interface
opticut <-
function(formula, data, strata, dist="gaussian",
comb=c("rank", "all"), cl=NULL, ...)
{
    comb <- match.arg(comb)
    if (missing(data))
        data <- parent.frame()
    mf <- match.call(expand.dots = FALSE)
    mm <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, mm)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.response(mf, "numeric")
    Y <- data.matrix(Y)
    if (is.null(colnames(Y)))
        colnames(Y) <- paste("Species", seq_len(ncol(Y)))
    if (any(duplicated(colnames(Y))))
        stop("Duplicate colnames found in LHS")
    ff <- formula
    ff[[2]] <- NULL
    mt <- terms(ff, data = data)
    X <- model.matrix(mt, mf)
    if (missing(strata))
        stop("strata is missing")
    if (is.null(dim(strata))) {
        if (comb == "rank") {
            Z <- droplevels(as.factor(strata)) # factor
            levels(Z) <- gsub(" +", "", levels(Z))
        }
        if (comb == "all") {
            Z <- allComb(strata) # matrix
        }
    } else {
        Z <- as.matrix(strata) # matrix
        comb <- NA # user supplied matrix, not checked
    }
    Y <- data.matrix(Y)
    if (!is.function(dist))
        if (dist=="rspf" && ncol(Y) > 1L)
            stop("rspf is only available for single species in RHS")

    ## sequential
    if (is.null(cl)) {
        ## show progress bar
        if (ncol(Y) > 1L && interactive()) {
            res <- pbapply::pbapply(Y, 2, function(yy, ...)
                opticut1(Y=yy, X=X, Z=Z, dist=dist, ...), ...)
        ## do not show progress bar
        } else {
            res <- apply(Y, 2, function(yy, ...)
                opticut1(Y=yy, X=X, Z=Z, dist=dist, ...), ...)
        }
    ## parallel
    } else {
        ## snow type cluster
        if (inherits(cl, "cluster")) {

            #parallel::clusterExport(cl, c("opticut1",".opticut1",
            #    "checkComb","allComb","kComb","rankComb","oComb"))
            #parallel::clusterEvalQ(cl, library(opticut))
            parallel::clusterEvalQ(cl, requireNamespace("opticut"))
            e <- new.env()
            assign("dist", dist, envir=e)
            assign("X", X, envir=e)
            assign("Z", X, envir=e)
            parallel::clusterExport(cl, c("X","Z","dist"), envir=e)
            res <- parallel::parApply(cl, Y, 2, function(yy, ...)
                opticut1(Y=yy, X=X, Z=Z, dist=dist, ...), ...)
            #parallel::clusterEvalQ(cl, rm(list=c("opticut1",".opticut1",
            #    "X","Z","dist",
            #    "checkComb","allComb","kComb","rankComb","oComb")))
        ## forking
        } else {
            if (cl < 2)
                stop("cl must be at least 2 for forking")
            res <- parallel::mclapply(1:ncol(Y), function(i, ...)
                opticut1(Y=Y[,i], X=X, Z=Z, dist=dist, ...), ...)
        }
    }
    out <- list(call=match.call(),
        species=res,
        X=X,
        Y=Y,
        strata=Z,
        dist=dist,
        comb=comb)
    class(out) <- "opticut"
    out
}

.parseAssoc <- function(x) {
    LRc <- rep(1L, nrow(x))
    LRc[x$logLR > 2] <- 2L
    LRc[x$logLR > 8] <- 3L
    Sign <- c("-","0","+")[x$assoc + 2L]
    Assoc <- character(nrow(x))
    for (i in 1:length(Assoc))
        Assoc[i] <- paste0(rep(Sign[i], LRc[i]), collapse="")
    Assoc[x$assoc == 0] <- "0"
    factor(Assoc, levels=c("---","--","-","0","+","++","+++"))
}

print.opticut1 <- function(x, cut=2, sort=TRUE, digits, ...) {
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    xx <- x
    if (sort)
        xx <- xx[order(xx$I, decreasing=TRUE),]
    xx$assoc <- .parseAssoc(xx)
    if (any(xx$logLR >= cut)) {
        SHOW <- which(xx$logLR >= cut)
        tmp <- if (length(SHOW) > 1L)
            "Best supported models" else "Best supported model"
        TXT <- paste0(tmp, " with logLR >= ",
            format(cut, digits = digits), ":")
    } else {
        SHOW <- 1L
        TXT <- paste0("Best supported model:")
    }
    xx <- xx[SHOW,,drop=FALSE]
    cat("Univariate opticut results, comb = ", attr(x, "comb"),
        ", dist = ", attr(x, "dist"),
        "\nI = ", format(xx[1L,"I"], digits = digits),
        "; w = ", format(xx[1L,"w"], digits = digits),
        "; H = ", format(attr(x, "H"), digits = digits),
        "; logL_null = ", format(attr(x, "logL_null"), digits = digits),
        "\n\n", TXT, "\n", sep="")
    print.data.frame(xx, digits=digits, ...)
    DROP <- nrow(x) - nrow(xx)
    if (DROP > 0) {
        cat(nrow(x), " binary ",
            ifelse(nrow(x) > 1, "splits", "split"),
            " (", DROP,
            ifelse(DROP > 1, " models", " model"),
            " not shown)\n", sep="")
    } else {
        cat(nrow(x), "binary",
            ifelse(nrow(x) > 1, "splits", "split"), "\n")
    }
    cat("\n")
    invisible(x)
}


plot.opticut1 <-
function(x, cut=2, ylim=c(-1,1),
ylab="Model weight * Association", xlab="Partitions", ...)
{
    w <- x$w * x$assoc
    names(w) <- rownames(x)
    if (!any(x$logLR >= cut)) {
        warning("All logLR < cut: cut ignored")
    } else {
        w <- w[x$logLR >= cut]
    }
    COL <- c(colorRampPalette(c("red","yellow"))(10),
        colorRampPalette(c("yellow","green"))(10))
    br <- seq(-1, 1, 0.1)
    barplot(rep(0, length(w)), width=1, space=0,
        col=COL[as.integer(cut(w, breaks=seq(-1, 1, 0.1)))],
        ylim=ylim, xlab=xlab, ylab=ylab, ...)
    lines(rep(which.max(abs(w))-0.5, 2), c(-1,1), col="grey", lwd=2)
    barplot(w, width=1, space=0, #border=NA,
        col=COL[as.integer(cut(w, breaks=seq(-1, 1, 0.1)))],
        ylim=ylim, xlab="", ylab="", add=TRUE, ...)
    abline(0,0)
    box()
    invisible(x)
}

print.opticut <- function(x, digits, ...) {
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    cat("Multivariate opticut results, comb = ", x$comb, ", dist =", x$dist, "\n", sep="")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
    cat(length(x$species), "species, ")
    nstr <- if (is.factor(x$strata))
        nlevels(x$strata) else ncol(x$strata)
    cat(nstr, ifelse(nstr > 1, "binary splits\n", "binary split\n"))
    cat("\n")
    invisible(x)
}

print.summary.opticut <- function(x, digits, ...) {
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    cat("Multivariate opticut results, comb = ", x$comb, ", dist =", x$dist, "\n", sep="")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
    print(format.data.frame(x$summary, digits=digits))
    cat(x$nsplit, "binary", ifelse(x$nsplit > 1, "splits\n", "split\n"))
    if (x$missing)
        cat(x$missing, "species not shown\n")
    cat("\n")
    invisible(x)
}

summary.opticut <- function(object, cut=2, sort=TRUE, ...) {
    spp <- lapply(object$species, function(z)
        as.matrix(z[order(z$w, decreasing=TRUE)[1L],]))
    sppmat <- t(sapply(spp, function(z) as.matrix(z)))
    hab <- sapply(spp, rownames)
    #hab <- factor(hab, levels=unique(hab))
    colnames(sppmat) <- colnames(object$species[[1L]])
    res <- data.frame(split=hab, sppmat)
    res$assoc <- .parseAssoc(res)
    if (sort)
        res <- res[order(res$split, 1-res$I, decreasing=FALSE),]
    res <- res[res$logLR >= min(max(res$logLR), cut),,drop=FALSE]
    res$logL <- NULL
    object$summary <- res
    object$nsplit <- if (is.factor(object$strata))
        nlevels(object$strata) else ncol(object$strata)
    object$missing <- length(object$species) - nrow(res)
    class(object) <- c("summary.opticut")
    object
}

plot.opticut <-
function(x, which=NULL, cut=2, sort=TRUE, las=1,
ylab="Model weight * Association", xlab="Partitions", ...)
{
    if (!is.null(which) && length(which) == 1L) {
        plot(x$species[[which]],
            cut=cut, ylab=ylab, xlab=xlab, ...)
    } else {
        if (x$comb == "rank")
            stop("Plot single species (comb='rank'): specify 'which' argument")
        if (!is.null(which) && length(which) > 1L)
            x$species <- x$species[which]
        COL <- c(colorRampPalette(c("red","yellow"))(10),
            colorRampPalette(c("yellow","green"))(10))
        br <- seq(-1, 1, 0.1)
        xx <- summary(x, cut=cut, sort=sort)
        nsplit <- xx$nsplit
        nspp <- nrow(xx$summary)
        #spp <- xx$species[rownames(xx$summary)]
        ww <- sapply(xx$species[rownames(xx$summary)], "[[", "w")
        ss <- sapply(xx$species[rownames(xx$summary)], "[[", "assoc")
        ss[ss==0] <- 1
        ww <- ww * ss
        rownames(ww) <- rownames(xx$species[[1]])
        colnames(ww) <- rownames(xx$summary)
        llr <- sapply(xx$species[rownames(xx$summary)], "[[", "logLR")
        ww[llr < cut] <- 0
        ww <- ww[rowSums(ww) != 0,,drop=FALSE]
        op <- par(las=las)
        plot(0, xlim=c(0, nrow(ww)), ylim=c(ncol(ww),0),
            type="n", axes=FALSE, ann=FALSE, ...)
        title(ylab=ylab, xlab=xlab)
        axis(2, at=1:ncol(ww)-0.5,
            labels=colnames(ww), tick=TRUE)
        axis(1, at=1:nrow(ww)-0.5,
            labels=rownames(ww), tick=TRUE)
        abline(h=1:ncol(ww)-0.5)
        abline(v=0:nrow(ww), col="lightgrey")
        for (i in 1:ncol(ww)) {
            lines(rep(which.max(abs(ww[,i])), 2)-0.5, c(-0.5,0.5)+i-0.5,
                col="grey", lwd=2)
            for (j in 1:nrow(ww)) {
                h <- - ww[j,i] * 0.45
                polygon(c(0,1,1,0)+j-1, c(0,0,h,h)+i-0.5,
                    col=COL[as.integer(cut(-h, breaks=seq(-1, 1, 0.1)))])
                    #col=grey(1-abs(ww[j,i])))
            }
        }
        box()
        par(op)
        invisible(ww)
    }
}

bestpart <- function (object, ...)
    UseMethod("bestpart")

bestpart.opticut <-
function (object, ...)
{
    out <- list()
    if (object$comb == "rank") {
        for (spp in names(object$species)) {
            i <- rownames(object$species[[spp]])[which.max(object$species[[spp]]$logLR)]
            out[[spp]] <- ifelse(object$strata %in% strsplit(i, " ")[[1L]], 1L, 0L)
        }
        out <- do.call(cbind, out)
        rownames(out) <- object$strata
    } else {
        for (spp in names(object$species)) {
            i <- rownames(object$species[[spp]])[which.max(object$species[[spp]]$logLR)]
            out[[spp]] <- object$strata[,i]
        }
        out <- do.call(cbind, out)
    }
    out
}

.extractOpticut <-
function (object, which=NULL, boot=FALSE,
internal=TRUE, best=TRUE, ...)
{
    if (is.null(which))
        which <- names(object$species)
    bp <- bestpart(object)
    spp <- names(object$species)
    names(spp) <- spp
    spp <- spp[which]
    n <- NROW(object$Y)
    if (is.logical(boot)) {
        j <- if (boot)
            sample.int(n, replace=TRUE) else seq_len(n)
    } else {
        j <- boot
    }
    out <- list()
    for (i in spp) {
        if (internal) {
            if (!best)
                stop("use best=TRUE when internal=TRUE")
            out[[i]] <- .opticut1(
                Y=object$Y[j,i],
                X=object$X[j,],
                Z1=bp[j,i],
                dist=object$dist, ...)
        } else {
            if (best) {
                out[[i]] <- opticut1(
                    Y=object$Y[j,i,drop=TRUE],
                    X=object$X[j,],
                    Z=bp[j,i,drop=FALSE],
                    dist=object$dist, ...)
            } else {
                zz <- object$strata
                if (is.null(dim(zz))) {
                    zz <- zz[j]
                } else {
                    zz <- zz[j,]
                }
                out[[i]] <- opticut1(
                    Y=object$Y[j,i,drop=TRUE],
                    X=object$X[j,],
                    Z=zz,
                    dist=object$dist, ...)
            }
        }
    }
    out
}

bestmodel <- function (object, ...)
    UseMethod("bestmodel")

bestmodel.opticut <-
function (object, which=NULL, ...)
{
    .extractOpticut(object, which,
        boot=FALSE,
        internal=TRUE,
        full_model=TRUE, ...)
}

uncertainty <- function (object, ...)
    UseMethod("uncertainty")

uncertainty.opticut <-
function (object, which=NULL,
type=c("asymp", "boot", "multi"), B=99, ...)
{
    type <- match.arg(type)
    B <- as.integer(B)
    if (B < 1)
        stop("B must be > 0")
    linkinv <- .opticut1(
        Y=object$Y[,1L],
        X=object$X,
        Z1=NULL,
        dist=object$dist, ...)$linkinv
    m <- .extractOpticut(object, which,
        boot=FALSE,
        internal=TRUE,
        full_model=TRUE,
        best=TRUE, ...)
    spp <- names(m)
    out <- list()
    if (type == "asymp") {
        for (i in spp) {
            k <- which.max(object$species[[i]]$logLR)
            bm <- rownames(object$species[[i]])[k]
            m1 <- m[[i]]
            cf <- MASS::mvrnorm(B, coef(m1), vcov(m1))[,c(1L, 2L)]
            cf <- rbind(coef(m1)[c(1L, 2L)], cf)
            cf0 <- linkinv(cf[,1L])
            cf1 <- linkinv(cf[,1L] + cf[,2L])
            I <- 1 - (pmin(cf0, cf1) / pmax(cf0, cf1))
            out[[i]] <- data.frame(best=bm, I=I, mu0=cf0, mu1=cf1)
        }
    }
    if (type == "boot") {
        for (i in spp) {
            k <- which.max(object$species[[i]]$logLR)
            bm <- rownames(object$species[[i]])[k]
            cf <- t(pbapply::pbsapply(seq_len(B), function(z) {
                .extractOpticut(object, i,
                    boot=TRUE,
                    internal=TRUE,
                    full_model=FALSE,
                    best=TRUE, ...)[[1L]]$coef[c(1L, 2L)]
            }))
            cf <- rbind(coef(m[[i]])[c(1L, 2L)], cf)
            cf0 <- linkinv(cf[,1L])
            cf1 <- linkinv(cf[,1L] + cf[,2L])
            I <- 1 - (pmin(cf0, cf1) / pmax(cf0, cf1))
            out[[i]] <- data.frame(best=bm, I=I, mu0=cf0, mu1=cf1)
        }
    }
    if (type == "multi") {
        for (i in spp) {
            if (object$comb == "all")
                stop("comb='all' is no good, use 'rank' instead")
            k <- which.max(object$species[[i]]$logLR)
            bm <- character(B + 1L)
            bm[1L] <- rownames(object$species[[i]])[k]
            mat <- matrix(NA, B + 1L, 3)
            colnames(mat) <- c("I", "mu0", "mu1")
            tmp <- as.numeric(object$species[[i]][k, -1L])
            names(tmp) <- colnames(object$species[[i]])[-1L]
            mat[1L, ] <- tmp[c("I", "mu0", "mu1")]
            pb <- pbapply::startpb(0, B)
            on.exit(pbapply::closepb(pb))
            for (j in seq_len(B)) {
                ## Z is factor, thus 'rank' applied
                mod <- .extractOpticut(object, i,
                    boot=TRUE,
                    internal=FALSE,
                    best=FALSE, ...)[[1L]]
                k <- which.max(mod$logLR)
                bm[j + 1L] <- rownames(mod)[k]
                tmp <- as.numeric(mod[k, -1L])
                names(tmp) <- colnames(mod)[-1L]
                mat[j + 1L, ] <- tmp[c("I", "mu0", "mu1")]
                pbapply::setpb(pb, j)
            }
            out[[i]] <- data.frame(best=bm, mat)
        }
    }
    out
}

