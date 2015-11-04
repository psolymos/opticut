## -------------- OptiCut -----------

## higher than kmax is complement, 
## e.g. 100 is same as 011 for our purposes
## this returns a 'contrast' matrix corresponding to
## all possible binary partitions of the factor levels n
allComb <- 
function(n) 
{
    n <- as.integer(n)
    if (n < 2)
        stop("n must be at least 2")
    kmax <- floor(n/2)
    if (kmax > n)
        stop("k must not be higher than n")
    s <- seq_len(n)
    clist <- lapply(seq_len(kmax), function(k) combn(n, k))
    ## if kmax is even, take care of cases like
    ## 1100 and 0011
    if (kmax == n/2) {
        COL <- seq_len(ncol(clist[[kmax]])/2)
        clist[[kmax]] <- clist[[kmax]][,COL, drop=FALSE]
    }
    m <- sapply(clist, ncol)
    out <- matrix(0L, n, sum(m))
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
modelComb <- 
function(x, collapse = " ")
{
    f <- droplevels(as.factor(x))
    LEVELS <- gsub("\\s", "", levels(f))
    i <- as.integer(f)
    n <- max(i)
    s <- seq_len(n)
    ac <- allComb(n)
    LABELS <- apply(ac, 2, function(z) 
        paste(LEVELS[as.logical(z)], collapse=collapse))
    out <- apply(ac, 2, function(z) z[match(i, s)])
    rownames(out) <- f
    colnames(out) <- LABELS
    out
}

## this checks a design matrix for complementary rows
## e.g. 1100 vs 0011
checkModelComb <- function(x) {
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
    attr(out, "comp") <- cbind(i=row(mat)[which(mat[upper.tri(mat)])], 
        j=col(mat)[which(mat[upper.tri(mat)])])
    attr(out, "same") <- cbind(i=row(mat)[which(mat[lower.tri(mat)])], 
        j=col(mat)[which(mat[lower.tri(mat)])])
    out
}

if (FALSE) {
allComb(2)
allComb(3)
allComb(4)
x <- sample(LETTERS[1:4], 10, TRUE)
modelComb(x, "_")
x <- sample(LETTERS[1:5], 100, TRUE)
modelComb(x, "_")
checkModelComb(modelComb(x, "_"))
## check if this is working OK
f <- function(n) {
    sum(sapply(1:(n-1), function(z) choose(n, z)))
}
x1 <- sapply(2:10, function(z) ncol(allComb(z))*2)
x2 <- sapply(2:10, f)
all(x1==x2)
}

## Z1 is NULL or a single column from Z
.opticut1 <- 
function(Y, X, Z1=NULL, 
dist="gaussian", linkinv, ...)
{
    if (missing(linkinv))
        linkinv <- is.null(Z1)
    if (is.null(Z1)) {
        XX <- as.data.frame(X)
    } else {
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
            mod <- glm(Y ~ .-1, data=XX, family=gaussian(link), ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- family(mod)$linkinv
        }
        if (dist == "poisson") {
            link <- list(...)$link
            if (is.null(link))
                link <- "log"
            mod <- glm(Y ~ .-1, data=XX, family=poisson(link), ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- family(mod)$linkinv
    #        se <- sqrt(diag(vcov(mod)))
        }
        if (dist == "binomial") {
            link <- list(...)$link
            if (is.null(link))
                link <- "logit"
            mod <- glm(Y ~ .-1, data=XX, family=binomial(link), ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- family(mod)$linkinv
    #        se <- sqrt(diag(vcov(mod)))
        }
        if (dist == "negbin") {
            require(MASS)
            mod <- glm.nb(Y ~ .-1, data=XX, ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- family(mod)$linkinv
    #        se <- sqrt(diag(vcov(mod)))
        }
        if (dist == "beta") {
            require(betareg)
            mod <- betareg(Y ~ .-1, data=XX, ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- mod$link$mean$linkinv
    #        se <- sqrt(diag(vcov(mod)))
        }
        if (dist == "zip") {
            require(pscl)
            mod <- zeroinfl(Y ~ .-1 | 1, data=XX, dist="poisson", ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- mod$linkinv
    #        se <- sqrt(diag(vcov(mod)))
        }
        if (dist == "zinb") {
            require(pscl)
            mod <- zeroinfl(Y ~ .-1 | 1, data=XX, dist="negbin", ...)
            cf <- coef(mod)
            ll <- as.numeric(logLik(mod))
            linv <- mod$linkinv
    #        se <- sqrt(diag(vcov(mod)))
        }
        if (dist == "ordered") {
            if (!is.null(list(...)$method))
                if (list(...)$method != "logistic")
                    stop("only logisic model allowed for ordered")
            require(MASS)
            Y <- as.ordered(Y)
            if (nlevels(Y) > 2) { # ordinal
                ## need to keep the intercept
                if (ncol(XX)==1) {
                    if (!is.null(list(...)$data))
                        stop("data argument should not be provided as part of ...")
                    mod <- polr(Y ~ 1, method="logistic", ...)
                } else {
                    mod <- polr(Y ~ ., data=data.frame(XX[,-1,drop=FALSE]), 
                        method="logistic", ...)
                }
                cf <- c(0, coef(mod))
    #            se <- rep(NA, length(cf)) # c(0, sqrt(diag(vcov(mod))))
                #se <- sqrt(diag(vcov(mod)))[1:length(cf)]
            } else {
                mod <- glm(Y ~ .-1, data=XX, family=binomial("logit"), ...)
                cf <- coef(mod)
    #            se <- sqrt(diag(vcov(mod)))
            }
            ll <- as.numeric(logLik(mod))
            linv <- binomial("logit")$linkinv
        }
        if (dist == "rspf") {
            library(ResourceSelection)
                m <- list(...)$m
            if (is.null(m))
                stop("'m' must be provided, see ?rspf")
            link <- list(...)$link
            if (is.null(link))
                link <- "logit"
    #        mod <- ResourceSelection:::rsf.fit(X=data.matrix(XX), Y=Y, m=m, link=link, ...)
            if (link == "log") {
                stop("rsf (rspf with log link) not implemented")
                #mod <- rsf(Y ~ .-1, data=XX, ...)
            } else {
                mod <- rspf(Y ~ ., data=XX[,-1,drop=FALSE], link=link, ...)
            }
            cf <- mod$coefficients
            ll <- as.numeric(mod$loglik)
            linv <- binomial(link)$linkinv
    #        se <- mod$std.error
        }
        if (!linkinv)
            linv <- NULL
    #    out <- list(coef=cf, se=se, logLik=ll, linkinv=linv)
        out <- list(coef=cf, logLik=ll, linkinv=linv)
    } else {
        out <- dist(Y, XX, linkinv, ...)
    }
    out
}

## Y is abundance vector
## X is model matrix for nuisance variables
## Z is design matrix for binary splits
opticut1 <- 
function(Y, X, Z, dist="gaussian", ...)
{
    #dist <- match.arg(dist)
    X <- data.matrix(X)
    if (is.null(rownames(X)))
        rownames(X) <- seq_len(nrow(X))
    Z <- data.matrix(Z)
    if (is.null(colnames(Z)))
        colnames(Z) <- paste0("split.", seq_len(ncol(Z)))
    if (!checkModelComb(Z))
        stop("complementary design variables found")
    if (length(unique(c(length(Y), nrow(X), nrow(Z)))) > 1)
        stop("dimension mismatch")
    N <- ncol(Z)
    res0 <- .opticut1(Y, X, Z1=NULL, dist=dist, ...)
    cf <- matrix(0, N, length(res0$coef)+1)
    rownames(cf) <- colnames(Z)
#    se <- cf
    ll <- numeric(N)
    names(ll) <- colnames(Z)
    for (i in seq_len(N)) {
        res <- .opticut1(Y, X, Z1=Z[,i], dist=dist, ...)
        cf[i,] <- res$coef
#        se[i,] <- res$se
        ll[i] <- res$logLik
    }
    dll <- ll - max(ll)
    w <- exp(dll) / sum(exp(dll))
    cfnull <- res0$linkinv(res0$coef[1L])
    cf0 <- res0$linkinv(cf[,1L])
    cf1 <- res0$linkinv(cf[,1L] + cf[,2L])
    h <- sign(cf[,2L])
    out <- data.frame(assoc=h, w=w, 
        null=cfnull,
        mu0=cf0, mu1=cf1, 
        logL=ll, logLR=ll-res0$logLik)
    rownames(out) <- colnames(Z)
    attr(out, "logL_null") <- res0$logLik
#    attr(out, "se_null") <- res0$se[1L]
#    attr(out, "se") <- cbind(se0=se[,1L], se1=se[,2L])
    attr(out, "H") <- sum(w^2)
    attr(out, "dist") <- if (is.function(dist))
        deparse(substitute(dist)) else dist
    class(out) <- c("opticut1", "data.frame")
    out    
}

## this is the main user interface
opticut <- 
function(formula, data, strata, dist="gaussian", cl=NULL, ...)
{
    if (dist=="rspf")
        stop("rspf is available for single species: use opticut1")

    if (missing(data)) 
        data <- parent.frame()
    mf <- match.call(expand.dots = FALSE)
    mm <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, mm)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.response(mf, "numeric")
    if (any(duplicated(colnames(Y))))
        stop("Duplicate colnames found in LHS")
    ff <- formula
    ff[[2]] <- NULL
    mt <- terms(ff, data = data)
    X <- model.matrix(mt, mf)
    if (missing(strata))
        stop("strata is missing")
    if (is.null(dim(strata))) {
        Z <- modelComb(strata)
    } else {
        Z <- strata
    }
    Y <- data.matrix(Y)

    ## sequential
    if (is.null(cl)) {
        ## show prograss bar
        if (ncol(Y) > 1 && interactive() && require(pbapply)) {
            res <- pbapply(Y, 2, function(yy, ...) 
                opticut1(Y=yy, X=X, Z=Z, dist=dist, ...))
        ## do not show prograss bar
        } else {
            res <- apply(Y, 2, function(yy, ...) 
                opticut1(Y=yy, X=X, Z=Z, dist=dist, ...))
        }
    ## parallel
    } else {
        library(parallel)
        ## snow type cluster
        if (inherits(cl, "cluster")) {
            clusterExport(cl, c("opticut1",".opticut1",
                "checkModelComb","allComb","modelComb"))
            e <- new.env()
            assign("dist", dist, envir=e)
            assign("X", X, envir=e)
            assign("Z", X, envir=e)
            clusterExport(cl, c("X","Z","dist"), envir=e)
            res <- parApply(cl, Y, 2, function(yy, ...) 
                opticut1(Y=yy, X=X, Z=Z, dist=dist, ...))
            clusterEvalQ(cl, rm(list=c("opticut1",".opticut1",
                "X","Z","dist",
                "checkModelComb","allComb","modelComb")))
        ## forking
        } else {
            if (cl < 2)
                stop("cl must be at least 2 for forking")
            res <- mclapply(1:ncol(Y), function(i, ...) 
                opticut1(Y=Y[,i], X=X, Z=Z, dist=dist, ...))
        }
    }

    out <- list(call=match.call(),
        species=res,
        strata=Z,
        dist=dist)
    class(out) <- "opticut"
    out
}

parseAssoc <- function(x) {
    LRc <- rep(1, nrow(x))
    LRc[x$logLR > 2] <- 2
    LRc[x$logLR > 8] <- 3
    Sign <- c("-","0","+")[x$assoc+2]
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
        xx <- xx[order(xx$w, decreasing=TRUE),]
    xx$assoc <- parseAssoc(xx)
    if (any(xx$logLR >= cut)) {
        SHOW <- which(xx$logLR >= cut)
        tmp <- if (length(SHOW) > 1)
            "Best supported models" else "Best supported model"
        TXT <- paste0(tmp, " with logLR >= ",
            format(cut, digits = digits), ":")
    } else {
        SHOW <- 1L
        TXT <- paste0("Best supported model:")
    }
    xx <- xx[SHOW,,drop=FALSE]
    cat("Univariate opticut results, dist = ", attr(x, "dist"),
        "\nw = ",
        format(xx[1,"w"], digits = digits),
        "; H = ", format(attr(x, "H"), digits = digits), 
        "; logL_null = ", format(attr(x, "logL_null"), digits = digits), 
        "\n\n", TXT, "\n", sep="")
    print.data.frame(xx, ...)
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
#    cs <- cumsum(w[order(w, decreasing=TRUE)])
#    r <- cs <= coverage
#    r[1] <- TRUE
#    w <- w[names(w) %in% names(r[r])]
    w <- w[x$logLR >= cut]
    COL <- c(colorRampPalette(c("red","yellow"))(10),
        colorRampPalette(c("yellow","green"))(10))
    br <- seq(-1, 1, 0.1)
    barplot(rep(0, length(w)), width=1, space=0, 
        col=COL[as.integer(cut(w, breaks=seq(-1, 1, 0.1)))],
        ylim=ylim, xlab=xlab, ylab=ylab, ...)
    lines(rep(which.max(abs(w))-0.5, 2), c(-1,1), col="grey", lwd=2)
    barplot(w, width=1, space=0, #border=NA,
        #col=rev(heat.colors(100))[floor(w*100)+1],
        #col=grey(1-abs(w)),
        col=COL[as.integer(cut(w, breaks=seq(-1, 1, 0.1)))],
        ylim=ylim, xlab="", ylab="", add=TRUE, ...)
    abline(0,0)
    #points(which.max(x$w)-0.5, 0)
    box()
    invisible(x)
}

print.opticut <- function(x, digits, ...) {
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    cat("Multivariate opticut results, dist =", x$dist, "\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
            "\n\n", sep = "")
    cat(length(x$species), "species, ")
    cat(ncol(x$strata), 
        ifelse(ncol(x$strata) > 1, "binary splits\n", "binary split\n"))
    cat("\n")
    invisible(x)
}

print.summary.opticut <- function(x, digits, ...) {
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    cat("Multivariate opticut results, dist =", x$dist, "\n")
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
    res$assoc <- parseAssoc(res)
    if (sort)
        res <- res[order(res$split, 1-res$w, decreasing=FALSE),]
    res <- res[res$logLR >= min(max(res$logLR), cut),,drop=FALSE]
    res$logL <- NULL
    object$summary <- res
    object$nsplit <- ncol(object$strata)
    object$missing <- length(object$species) - nrow(res)
    class(object) <- c("summary.opticut")
    object
}

plot.opticut <- 
function(x, what=NULL, cut=2, sort=TRUE, coverage=0.95, las=1, ...) 
{
    if (!is.null(what)) {
        plot(x$species[[what]])
    } else {
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


Lc_cut <-
function (x)
{
    if (is.factor(x))
        stop("x must be numeric")
    if (any(x < 0))
        stop("x must be non-negative")
    ## modeling after 'ineq'
    o <- order(x)
    z <- x[o]
    p <- seq_len(length(z))/sum(length(z))
    L <- cumsum(z)/sum(z)
    p <- c(0, p)
    L <- c(0, L)
    ## approximating 1st derivative for Youden index
    J <- p - L
    ## Gini index, after 'ineq'
    G <- sum(z * 1:length(z))
    G <- 2 * G/(length(z) * sum(z))
    G <- G - 1 - (1/length(z))
    ##
    m1 <- which.max(J)

    list(
        ## threshold (x corresponding to L)
        x=unname(z[m1]), 
        ## cumulative y
        L=unname(L[m1+1]), 
        ## cumulative x
        p=unname(p[m1+1]), 
        ## asymmetry index
        S=unname(L[m1+1] + p[m1+1]), 
        ## Gini index
        G=G, 
        ## Youden index
        J=max(J))
}

lorenzComb <- 
function(x, g, cut=0.5)
{
    Freq <- table(g=as.factor(g), 
        lc=ifelse(x >= Lc_cut(x)$x, 1L, 0L))
    Prob <- Freq[,"1"] / rowSums(Freq)
    ## high valued classes
    Hi <- names(Prob)[Prob >= cut]
    structure(ifelse(g %in% Hi, 1L, 0L),  x1=Hi)
}

rankComb <-
function(Y, X, g, dist="gaussian", ...) 
{
    g <- as.factor(x0)
    Z <- model.matrix(~g)
    m <- .opticut1(Y, X, Z1=Z[,-1,drop=FALSE], 
        linkinv=TRUE, dist=dist, ...)
    x <- rank(-c(m$coef[1], m$coef[1] + m$coef[2:ncol(Z)]))
    names(x) <- levels(g)
    o <- x[order(x)]
    ## all levels is H0, not needed (thus -1)
    Z1 <- matrix(0, length(g), nlevels(g)-1)
    colnames(Z1) <- 1:(nlevels(g)-1)
    for (i in seq_len(length(o)-1)) {
        Z1[g %in% names(o)[1:i], i] <- 1
        colnames(Z1)[i] <- paste(names(o)[1:i], collapse=" ")
    }
    Z1
}


