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

print.opticut1 <- function(x,
cut=getOption("ocoptions")$cut, sort=getOption("ocoptions")$sort,
digits, ...)
{
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

summary.opticut <- function(object,
cut=getOption("ocoptions")$cut, sort=getOption("ocoptions")$sort, ...)
{
    spp <- lapply(object$species, function(z)
        as.matrix(z[order(z$w, decreasing=TRUE)[1L],]))
    sppmat <- t(sapply(spp, function(z) as.matrix(z)))
    hab <- sapply(spp, rownames)
    #hab <- factor(hab, levels=unique(hab))
    colnames(sppmat) <- colnames(object$species[[1L]])
    res <- data.frame(split=hab, sppmat)
    res$assoc <- .parseAssoc(res)
    res$logL <- NULL
    bp <- bestpart(object)
    bp <- mefa4::nonDuplicated(bp, rownames(bp), TRUE)
    sgn <- sign(c(-3, -2, -1, 0, 1, 2, 3)[as.integer(res$assoc)])
    lab1 <- character(ncol(bp))
    lab0 <- character(ncol(bp))
    for (i in seq_len(ncol(bp))) {
        if (sgn[i] < 0)
            bp[,i] <- 1 - bp[,i]
        lab1[i] <- paste(rownames(bp)[bp[,i] == 1],
            collapse=getOption("ocoptions")$collapse)
        lab0[i] <- paste(rownames(bp)[bp[,i] == 0],
            collapse=getOption("ocoptions")$collapse)
    }
    oo <- if (sort)
        order(-rowSums(bp), rownames(bp)) else order(rownames(bp))
    bp <- bp[oo,,drop=FALSE]
    o <- order(colSums(bp), lab1, res$I, decreasing=TRUE)
    if (sort) {
        res <- res[o,,drop=FALSE]
        lab1 <- lab1[o]
        lab0 <- lab0[o]
        bp <- bp[,o,drop=FALSE]
    }
    keep <- res$logLR >= min(max(res$logLR), cut)
    object$summary <- res[keep,,drop=FALSE]
    object$highmat <- t(bp[,keep,drop=FALSE])
    object$highlabel <- lab1[keep]
    object$lowlabel <- lab0[keep]
    object$nsplit <- if (is.factor(object$strata))
        nlevels(object$strata) else ncol(object$strata)
    object$missing <- length(object$species) - nrow(res)
    class(object) <- c("summary.opticut")
    object
}

## generic for plotting model weights
wplot <- function (x, ...)
    UseMethod("wplot")

## plotting model weights, single species
wplot.opticut1 <-
function(x, cut=getOption("ocoptions")$cut, ylim=c(-1,1), las=1,
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
    op <- par(las=las)
    on.exit(par(op))
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

## plotting model weights, multi species
wplot.opticut <-
function(x, which=NULL,
cut=getOption("ocoptions")$cut, sort=getOption("ocoptions")$sort, las=1,
ylab="Model weight * Association", xlab="Partitions", ...)
{
    if (!is.null(which) && length(which) == 1L) {
        wplot.opticut1(x$species[[which]],
            cut=cut, ylab=ylab, xlab=xlab, ...)
    } else {
        if (is.na(x$comb))
            stop("Plot single species (user defined strata matrix, comb=NA):\nspecify 'which' argument")
        if (x$comb == "rank")
            stop("Plot single species (comb='rank'):\nspecify 'which' argument")
        if (!is.null(which) && length(which) > 1L)
            x$species <- x$species[which]
        COL <- c(colorRampPalette(c("red","yellow"))(10),
            colorRampPalette(c("yellow","green"))(10))
        br <- seq(-1, 1, 0.1)
        xx <- summary(x, cut=cut, sort=sort)

        if (nrow(xx$summary) < 2) {
            wplot.opticut1(x$species[[rownames(xx$summary)]],
                cut=cut, ylab=ylab, xlab=xlab, ...)
        } else {
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
            on.exit(par(op))
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
            invisible(ww)
        }
    }
}


plot.opticut <-
function(x, show_I=TRUE, show_S=TRUE, hr=TRUE,
cut=getOption("ocoptions")$cut, sort=getOption("ocoptions")$sort, las=1,
ylab="Relative abundance", xlab="Partitions",
palette=colorRampPalette(c("blue","green","red")),
mar=c(5, 4, 4, 4) + 0.1, ...)
{
    xx <- summary(x, cut=cut, sort=sort)
    bp <- xx$highmat
    iv <- xx$summary[,c("I", "mu0", "mu1")]
    iv$h0 <- pmin(iv$mu0, iv$mu1) / (iv$mu0 + iv$mu1)
    iv$h0[iv$h0 < 0.01] <- 0.01
    iv$h1 <- pmax(iv$mu0, iv$mu1) / (iv$mu0 + iv$mu1)
    n <- nrow(bp)
    p <- ncol(bp)
    op <- par(las=las, mar=mar)
    on.exit(par(op))
    plot(0, xlim=c(0, p), ylim=c(n, 0),
        type="n", axes=FALSE, ann=FALSE, ...)
    title(ylab=ylab, xlab=xlab)
    axis(1, at=seq_len(ncol(bp))-0.5,
        labels=colnames(bp), tick=TRUE)
    axis(2, at=seq_len(n)-0.5,
        labels=rownames(bp), tick=TRUE)
    if (show_S)
        axis(3, at=seq_len(ncol(bp))-0.5,
            labels=colSums(bp), tick=FALSE)
    if (show_I)
        axis(4, at=seq_len(n)-0.5,
            labels=format(round(iv$I, 2), nsmall = 2), tick=FALSE)
    Cols <- palette(101)
    Br <- seq(0, 1, length=101)
    if (hr)
        abline(h=1:n-0.5, col=Cols[1L])
    #abline(v=0:p, col=Cols[1L])
    for (i in seq_len(n)) {
        for (j in seq_len(p)) {
            h <- if (bp[i,j] == 1)
                iv$h1[i] else iv$h0[i]
            #Col <- grey(1-(h/1.2))
            #Col <- grey(1-h)
            Col <- Cols[which.min(h >= Br)[1L]]
            polygon(c(0,1,1,0)+j-1, 0.45*c(-h,-h,h,h)+i-0.5,
                col=Col, border=NA)
        }
    }
    box(col="grey")
    invisible(x)
}
