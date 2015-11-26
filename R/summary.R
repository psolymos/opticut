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

