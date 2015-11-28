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

print.opticut1 <- function(x, cut, sort, digits, ...)
{
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    xx <- x
    xx$assoc <- .parseAssoc(xx)
    xx <- xx[, c("assoc", "I", "mu0", "mu1", "logLR", "w")]
    if (sort)
        xx <- xx[order(xx$logLR, decreasing=TRUE),]
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
    cat("Multivariate opticut results, comb = ", x$comb, ", dist = ", x$dist,
        "\n", sep="")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
    cat(length(x$species), "species, ")
    cat(x$nsplit, ifelse(x$nsplit > 1, "binary splits\n", "binary split\n"))
    cat("\n")
    invisible(x)
}

print.summary.opticut <- function(x, cut, sort, digits, ...) {
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    xx <- x$summary[, c("split", "assoc", "I", "mu0", "mu1", "logLR", "w")]
    if (sort) {
        xx <- xx[attr(x$bestpart, "row.order"),]
    }
    xx <- xx[xx$logLR >= cut, , drop=FALSE]
    Missing <- nrow(x$summary) - nrow(xx)
    cat("Multivariate opticut results, comb = ", x$comb, ", dist = ", x$dist,
        "\n", sep="")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n\n", sep = "")
    print(format.data.frame(xx, digits=digits), ...)
    cat(x$nsplit, "binary", ifelse(x$nsplit > 1, "splits\n", "split\n"))
    if (Missing)
        cat(Missing, "species not shown\n")
    cat("\n")
    invisible(x)
}

summary.opticut <- function(object, ...)
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
    bp <- t(bp[order(rownames(bp)),,drop=FALSE])
    attr(bp, "col.order") <- order(-colSums(bp), colnames(bp))
    attr(bp, "row.order") <- order(ncol(bp) - rowSums(bp),
        lab1, 1 - res$I, decreasing=FALSE)
    res$lablo <- lab0
    res$labhi <- lab1
    object$summary <- res
    object$bestpart <- bp
    object$species <- NULL
    class(object) <- c("summary.opticut")
    object
}

print.uncertainty1 <-
function(x, ...)
{
    cat("Univariate opticut uncertainty results",
        ", type = ", attr(x, "type"), ", B = ", attr(x, "B"),
        "\n\n", sep="")
    print(summary(x), ...)
    invisible(x)
}

## summary:
## - highest % split label
## - selection freq
## - I
## - confint
## plot also uses level argument
summary.uncertainty <-
function(object, level=0.95, ...)
{
    object$uctab <- NA
    class(object) <- "summary.uncertainty"
    object
}
print.summary.uncertainty <-
function(x, ...)
{
    cat("Multivariate opticut uncertainty results",
        ", type = ", attr(x, "type"), ", B = ", attr(x, "B"),
        "\n\n", sep="")
    print(x$uctab[,c()], ...)
    invisible(x)
}

