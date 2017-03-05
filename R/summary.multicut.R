summary.multicut <- function(object, ...)
{

    bp <- .lorenz_cut(object, type="lc", force=FALSE)
    lab1 <- character(ncol(bp))
    lab0 <- character(ncol(bp))
    for (i in seq_len(ncol(bp))) {
        lab1[i] <- paste(rownames(bp)[bp[,i] == 1],
            collapse=getOption("ocoptions")$collapse)
        lab0[i] <- paste(rownames(bp)[bp[,i] == 0],
            collapse=getOption("ocoptions")$collapse)
    }
    logLR <- sapply(object$species, "[[", "logLR")
    fit <- fitted(object)
    lc <- t(apply(fit, 2, function(z) summary(lorenz(z))))
    res <- data.frame(
        split=lab1,
        assoc=.parseAssoc(data.frame(logLR=logLR, assoc=1)),
        lc[,c("J", "G")],
        null=sapply(object$species, "[[", "null"),
        logLR=logLR,
        logL=sapply(object$species, "[[", "logL"),
        logL_null=sapply(object$species, function(z) attr(z, "logL_null")))
    bp <- t(bp)
    attr(bp, "col.order") <- order(-colSums(bp), colnames(bp))
    attr(bp, "row.order") <- order(ncol(bp) - rowSums(bp),
        lab1, 1 - ifelse(is.na(res$G), 0, res$G), decreasing=FALSE)
    res$lablo <- lab0
    res$labhi <- lab1
    object$summary <- res
    object$bestpart <- bp
    object$mu <- t(sapply(object$species, "[[", "mu"))
    object$species <- NULL
    class(object) <- c("summary.multicut")
    object
}
