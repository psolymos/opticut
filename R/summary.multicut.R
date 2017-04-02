summary.multicut <-
function(object, ...)
{
    #bp <- .lc_cut(object, fix_fitted=getOption("ocoptions")$fix_fitted)
    bp <- sapply(object$species, "[[", "bestpart")
    lab1 <- character(ncol(bp))
    lab0 <- character(ncol(bp))
    for (i in seq_len(ncol(bp))) {
        lab1[i] <- paste(rownames(bp)[bp[,i] == 1],
            collapse=getOption("ocoptions")$collapse)
        lab0[i] <- paste(rownames(bp)[bp[,i] == 0],
            collapse=getOption("ocoptions")$collapse)
    }
    logLR <- sapply(object$species, "[[", "logLR")
    res <- data.frame(
        split=lab1,
        assoc=.parseAssoc(data.frame(logLR=logLR, assoc=1)),
        I=sapply(object$species, "[[", "I"),
        null=sapply(object$species, "[[", "null"),
        logLR=logLR,
        logL=sapply(object$species, "[[", "logL"),
        logL_null=sapply(object$species, function(z) attr(z, "logL_null")))
    bp <- t(bp)
    attr(bp, "col.order") <- order(-colSums(bp), colnames(bp))
    attr(bp, "row.order") <- order(ncol(bp) - rowSums(bp),
        lab1, 1 - ifelse(is.na(res$I), 0, res$I), decreasing=FALSE)
    res$lablo <- lab0
    res$labhi <- lab1
    object$summary <- res
    object$bestpart <- bp
    object$mu <- t(sapply(object$species, "[[", "mu"))
    object$species <- NULL
    class(object) <- c("summary.multicut")
    attr(object, "collapse") <- getOption("ocoptions")$collapse
    object
}
