summary.multicut <- function(object, ...)
{

    bp <- bestpart(object)
    bp <- mefa4::groupSums(bp, 1, rownames(bp))
    bp <- bp[levels(strata(object)),,drop=FALSE]
    bp <- bp / as.numeric(table(strata(object)))
    bp <- ifelse(bp > 0.5, 1, 0)
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
#    tmp <- t(apply(object$mu, 1, function(z)
#        ifelse((z-min(z)) / max(z-min(z)) > 0.5, 1, 0)))
#    object$col.order <- order(-colSums(tmp), colnames(tmp))
#    object$row.order <- order(ncol(tmp) - rowSums(tmp), decreasing=FALSE)
    class(object) <- c("summary.multicut")
    object
}
