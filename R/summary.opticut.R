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
    #bp <- mefa4::nonDuplicated(bp, rownames(bp), TRUE)
    bp <- bp[!duplicated(rownames(bp)),,drop=FALSE]
    bp <- bp[levels(strata(object)),]
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
#    bp <- t(bp[order(rownames(bp)),,drop=FALSE])
    bp <- t(bp)
    attr(bp, "col.order") <- order(-colSums(bp), colnames(bp))
    attr(bp, "row.order") <- order(ncol(bp) - rowSums(bp),
        lab1, 1 - ifelse(is.na(res$I), 0, res$I), decreasing=FALSE)
    res$lablo <- lab0
    res$labhi <- lab1
    object$summary <- res
    object$bestpart <- bp
    object$species <- NULL
    class(object) <- c("summary.opticut")
    object
}
