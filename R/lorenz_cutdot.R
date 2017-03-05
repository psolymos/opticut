.lorenz_cut <-
function(object, force=FALSE)
{
    bp <- bestpart(object)
    pt <- attr(bp, "p")
    bp <- mefa4::groupSums(bp, 1, rownames(bp))
    bp <- bp[levels(strata(object)),,drop=FALSE]
    ntot <- as.numeric(table(strata(object))[levels(strata(object))])
    bp <- bp / ntot
    out <- bp
    for (j in colnames(bp)) {
        tmp <- ifelse(bp[,j] >= pt[j], 1, 0)
        if (sum(tmp)==0 && force) {
            tmp[bp[,j] >= bp[which.max(bp[,j]),j]] <- 1
        }
        out[,j] <- tmp
    }
    out
}
