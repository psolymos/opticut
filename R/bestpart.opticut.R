## this returns an n x #spp matrix
bestpart.opticut <-
function (object, ...)
{
    out <- list()
    if (!is.na(object$comb) && object$comb == "rank") {
        for (spp in names(object$species)) {
            obj <- object$species[[spp]]
            i <- rownames(obj)[which.max(obj$logLR)]
            ## collapse value is taken from object
            ## so that post-hoc changes are not in effect
            ## avoid regexp
            out[[spp]] <- ifelse(object$strata %in%
                strsplit(i, object$collapse, fixed=TRUE)[[1L]], 1L, 0L)
        }
        out <- do.call(cbind, out)
        rownames(out) <- object$strata
    } else {
        for (spp in names(object$species)) {
            obj <- object$species[[spp]]
            i <- rownames(obj)[which.max(obj$logLR)]
            out[[spp]] <- object$strata[,i]
        }
        out <- do.call(cbind, out)
    }
    out
}
