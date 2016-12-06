## this returns an n x #spp matrix
bestpart.opticut <-
function (object, pos_only=FALSE, ...)
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
    } else {
        for (spp in names(object$species)) {
            obj <- object$species[[spp]]
            i <- rownames(obj)[which.max(obj$logLR)]
            ## comb!=rank does not recognize assoc (+/-)
            out[[spp]] <- object$strata[,i]
            if (pos_only && obj[i,"assoc"] < 0) {
                out[[spp]] <- 1 - out[[spp]]
            }
        }
        out <- do.call(cbind, out)
    }
    rownames(out) <- strata(object)
    out
}
