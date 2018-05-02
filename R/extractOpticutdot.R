.extractOpticut <-
function (object, which=NULL, boot=FALSE,
internal=TRUE, best=TRUE, Z=NULL, ...)
{
    if (is.null(which))
        which <- names(object$species)
    bp <- bestpart(object)
    spp <- names(object$species)
    names(spp) <- spp
    spp <- spp[which]
    n <- NROW(object$Y)
    if (is.logical(boot)) {
        j <- if (boot)
            sample.int(n, replace=TRUE) else seq_len(n)
    } else {
        j <- boot # boot can be the resampling vector
    }
    out <- list()
    for (i in spp) {
        if (internal) {
            if (!best)
                stop("use best=TRUE when internal=TRUE")
            ZZ <- if (is.null(Z))
                bp[j,i] else Z[j]
            out[[i]] <- .opticut1(
                Y=object$Y[j,i],
                X=object$X[j,],
                Z1=ZZ,
                dist=object$dist, ...)
        } else {
            if (!is.null(Z))
                stop("Z used only for best=TRUE & internal=TRUE")
            if (best) {
                out[[i]] <- opticut1(
                    Y=object$Y[j,i,drop=TRUE],
                    X=object$X[j,],
                    Z=bp[j,i,drop=FALSE],
                    dist=object$dist, ...)
            } else {
                zz <- object$strata
                if (is.null(dim(zz))) {
                    zz <- zz[j]
                } else {
                    zz <- zz[j,]
                }
                out[[i]] <- opticut1(
                    Y=object$Y[j,i,drop=TRUE],
                    X=object$X[j,],
                    Z=zz,
                    dist=object$dist, ...)
            }
        }
    }
    out
}
