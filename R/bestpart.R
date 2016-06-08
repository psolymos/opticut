bestpart <- function (object, ...)
    UseMethod("bestpart")

## this returns an n x #spp matrix
bestpart.opticut <-
function (object, ...)
{
    out <- list()
    if (object$comb == "rank") {
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

.extractOpticut <-
function (object, which=NULL, boot=FALSE,
internal=TRUE, best=TRUE, ...)
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
            out[[i]] <- .opticut1(
                Y=object$Y[j,i],
                X=object$X[j,],
                Z1=bp[j,i],
                dist=object$dist, ...)
        } else {
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

bestmodel <- function (object, ...)
    UseMethod("bestmodel")

bestmodel.opticut <-
function (object, which=NULL, ...)
{
    .extractOpticut(object, which,
        boot=FALSE,
        internal=TRUE,
        full_model=TRUE, ...)
}

bestmodel.optilevel <-
function (object, ...)
{
    NULL
}
