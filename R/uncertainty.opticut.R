uncertainty.opticut <-
function (object, which=NULL,
type=c("asymp", "boot", "multi"), B=99, cl=NULL, ...)
{
    ## sanity checks
    type <- match.arg(type)
    if (type == "multi" && is.na(object$comb))
        stop("Custom conbinations incompatible with type='multi':",
                "\nuse comb='rank' instead")
    if (type == "multi" && object$comb == "all")
        stop("comb='all' incompatible with type='multi':",
                "\nuse comb='rank' instead")
    if (length(B) < 2) {
        B <- as.integer(B)
        if (B < 1)
            stop("Are you kidding? B must be > 0")
        niter <- B
    } else {
        if (type == "asymp")
            stop("B must be a single integer")
        if (is.null(dim(B)))
            stop("B must be a matrix-like object")
        B <- as.matrix(B)
## this is not satisfied for jackknife
#        if (nrow(B) != NROW(object$Y))
#            stop("rows in B must match the length of observations")
        if (ncol(B) < 2)
            stop("Are you kidding? ncol(B) must be > 0")
        niter <- ncol(B)
    }
    ## subset
    object <- subset(object, subset=which)
    spp <- names(object$species)
#    names(spp) <- spp

    ## template for return value
    out <- summary(object)
    out$B <- niter
    out$type <- type
    class(out) <- c("uncertainty_opti", "uncertainty")

    if (inherits(cl, "cluster")) {
        parallel::clusterEvalQ(cl, library(opticut))
        e <- new.env()
        assign("object", object, envir=e)
        assign("type", type, envir=e)
        assign("B", B, envir=e)
        parallel::clusterExport(cl, c("object","type","B"), envir=e)
        on.exit(parallel::clusterEvalQ(cl, rm(list=c("object","type","B"))), add=TRUE)
        on.exit(parallel::clusterEvalQ(cl, detach(package:opticut)), add=TRUE)
    }
    res <- pbapply::pblapply(spp, function(i, ...)
        .uncertaintyOpticut1(object=object, i, type=type, B=B,
            pb = FALSE, ...), cl=cl, ...)

    names(res) <- spp
    out$uncertainty <- res
    out
}
