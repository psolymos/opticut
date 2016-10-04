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
        if (nrow(B) != length(object$Y))
            stop("rows in B must match the length of observations")
        if (ncol(B) < 2)
            stop("Are you kidding? ncol(B) must be > 0")
        niter <- ncol(B)
    }
    ## subset
    if (is.null(which))
        which <- names(object$species)
    spp <- names(object$species)
    names(spp) <- spp
    spp <- spp[which]
    ## subset object according to which
    object$species <- object$species[spp]

    ## template for return value
    out <- summary(object)
    out$B <- niter
    out$type <- type
    out$Y <- out$Y[,spp,drop=FALSE]
    class(out) <- "uncertainty"

    ## sequential
    if (is.null(cl)) {
        ## show progress bar by species
        if (length(spp) > 1L && interactive()) {
            res <- pbapply::pblapply(spp, function(i, ...)
                .uncertaintyOpticut1(object=object, i, type=type, B=B,
                    pb = FALSE, ...), ...)
        ## show progress bar when interactive (by B iterations)
        } else {
            res <- lapply(spp, function(i, ...)
                .uncertaintyOpticut1(object=object, i, type=type, B=B,
                    pb = interactive(), ...), ...)
        }
    ## parallel (not showing progress)
    } else {
        ## snow type cluster
        if (inherits(cl, "cluster")) {
            if (length(cl) < 2)
                stop("Are you kidding? Set cl to utilize at least 2 workers.")
            parallel::clusterEvalQ(cl, library(opticut))
            e <- new.env()
            assign("object", object, envir=e)
            assign("type", type, envir=e)
            assign("B", B, envir=e)
            parallel::clusterExport(cl, c("object","type","B"), envir=e)
            res <- parallel::parLapply(cl, spp, function(i, ...)
                .uncertaintyOpticut1(object=object, i, type=type, B=B,
                    pb = FALSE, ...), ...)
            parallel::clusterEvalQ(cl, rm(list=c("object","type","B")))
            parallel::clusterEvalQ(cl, detach(package:opticut))
        ## forking
        } else {
            if (.Platform$OS.type == "windows" && cl != 1)
                stop("Did you know that forking (cl > 1) does not work on Windows?",
                     "Try cl as a cluster instead, see ?makeCluster.")
            if (cl < 2)
                stop("Are you kidding? Set cl to utilize at least 2 workers.")
            res <- parallel::mclapply(spp, function(i, ...)
                .uncertaintyOpticut1(object=object, i, type=type, B=B,
                    pb = FALSE, ...), mc.cores = cl, ...)
        }
    }
    out$uncertainty <- res
    out
}
