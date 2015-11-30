## internal wrapper to do single species
.uncertaintyOpticut1 <-
function (object, which=NULL,
type=c("asymp", "boot", "multi"), B=99, pb=FALSE, ...)
{
    if (missing(which))
        stop("specify which argument")
    if (!length(which))
        stop("which argument must have length 1")
    linkinv <- .opticut1(
        Y=object$Y[,1L],
        X=object$X,
        Z1=NULL,
        dist=object$dist, ...)$linkinv
    m1 <- .extractOpticut(object, which,
        boot=FALSE,
        internal=TRUE,
        full_model=TRUE,
        best=TRUE, ...)[[1L]]
    obj <- object$species[[which]]
    k <- which.max(obj$logLR)
    bm <- rownames(obj)[k]
    n <- length(object$Y)
    if (type == "asymp") {
        if (length(B) > 1)
            stop("Provide single integer for B.")
        niter <- B
        bm <- rownames(obj)[k]
        cf <- MASS::mvrnorm(niter, coef(m1), vcov(m1))[,c(1L, 2L)]
        cf <- rbind(coef(m1)[c(1L, 2L)], cf)
        cf0 <- linkinv(cf[,1L])
        cf1 <- linkinv(cf[,1L] + cf[,2L])
        I <- 1 - (pmin(cf0, cf1) / pmax(cf0, cf1))
        out <- data.frame(best=bm, I=I, mu0=cf0, mu1=cf1)
    } else {
        if (length(B) == 1) {
            BB <- replicate(B, sample.int(n, replace=TRUE))
            niter <- B
        } else {
            BB <- B
            niter <- ncol(B)
        }
    }
    if (type == "boot") {
        bm <- rownames(obj)[k]
        cf <- if (pb) {
            t(pbapply::pbapply(BB, 2, function(z) {
                .extractOpticut(object, which,
                    boot=z,
                    internal=TRUE,
                    full_model=FALSE,
                    best=TRUE, ...)[[1L]]$coef[c(1L, 2L)]
            }))
#            t(pbapply::pbsapply(seq_len(B), function(z) {
#                .extractOpticut(object, which,
#                    boot=TRUE,
#                    internal=TRUE,
#                    full_model=FALSE,
#                    best=TRUE, ...)[[1L]]$coef[c(1L, 2L)]
#            }))
        } else {
            t(sapply(BB, function(z) {
                .extractOpticut(object, which,
                    boot=z,
                    internal=TRUE,
                    full_model=FALSE,
                    best=TRUE, ...)[[1L]]$coef[c(1L, 2L)]
            }))
#            t(sapply(seq_len(B), function(z) {
#                .extractOpticut(object, which,
#                    boot=TRUE,
#                    internal=TRUE,
#                    full_model=FALSE,
#                    best=TRUE, ...)[[1L]]$coef[c(1L, 2L)]
#            }))
        }
        cf <- rbind(coef(m1)[c(1L, 2L)], cf)
        cf0 <- linkinv(cf[,1L])
        cf1 <- linkinv(cf[,1L] + cf[,2L])
        I <- 1 - (pmin(cf0, cf1) / pmax(cf0, cf1))
        out <- data.frame(best=bm, I=I, mu0=cf0, mu1=cf1)
    }
    if (type == "multi") {
        bm <- character(niter + 1L)
        bm[1L] <- rownames(obj)[k]
        mat <- matrix(NA, niter + 1L, 3)
        colnames(mat) <- c("I", "mu0", "mu1")
        tmp <- as.numeric(obj[k, -1L])
        names(tmp) <- colnames(obj)[-1L]
        mat[1L, ] <- tmp[c("I", "mu0", "mu1")]
        if (pb) {
            pbar <- pbapply::startpb(0, niter)
            on.exit(pbapply::closepb(pbar))
        }
        for (j in seq_len(niter)) {
            ## Z is factor, thus 'rank' applied
            mod <- .extractOpticut(object, which,
                boot=BB[,j],
                internal=FALSE,
                best=FALSE, ...)[[1L]]
#            mod <- .extractOpticut(object, which,
#                boot=TRUE,
#                internal=FALSE,
#                best=FALSE, ...)[[1L]]
            k <- which.max(mod$logLR)
            bm[j + 1L] <- rownames(mod)[k]
            tmp <- as.numeric(mod[k, -1L])
            names(tmp) <- colnames(mod)[-1L]
            mat[j + 1L, ] <- tmp[c("I", "mu0", "mu1")]
            if (pb)
                pbapply::setpb(pbar, j)
        }
        out <- data.frame(best=bm, mat)
        attr(out, "est") <- attr(obj, "est")
    }
    class(out) <- c("uncertainty1", "data.frame")
    attr(out, "B") <- niter
    attr(out, "type") <- type
    attr(out, "collapse") <- object$collapse
    out
}
#str(.uncertaintyOpticut1(oc, 1, type="asymp", B=1000))
#str(.uncertaintyOpticut1(oc, 1, type="boot", B=20))
#str(.uncertaintyOpticut1(oc, 1, type="multi", B=20))

uncertainty <- function (object, ...)
    UseMethod("uncertainty")

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
                stop("Unfortunately forking (cl > 1) does not work on Windows:",
                     "try cl as a cluster instead.")
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
