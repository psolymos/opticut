bestpart <- function (object, ...)
    UseMethod("bestpart")

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
            out[[spp]] <- ifelse(object$strata %in%
                strsplit(i, object$collapse)[[1L]], 1L, 0L)
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
        j <- boot
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

uncertainty <- function (object, ...)
    UseMethod("uncertainty")

uncertainty.opticut <-
function (object, which=NULL,
type=c("asymp", "boot", "multi"), B=99, ...)
{
    type <- match.arg(type)
    B <- as.integer(B)
    if (B < 1)
        stop("B must be > 0")
    linkinv <- .opticut1(
        Y=object$Y[,1L],
        X=object$X,
        Z1=NULL,
        dist=object$dist, ...)$linkinv
    m <- .extractOpticut(object, which,
        boot=FALSE,
        internal=TRUE,
        full_model=TRUE,
        best=TRUE, ...)
    spp <- names(m)
    out <- list()
    if (type == "asymp") {
        for (i in spp) {
            k <- which.max(object$species[[i]]$logLR)
            bm <- rownames(object$species[[i]])[k]
            m1 <- m[[i]]
            cf <- MASS::mvrnorm(B, coef(m1), vcov(m1))[,c(1L, 2L)]
            cf <- rbind(coef(m1)[c(1L, 2L)], cf)
            cf0 <- linkinv(cf[,1L])
            cf1 <- linkinv(cf[,1L] + cf[,2L])
            I <- 1 - (pmin(cf0, cf1) / pmax(cf0, cf1))
            out[[i]] <- data.frame(best=bm, I=I, mu0=cf0, mu1=cf1)
        }
    }
    if (type == "boot") {
        for (i in spp) {
            k <- which.max(object$species[[i]]$logLR)
            bm <- rownames(object$species[[i]])[k]
            cf <- t(pbapply::pbsapply(seq_len(B), function(z) {
                .extractOpticut(object, i,
                    boot=TRUE,
                    internal=TRUE,
                    full_model=FALSE,
                    best=TRUE, ...)[[1L]]$coef[c(1L, 2L)]
            }))
            cf <- rbind(coef(m[[i]])[c(1L, 2L)], cf)
            cf0 <- linkinv(cf[,1L])
            cf1 <- linkinv(cf[,1L] + cf[,2L])
            I <- 1 - (pmin(cf0, cf1) / pmax(cf0, cf1))
            out[[i]] <- data.frame(best=bm, I=I, mu0=cf0, mu1=cf1)
        }
    }
    if (type == "multi") {
        for (i in spp) {
            if (object$comb == "all")
                stop("comb='all' is no good, use 'rank' instead")
            k <- which.max(object$species[[i]]$logLR)
            bm <- character(B + 1L)
            bm[1L] <- rownames(object$species[[i]])[k]
            mat <- matrix(NA, B + 1L, 3)
            colnames(mat) <- c("I", "mu0", "mu1")
            tmp <- as.numeric(object$species[[i]][k, -1L])
            names(tmp) <- colnames(object$species[[i]])[-1L]
            mat[1L, ] <- tmp[c("I", "mu0", "mu1")]
            pb <- pbapply::startpb(0, B)
            on.exit(pbapply::closepb(pb))
            for (j in seq_len(B)) {
                ## Z is factor, thus 'rank' applied
                mod <- .extractOpticut(object, i,
                    boot=TRUE,
                    internal=FALSE,
                    best=FALSE, ...)[[1L]]
                k <- which.max(mod$logLR)
                bm[j + 1L] <- rownames(mod)[k]
                tmp <- as.numeric(mod[k, -1L])
                names(tmp) <- colnames(mod)[-1L]
                mat[j + 1L, ] <- tmp[c("I", "mu0", "mu1")]
                pbapply::setpb(pb, j)
            }
            out[[i]] <- data.frame(best=bm, mat)
        }
    }
    out
}
