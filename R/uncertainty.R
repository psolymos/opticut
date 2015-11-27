## internal wrapper to do single species
.uncertaintyOpticut1 <-
function (object, which=NULL,
type=c("asymp", "boot", "multi"), B=99, ...)
{
#    type <- match.arg(type)
#    B <- as.integer(B)
#    if (B < 1)
#        stop("B must be > 0")
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
#    i <- which
    obj <- object$species[[which]]
    k <- which.max(obj$logLR)
    bm <- rownames(obj)[k]
#    spp <- names(m)
#    out <- list()
#    for (i in spp) {
        if (type == "asymp") {
#            k <- which.max($logLR)
            bm <- rownames(obj)[k]
#            m1 <- m[[i]]
            cf <- MASS::mvrnorm(B, coef(m1), vcov(m1))[,c(1L, 2L)]
            cf <- rbind(coef(m1)[c(1L, 2L)], cf)
            cf0 <- linkinv(cf[,1L])
            cf1 <- linkinv(cf[,1L] + cf[,2L])
            I <- 1 - (pmin(cf0, cf1) / pmax(cf0, cf1))
#            out[[i]] <- data.frame(best=bm, I=I, mu0=cf0, mu1=cf1)
            out <- data.frame(best=bm, I=I, mu0=cf0, mu1=cf1)
        }
        if (type == "boot") {
#            k <- which.max(obj$logLR)
            bm <- rownames(obj)[k]
            cf <- t(pbapply::pbsapply(seq_len(B), function(z) {
                .extractOpticut(object, which,
                    boot=TRUE,
                    internal=TRUE,
                    full_model=FALSE,
                    best=TRUE, ...)[[1L]]$coef[c(1L, 2L)]
            }))
#            cf <- rbind(coef(m[[i]])[c(1L, 2L)], cf)
            cf <- rbind(coef(m1)[c(1L, 2L)], cf)
            cf0 <- linkinv(cf[,1L])
            cf1 <- linkinv(cf[,1L] + cf[,2L])
            I <- 1 - (pmin(cf0, cf1) / pmax(cf0, cf1))
#            out[[i]] <- data.frame(best=bm, I=I, mu0=cf0, mu1=cf1)
            out <- data.frame(best=bm, I=I, mu0=cf0, mu1=cf1)
        }
        if (type == "multi") {
#            if (object$comb == "all")
#                stop("comb='all' incompatible with type='multi':",
#                    "\nuse comb='rank' instead")
#            k <- which.max(obj$logLR)
            bm <- character(B + 1L)
            bm[1L] <- rownames(obj)[k]
            mat <- matrix(NA, B + 1L, 3)
            colnames(mat) <- c("I", "mu0", "mu1")
            tmp <- as.numeric(obj[k, -1L])
            names(tmp) <- colnames(obj)[-1L]
            mat[1L, ] <- tmp[c("I", "mu0", "mu1")]
#            pb <- pbapply::startpb(0, B)
#            on.exit(pbapply::closepb(pb))
            for (j in seq_len(B)) {
                ## Z is factor, thus 'rank' applied
                mod <- .extractOpticut(object, which,
                    boot=TRUE,
                    internal=FALSE,
                    best=FALSE, ...)[[1L]]
                k <- which.max(mod$logLR)
                bm[j + 1L] <- rownames(mod)[k]
                tmp <- as.numeric(mod[k, -1L])
                names(tmp) <- colnames(mod)[-1L]
                mat[j + 1L, ] <- tmp[c("I", "mu0", "mu1")]
#                pbapply::setpb(pb, j)
            }
#            out[[i]] <- data.frame(best=bm, mat)
            out <- data.frame(best=bm, mat)
        }
#    }
    class(out) <- "uncertainty1"
    attr(out, "B") <- B
    attr(out, "type") <- type
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
    if (type == "multi" && object$comb == "all")
        stop("comb='all' incompatible with type='multi':",
                "\nuse comb='rank' instead")
    object$species <- object$species[spp]
    out <- summary(object)
    out$B <- B
    out$type <- type
    class(out) <- "uncertainty"

## bring in here cl stuff & .uncertaintyOpticut1
## plus use this:
    if (is.null(which))
        which <- names(object$species)
    bp <- bestpart(object)
    spp <- names(object$species)
    names(spp) <- spp
    spp <- spp[which]


    out
}
