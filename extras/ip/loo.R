## LOO

## generic functions
loso <- function (object, ...)
    UseMethod("loso")
loto <- function (object, ...)
    UseMethod("loto")
lotso <- function (object, ...)
    UseMethod("lotso")

## refit model with -i data and IP for i
.loso1 <- function(i, object, ...)
{
    pbo <- pboptions(type="none")
    on.exit(pboptions(pbo))
    ivec <- seq_len(nobs(object))
    sset <- which(ivec != i)
    ## cl must be NULL for parallel
    o <- update(object, sset=which(ivec != i), cl=NULL)
    Ynew <- object$Y[ivec == i,,drop=FALSE]
    Xnew <- object$X[ivec == i,,drop=FALSE]
    if (ncol(object$X) < 2L)
        Xnew <- NULL
    ip <- ipredict(o,
        ynew=Ynew, xnew=Xnew,
        method="analytic", cl=NULL, ...)
}

## internal function for leave-one-site-out
## with true LOO or without
.loso <-
function(object, cl=NULL, ..., refit=FALSE)
{
    g0 <- strata(object)
    if (refit) {
        ## LOO used for IP
        iplist <- pblapply(seq_len(nobs(object)),
            .loso1, object=object, cl=cl, ...)
        ip <- iplist[[1L]]
        ip$ynew <- object$Y
        ip$xnew <- object$X
        ip$gnew <- sapply(iplist, "[[", "gnew")
        ip$results$loglik_species <- array(0,
            c(dim(object$Y), nlevels(strata(object))))
        dimnames(ip$results$loglik_species) <- c(dimnames(object$Y),
            list(levels(g0)))
        ip$results$loglik <- array(0,
            c(nrow(object$Y), nlevels(strata(object))))
        dimnames(ip$results$loglik) <- list(rownames(object$Y),
            levels(g0))
        for (i in seq_len(length(iplist))) {
            ip$results$loglik_species[i,,] <-
                iplist[[i]]$results$loglik_species[1L,,]
            ip$results$loglik[i,] <-
                iplist[[i]]$results$loglik[1L,]
        }
    } else {
        ## no LOO used for IP
        ip <- ipredict(object, ynew=object$Y,
            xnew=if (ncol(object$X) < 2L) NULL else object$X,
            method="analytic", cl=NULL, ...)
    }
    ip$strata <- g0
    ip$S <- ncol(object$Y)
    ip$multiclass <- multiclass(g0, ip$gnew)
    ip$kappa <- test_table(ip$multiclass$ctable)
    ip$refit <- refit
    ip
}

## methods
loso.opticut <- function(object, cl=NULL, ...)
    .loso(object, cl=cl, ..., refit=TRUE)
loso.multicut <- function(object, cl=NULL, ...)
    .loso(object, cl=cl, ..., refit=TRUE)

## internal function to calculate metrics for -j
## returns gnew integer vector
.loto1 <- function(i, ip) {
    ll <- ip$results$loglik_species[,-i,,drop=FALSE]
    LEV <- dimnames(ll)[[3L]]
    lls <- ip$results$loglik
    for (i in LEV) {
        lls[,i] <- rowSums(ll[,,i,drop=FALSE])
    }
    apply(lls, 1, which.max)
}

## internal function for calculating leave-one-taxon-out metrics
.loto <- function(object, cl=NULL, ..., refit=FALSE)
{
    ip <- .loso(object, cl=cl, ..., refit=refit)
    g0 <- ip$strata
    LEV <- levels(g0)
    gnew <- sapply(seq_len(ip$S), .loto1, ip=ip)
    ct0 <- ip$multiclass$ctable
    kappa <- apply(gnew, 2, function(z)
        kappacoef(predicted=ctable(g0, factor(LEV[z], LEV)), reference=ct0))
    ip$gnew_species <- gnew
    ip$kappa_species <- kappa
    ip
}
#plot(sort(kappa["k",]),type="b");abline(h=0)

## methods
loto.opticut <- function(object, cl=NULL, ...)
    .loto(object, cl=cl, ..., refit=FALSE)
loto.multicut <- function(object, cl=NULL, ...)
    .loto(object, cl=cl, ..., refit=FALSE)

## methods: leave-one-taxon-and-species-out

lotso.opticut <- function(object, cl=NULL, ...)
    .loto(object, cl=cl, ..., refit=TRUE)
lotso.multicut <- function(object, cl=NULL, ...)
    .loto(object, cl=cl, ..., refit=TRUE)

