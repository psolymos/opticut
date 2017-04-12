## LOO

## generic functions
loso <- function (object, ...)
    UseMethod("loso")
loto <- function (object, ...)
    UseMethod("loto")
lotso <- function (object, ...)
    UseMethod("lotso")

## refit model with -i data and IP for i
## i can be a vector for k-fold XV
.loso1 <- function(i, object, ...)
{
    pbo <- pboptions(type="none")
    on.exit(pboptions(pbo))
    ivec <- seq_len(nobs(object))
    leave_out <- ivec %in% i
    sset <- which(!leave_out)
    ## cl must be NULL for parallel
    o <- update(object, sset=sset, cl=NULL)
    Ynew <- object$Y[leave_out,,drop=FALSE]
    Xnew <- object$X[leave_out,,drop=FALSE]
    if (ncol(object$X) < 2L)
        Xnew <- NULL
    ip <- ipredict(o,
        ynew=Ynew, xnew=Xnew,
        method="analytic", cl=NULL, ...)
}

## internal function for leave-one-site-out
## with true LOO or without
## fold=NULL is leave one out (n-fold XV)
## when fold is a number (<n) it is used for k-fold XV
## !!! allow folds to be defined as input vector (is.list() or length>1)
.loso <-
function(object, fold=NULL, refit=FALSE, cl=NULL, ...)
{
    g0 <- strata(object)
    if (refit) {
        ## LOO used for IP
        Nobs <- nobs(object)
        if (is.null(fold)) {
            ii <- as.list(seq_len(Nobs))
        } else {
            if (is.list(fold) || length(fold) > 1L) {
                ii <- if (is.list(fold))
                    fold else lapply(unique(fold), function(z) which(fold == z))
            } else {
                if (fold < 2 || fold > Nobs)
                    stop("fold must be within 2 and nobs(object)")
                REP <- sample(rep(seq_len(fold), ceiling(Nobs/fold))[seq_len(Nobs)])
                ii <- lapply(seq_len(fold), function(z) which(REP == z))
            }
        }
        iplist <- pblapply(ii,
            .loso1, object=object, cl=cl, ...)
        ip <- iplist[[1L]]
        ip$ynew <- object$Y
        ip$xnew <- object$X
        ip$gnew <- g0
        ip$gnew[] <- NA
        ip$results$loglik_species <- array(0,
            c(dim(object$Y), nlevels(strata(object))))
        dimnames(ip$results$loglik_species) <- c(dimnames(object$Y),
            list(levels(g0)))
        ip$results$loglik <- array(0,
            c(nrow(object$Y), nlevels(strata(object))))
        dimnames(ip$results$loglik) <- list(rownames(object$Y),
            levels(g0))
        for (i in seq_len(length(ii))) {
            s <- ii[[i]]
            ip$results$loglik_species[s,,] <-
                iplist[[i]]$results$loglik_species[seq_len(length(s)),,]
            ip$results$loglik[s,] <-
                iplist[[i]]$results$loglik[seq_len(length(s)),]
            ip$gnew[s] <- iplist[[i]]$gnew
        }
    } else {
        ## no LOO used for IP
        ii <- NULL
        ip <- ipredict(object, ynew=object$Y,
            xnew=if (ncol(object$X) < 2L) NULL else object$X,
            method="analytic", cl=NULL, ...)
    }
    ip$strata <- g0
    ip$S <- ncol(object$Y)
    ip$multiclass <- multiclass(g0, ip$gnew)
    ip$kappa <- test_table(ip$multiclass$ctable)
    #ip$refit <- refit
    ip$folds <- ii # is.null(ii) means refit
    ip
}

## methods
loso.opticut <- function(object, fold=NULL, cl=NULL, ...)
    .loso(object, fold=fold, refit=TRUE, cl=cl, ...)
loso.multicut <- function(object, fold=NULL, cl=NULL, ...)
    .loso(object, fold=fold, refit=TRUE, cl=cl, ..., refit=TRUE)

## internal function to calculate metrics for -j (-j can be a vector)
## returns gnew integer vector
.loto1 <- function(j, ip) {
    ll <- ip$results$loglik_species[,-j,,drop=FALSE]
    LEV <- dimnames(ll)[[3L]]
    lls <- ip$results$loglik
    for (i in LEV) {
        lls[,i] <- rowSums(ll[,,i,drop=FALSE])
    }
    apply(lls, 1, which.max)
}

## internal function for calculating leave-one-taxon-out metrics
.loto <- function(object, fold=NULL, refit=FALSE, cl=NULL, ...)
{
    ip <- .loso(object, fold=fold, refit=refit, cl=cl, ...)
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
    .loto(object, refit=FALSE, cl=cl, ...)
loto.multicut <- function(object, cl=NULL, ...)
    .loto(object, refit=FALSE, cl=cl, ...)

## methods: leave-one-taxon-and-species-out

lotso.opticut <- function(object, fold=NULL, cl=NULL, ...)
    .loto(object, fold=fold, refit=TRUE, cl=cl, ...)
lotso.multicut <- function(object, fold=NULL, cl=NULL, ...)
    .loto(object, fold=fold, refit=TRUE, cl=cl, ...)

