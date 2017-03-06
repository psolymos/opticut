multicut.default <-
function(Y, X, strata, dist="gaussian", sset=NULL, cl=NULL, ...)
{
    if (missing(strata))
        stop("It looks like that strata is missing")
    Y <- data.matrix(Y)
    if (is.null(colnames(Y)))
        colnames(Y) <- paste("Sp", seq_len(ncol(Y)))
    if (any(duplicated(colnames(Y)))) {
        warning("Duplicate column names found and renamed in LHS")
        colnames(Y) <- make.names(colnames(Y), unique = TRUE)
    }
    if (!all(colSums(abs(Y)) > 0)) {
        warning("Empty columns in Y were dropped")
        Y <- Y[,colSums(abs(Y)) > 0,drop=FALSE]
    }
    if (missing(X)) {
        X <- matrix(1, nrow(Y), 1L)
        rownames(X) <- rownames(Y)
        colnames(X) <- "(Intercept)"
    }

    if (any(is.na(Y)))
        stop("Y contains NA")
    if (any(is.na(X)))
        stop("X contains NA")
    if (any(is.na(strata)))
        stop("strata argument contains NA")

    if (is.ordered(strata)) {
        warning("ordering in strata ignored")
        class(strata) <- "factor"
    }
    ## factors are not coerced (level ordering remains intact)
    if (!is.factor(strata))
        strata <- as.factor(strata) # coerce to factor
    strata <- droplevels(strata) # drop unused levels
    Z <- strata

    if (!is.function(dist)) {
        dist <- .opticut_dist(dist, make_dist=TRUE)
        Dist <- strsplit(as.character(dist), ":", fixed=TRUE)[[1L]][1L]
        ## sanity check for rsf/rspf
        if (Dist %in% c("rsf", "rspf") && ncol(Y) > 1L)
            stop("'", Dist, "' is only available for single species in RHS")
    }

    if (ncol(Y) < 2L) {
        pbo <- pbapply::pboptions(type = "none")
        on.exit(pbapply::pboptions(pbo), add=TRUE)
    }
    if (inherits(cl, "cluster")) {
        parallel::clusterEvalQ(cl, library(opticut))
        e <- new.env()
        assign("dist", dist, envir=e)
        assign("X", X, envir=e)
        assign("Z", X, envir=e)
        assign("Y", Y, envir=e)
        assign("sset", sset, envir=e)
        parallel::clusterExport(cl, c("Y", "X","Z","dist"), envir=e)
        on.exit(parallel::clusterEvalQ(cl, rm(list=c("Y", "X","Z","dist"))), add=TRUE)
        on.exit(parallel::clusterEvalQ(cl, detach(package:opticut)), add=TRUE)
    }
    if (getOption("ocoptions")$try_error) {
        res <- pbapply::pblapply(seq_len(ncol(Y)), function(i, ...)
            try(multicut1(Y=Y[,i], X=X, Z=Z, dist=dist, sset=sset, ...)), cl=cl, ...)
        names(res) <- colnames(Y)
        Failed <- sapply(res, inherits, "try-error")
        failed <- names(res)[Failed]
        if (any(Failed)) {
            if (length(failed) == length(res))
                stop("Bad news: opticut failed for all species.")
            warning("Bad news: opticut failed for ", length(failed),
                " out of ", length(res), " species.")
        }
    } else {
        res <- pbapply::pblapply(seq_len(ncol(Y)), function(i, ...)
            multicut1(Y=Y[,i], X=X, Z=Z, dist=dist, sset=sset, ...), cl=cl, ...)
        names(res) <- colnames(Y)
        Failed <- logical(length(res))
        failed <- character(0)
    }

    NOBS <- if (is.null(sset))
        NROW(Y) else NROW(data.matrix(Y)[sset,,drop=FALSE])
    out <- list(call=match.call(),
        species=res[!Failed],
        X=X,
        Y=Y[,!Failed,drop=FALSE],
        strata=Z,
        nobs=NOBS,
        sset=sset,
        dist=dist,
        failed=failed)
    if (is.function(dist)) {
        attr(out$dist, "dist") <- deparse(substitute(dist))
        for (i in seq_len(length(out$species)))
            attr(out$species[[i]], "dist") <- deparse(substitute(dist))
    }
    class(out) <- "multicut"
    fit <- fitted(out)
    if (any(fit < 0)) {
        warning("Negative fitted values found for ",
            sum(colSums(fit < 0) > 0), " species.")
    }
    out
}
