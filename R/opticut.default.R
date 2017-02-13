opticut.default <-
function(Y, X, strata, dist="gaussian",
comb=c("rank", "all"), sset=NULL, cl=NULL, ...)
{
    comb <- match.arg(comb)
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

    if (is.null(dim(strata))) {
        if (nchar(getOption("ocoptions")$collapse) < 1)
            stop("nchar(getOption('ocoptions')$collapse) must be > 0")
        ## ordered treated as factor
        if (is.ordered(strata)) {
            warning("ordering in strata ignored")
            class(strata) <- "factor"
        }
        ## factors are not coerced (level ordering remains intact)
        if (!is.factor(strata))
            strata <- as.factor(strata) # coerce to factor
        strata <- droplevels(strata) # drop unused levels
        ## make syntactically valid names
        #levels(strata) <- make.names(levels(strata), unique = TRUE)
        ## make sure that collapse is not in levels
        if (any(grepl(getOption('ocoptions')$collapse, levels(strata), fixed=TRUE)))
            stop("Collapse value found in levels")
        if (comb == "rank") {
            Z <- strata
        }
        if (comb == "all") {
            Z <- allComb(strata) # matrix
        }
    } else {
        Z <- as.matrix(strata) # matrix
        if (getOption("ocoptions")$check_comb && !checkComb(Z))
            stop("Guess what! Complementary design variables found:\nuse 'checkComb'")
        #colnames(Z) <- make.names(colnames(Z), unique = TRUE)
        comb <- NA # user supplied matrix, not checked
    }

    if (!is.function(dist)) {
        Dist <- strsplit(as.character(dist), ":", fixed=TRUE)[[1L]][1L]
        Dist <- match.arg(Dist,
            c("gaussian","poisson","binomial","negbin",
            "beta","zip","zinb","ordered", "rsf", "rspf",
            "zip2", "zinb2"))
        ## sanity check for ordered/rsf/rspf
        if (Dist %in% c("ordered", "rsf", "rspf") && ncol(Y) > 1L)
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
            try(opticut1(Y=Y[,i], X=X, Z=Z, dist=dist, sset=sset, ...)), cl=cl, ...)
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
            opticut1(Y=Y[,i], X=X, Z=Z, dist=dist, sset=sset, ...), cl=cl, ...)
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
        nsplit=if (is.factor(Z)) # strata as factor implies K-1 splits
            (nlevels(Z) - 1L) else ncol(Z),
        dist=dist,
        comb=comb,
        scale=getOption("ocoptions")$scale,
        failed=failed,
        collapse=getOption("ocoptions")$collapse)
    if (is.function(dist)) {
        attr(out$dist, "dist") <- deparse(substitute(dist))
        for (i in which(!Failed)) {
            attr(res[[i]], "dist") <- deparse(substitute(dist))
        }
    }
    class(out) <- "opticut"
    out
}
