## this is the main user interface
opticut <-
function(formula, data, strata, dist="gaussian",
comb=c("rank", "all"), sset=NULL, cl=NULL, ...)
{
    comb <- match.arg(comb)
    if (missing(data))
        data <- parent.frame()
    if (missing(strata))
        stop("It looks like that strata is missing.")
    Strata <- deparse(substitute(strata))
    if (Strata %in% names(data))
        strata <- data[[Strata]]
    mf <- match.call(expand.dots = FALSE)
    mm <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, mm)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    Y <- model.response(mf, "numeric")
    Y <- data.matrix(Y)
    if (is.null(colnames(Y)))
        colnames(Y) <- paste("Species", seq_len(ncol(Y)))
    if (any(duplicated(colnames(Y))))
        warning("Duplicate column names found and renamed in LHS.")
    colnames(Y) <- make.names(colnames(Y))
    ff <- formula
    ff[[2]] <- NULL
    mt <- terms(ff, data = data)
    X <- model.matrix(mt, mf)

    if (is.null(dim(strata))) {
        if (nchar(getOption("ocoptions")$collapse) < 1)
            stop("nchar(getOption('ocoptions')$collapse) must be > 0")
        if (is.ordered(strata)) {
            warning("ordering in strata ignored")
            class(strata) <- "factor"
        }
        strata <- droplevels(as.factor(strata)) # factor
        ## make syntactically valid names
        #levels(strata) <- make.names(levels(strata), unique = TRUE)
        ## make sure that collapse is not in levels
        if (any(grepl(getOption('ocoptions')$collapse, levels(strata), fixed=TRUE)))
            stop("Collapse value found in levels.")
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
    Y <- data.matrix(Y)

    if (!is.function(dist)) {
        Dist <- strsplit(as.character(dist), ":", fixed=TRUE)[[1]][1]
        Dist <- match.arg(Dist,
            c("gaussian","poisson","binomial","negbin",
            "beta","zip","zinb","ordered", "rsf", "rspf",
            "zip2", "zinb2"))
        ## sanity check for rsf/rspf
        if (Dist %in% c("rsf", "rspf")) {
            if (ncol(Y) > 1L)
                stop("rsf/rspf is only available for single species in RHS")
            if (identical(as.character(ff[[2]]), "1"))
                stop("invalid formula, no covariates")
            factonly <- all(unique(sapply(mf, .MFclass)[-1]) %in% c("ordered", "factor"))
            if (factonly && dist == "rspf")
                stop("provide at least 1 continuous covariate for RSPF")
        }
        ## sanity check for rsf/rspf
        if (Dist == "ordered") {
            if (ncol(Y) > 1L)
                stop("ordered is only available for single species in RHS")
        }
    }

    ## sequential
    if (is.null(cl)) {
        ## show progress bar
        if (ncol(Y) < 2L) {
            pbo <- pboptions(type = "none")
            on.exit(pboptions(pbo))
        }
        if (getOption("ocoptions")$try_error) {
            res <- pbapply::pbapply(Y, 2, function(yy, ...)
                try(opticut1(Y=yy, X=X, Z=Z, dist=dist, sset=sset, ...)), ...)
        } else {
            res <- pbapply::pbapply(Y, 2, function(yy, ...)
                opticut1(Y=yy, X=X, Z=Z, dist=dist, sset=sset, ...), ...)
        }
    ## parallel
    } else {
        ## snow type cluster
        if (inherits(cl, "cluster")) {
            if (length(cl) < 2)
                stop("Are you kidding? Set cl to utilize at least 2 workers.")
            parallel::clusterEvalQ(cl, library(opticut))
            e <- new.env()
            assign("dist", dist, envir=e)
            assign("X", X, envir=e)
            assign("Z", X, envir=e)
            assign("Y", Y, envir=e)
            assign("sset", sset, envir=e)
            parallel::clusterExport(cl, c("Y", "X","Z","dist"), envir=e)
            if (getOption("ocoptions")$try_error) {
                res <- parallel::parLapply(cl, seq_len(ncol(Y)), function(i, ...)
                    try(opticut1(Y=Y[,i], X=X, Z=Z, dist=dist, sset=sset, ...)), ...)
            } else {
                res <- parallel::parLapply(cl, seq_len(ncol(Y)), function(i, ...)
                    opticut1(Y=Y[,i], X=X, Z=Z, dist=dist, sset=sset, ...), ...)
            }
            parallel::clusterEvalQ(cl, rm(list=c("Y", "X","Z","dist")))
            parallel::clusterEvalQ(cl, detach(package:opticut))
        ## forking
        } else {
            if (.Platform$OS.type == "windows" && cl != 1)
                stop("Did you know that forking (cl > 1) does not work on Windows?",
                     "Try cl as a cluster instead, see ?makeCluster.")
            if (cl < 2)
                stop("Are you kidding? Set cl to utilize at least 2 workers.")
            if (getOption("ocoptions")$try_error) {
                res <- parallel::mclapply(seq_len(ncol(Y)), function(i, ...)
                    try(opticut1(Y=Y[,i], X=X, Z=Z, dist=dist, sset=sset, ...)),
                    mc.cores = cl, ...)
            } else {
                res <- parallel::mclapply(seq_len(ncol(Y)), function(i, ...)
                    opticut1(Y=Y[,i], X=X, Z=Z, dist=dist, sset=sset, ...),
                    mc.cores = cl, ...)
            }
        }
    }
    if (getOption("ocoptions")$try_error) {
        Failed <- sapply(res, inherits, "try-error")
        failed <- names(res)[Failed]
        if (any(Failed)) {
            if (length(failed) == length(res))
                stop("Bad news: opticut failed for all species.")
            warning("Bad news: opticut failed for ", length(failed),
                " out of ", length(res), " species.")
        }
    } else {
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
        collapse=getOption("ocoptions")$collapse,
        dist=if (is.function(dist))
            deparse(substitute(dist)) else dist,
        comb=comb,
        failed=failed,
        collapse=getOption("ocoptions")$collapse)
    class(out) <- "opticut"
    out
}
