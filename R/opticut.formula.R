opticut.formula <-
function(formula, data, strata, dist="gaussian",
comb=c("rank", "all"), sset=NULL, cl=NULL, ...)
{
    if (missing(data))
        data <- parent.frame()
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

    opticut.default(Y=Y, X=X, strata=strata, dist=dist,
        comb=comb, sset=sset, cl=cl, ...)
}
