strata.opticut <-
function (object, ...) {
    if (!is.na(object$comb) && is.factor(object$strata))
        return(object$strata)
    if (!is.na(object$comb) && !is.factor(object$strata))
        return(as.factor(rownames(object$strata)))
    apply(object$strata, 1, paste, collapse="")
}
