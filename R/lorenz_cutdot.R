.lorenz_cut <-
function(object, type=c("max", "lc", "si"), force=FALSE)
{
    type <- match.arg(type)
    fit <- fitted(object)
    tmp <- lapply(seq_len(ncol(fit)), function(i) .lorenz_cut1(fit[,i], g=strata(object),
         fix_fitted=getOption("ocoptions")$fix_fitted))
    bp_max <- sapply(tmp, function(z) z[,"bp_max"])
    out <- switch(type,
        "max"=bp_max,
        "lc"=sapply(tmp, function(z) z[,"bp_lc"]),
        "si"=sapply(tmp, function(z) z[,"bp_si"]))
    colnames(out) <- colnames(fit)
    if (type != "max" && force)
        for (j in seq_len(ncol(fit)))
            out[bp_max[,j] > 0, j] <- 1
    out
}
