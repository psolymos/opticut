## this returns an n x #spp matrix
bestpart.multicut <-
function (object, ...)
{
    bp <- .lc_cut(object, fix_fitted=getOption("ocoptions")$fix_fitted)
    bp[match(strata(object), rownames(bp)),,drop=FALSE]
}
