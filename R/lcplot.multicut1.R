lcplot.multicut1 <-
function(x, ...)
{
    lc <- .lc_cut1(x$mu, x$n, fix_fitted=getOption("ocoptions")$fix_fitted)
    plot(lc, ...)
    invisible(lc)
}
