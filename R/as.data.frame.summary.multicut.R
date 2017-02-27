as.data.frame.summary.multicut <-
function(x, row.names = NULL, optional = FALSE, cut, sort, ...)
{
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    if (is.logical(sort)) {
        sort_r <- sort[1L]
        sort_c <- sort[1L]
    } else {
        sort_r <- 1 %in% sort
        sort_c <- 2 %in% sort
    }
    xx <- x$summary
    if (sort_r)
        xx <- xx[x$row.order, ]
    if (sort_c)
        xx <- xx[, x$col.order]
    xx <- xx[x$logLR >= cut, , drop = FALSE]
    as.data.frame(xx)
}
