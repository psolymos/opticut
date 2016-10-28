as.data.frame.summary.opticut <-
function(x, row.names = NULL, optional = FALSE, cut, sort, ...)
{
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    .summary_opticut(x, cut=cut, sort=sort)
}
