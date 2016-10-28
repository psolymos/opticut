as.data.frame.summary.uncertainty <-
function(x, row.names = NULL, optional = FALSE, sort, ...)
{
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    as.data.frame(.summary_uncertainty(x, sort=sort),
        row.names=row.names, optional=optional, ...)
}
