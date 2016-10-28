as.matrix.summary.uncertainty <-
function(x, sort, ...)
{
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    .summary_uncertainty(x, sort=sort)
}
