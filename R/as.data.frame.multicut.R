as.data.frame.multicut <-
function(x, row.names = NULL, optional = FALSE, cut, sort, ...)
{
    as.data.frame(summary(x), row.names=row.names, optional=optional,
        cut=cut, sort=sort, ...)
}
