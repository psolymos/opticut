as.data.frame.opticut <-
function(x, row.names = NULL, optional = FALSE, cut, sort, ...)
{
    as.data.frame(as.data.frame(x, cut=cut, sort=sort),
        row.names=row.names, optional=optional, ...)
}
