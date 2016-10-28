as.data.frame.uncertainty <-
function(x, row.names = NULL, optional = FALSE, sort, ...)
{
    as.data.frame(summary(x), row.names=row.names, optional=optional,
        sort=sort, ...)
}
