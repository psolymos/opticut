as.data.frame.uncertainty <-
function(x, row.names = NULL, optional = FALSE, sort, ...)
{
    as.data.frame(as.matrix(x, sort=sort),
        row.names=row.names, optional=optional, ...)
}
