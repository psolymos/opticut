as.matrix.uncertainty <-
function(x, sort, ...)
{
    as.matrix(summary(x, ...), sort=sort)
}
