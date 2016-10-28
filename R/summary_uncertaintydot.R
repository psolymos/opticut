.summary_uncertainty <- function(x, sort=TRUE)
{
    uct <- x$uctab
    if (sort)
        uct <- uct[order(uct$split, -uct$R, -uct$I),]
    uct
}
