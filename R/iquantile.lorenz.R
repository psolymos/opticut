iquantile.lorenz <-
function(x, values, type = c("L", "p"), ...)
{
    type <- match.arg(type)
    structure(sapply(values, function(z) min(x[x[,"x"] >= z, type])),
        names=paste0(format(values, trim=TRUE,
        digits = max(2L, getOption("digits")))))
}
