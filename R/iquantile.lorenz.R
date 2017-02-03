iquantile.lorenz <-
function(x, values, ...)
{
    structure(sapply(values, function(z) min(x[x[,"x"] >= z, "p"])),
        names=paste0(format(values, trim=TRUE,
        digits = max(2L, getOption("digits")))))
}

q <- seq(0, 1, 0.25)
qq <- quantile(x, q)
iquantile(x, qq)
