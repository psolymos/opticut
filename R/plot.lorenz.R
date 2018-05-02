plot.lorenz <-
function(x, type=c("L", "x"), tangent=NA, h=NA, v=NA, ...)
{
    type <- match.arg(type)
    xx <- if (type == "L")
        x[,c("p", "L")] else x[,c("p", "x")]
    ss <- summary(x)
    zz <- if (type == "L")
        list(x=ss["p(t)"], y=ss["L(t)"]) else list(x=ss["p(t)"], y=ss["x(t)"])
    if (type == "x" && !is.na(tangent)) {
        tangent <- NA
        warning("tangent cannot be plotted when type = 'x'")
    }
    plot(xx, type="l", xaxs = "i", yaxs = "i", ...)
    if (!is.na(tangent))
        abline(zz$y-zz$x, 1, col=tangent)
    if (!is.na(h))
        lines(c(0, zz$x), c(zz$y, zz$y), col=h)
    if (!is.na(v))
        lines(c(zz$x, zz$x), c(0, zz$y), col=v)
    invisible(x)
}
