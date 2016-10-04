plot.uncertainty1 <-
function(x, ...)
{
    dI <- density(x$I)
    d0 <- density(x$mu0)
    d1 <- density(x$mu1)
    rnk <- rank(colMeans(x[,c("mu0", "mu1")]))
    op <- par(mfrow=c(1, 2), ...)
    on.exit(par(op))
    plot(dI, xlim=c(0,1), lwd=2, col=1, main="I", xlab="",
        sub=paste0("type = ", attr(x, "type")))
    plot(d0, xlim=range(d0$x, d1$x), ylim=c(0, max(d0$y, d1$y)),
         lwd=2, main="mu0, mu1", col=occolors()(2)[rnk[1]],
         xlab="", sub=paste0("B = ", attr(x, "B")))
    lines(d1, lwd=2, col=occolors()(2)[rnk[2]])
    invisible(list(I=dI, m0=d0, m1=d1))
}
