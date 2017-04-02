lcplot.multicut1 <-
function(x,
ylab="Cumulative abundance", xlab="Strata",
bty="o", theme, ...)
{
    K <- length(x$mu)
    l <- lorenz(x$mu, x$n)
#    l <- .lc_cut1(x$mu, x$n,
#        fix_fitted=getOption("ocoptions")$fix_fitted, bp_only=FALSE)
#    bp <- attr(l, "bp")
    bp <- x$bestpart
    bp <- bp[match(rownames(l), names(bp))]
    bp[1] <- 0
    names(bp)[1] <- ""
    attr(l, "bp") <- bp
    COL <- occolors(theme)(100)
    h0 <- ifelse(bp == 1, 0.5*x$I+0.5, 0.5*(1-x$I))
    ColID <- as.integer(pmin(pmax(1, floor(h0 * 100) + 1), 100))
    plot(l, type="L", axes=FALSE, ann=FALSE, ...)
    #plot(l, type="L", axes=FALSE, ann=FALSE)
    att <- attr(l, "summary")
    abline(v=att["p[t]"], h=att["L[t]"], col="grey", lwd=0.5, lty=2)
    abline(0, 1, col="grey", lwd=0.5)
    for (i in seq_len(K)) {
        polygon(l[c(i, i+1, i+1, i), "p"],
            c(l[c(i, i+1), "L"], c(-1, -1)),
            col=COL[ColID[i+1]], border="white", lwd=0.5)
    }
    lines(l, col=1, lwd=2)
    #axis(2, seq(0, 1, 0.25))
    axis(2)
    axis(1, 0.5 * (l[-1,"p"]+l[-(K+1),"p"]), rownames(l)[-1L])
    axis(3)
    box(bty=bty)
    title(xlab=xlab, ylab=ylab)
    legend("topleft", pch=c("I", "G", "J"), bty="n",
        legend=paste("=", c(format(x$I, digits=2),
        format(att["G"], digits=2), format(att["J"], digits=2))))
    invisible(l)
}
