## plotting model weights, single species
wplot.opticut1 <-
function(x, cut, ylim=c(-1,1), las=1,
ylab="Model weight * Association", xlab="Partitions",
theme, mar=c(5, 4, 4, 4) + 0.1, bty="o", ...)
{
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    w <- x$w * x$assoc
    names(w) <- rownames(x)
    if (!any(x$logLR >= cut)) {
        warning("All logLR < cut: cut ignored")
    } else {
        w <- w[x$logLR >= cut]
    }
    COL <- occolors(theme)(20)
    br <- seq(-1, 1, 0.1)
    op <- par(las=las, mar=mar)
    on.exit(par(op))
    barplot(rep(0, length(w)), width=1, space=0,
        col=COL[as.integer(base::cut(w, breaks=seq(-1, 1, 0.1)))],
        ylim=ylim, xlab=xlab, ylab=ylab, ...)
    lines(rep(which.max(abs(w))-0.5, 2), c(-1,1), col="grey", lwd=2)
    barplot(w, width=1, space=0, #border=NA,
        col=COL[as.integer(base::cut(w, breaks=seq(-1, 1, 0.1)))],
        ylim=ylim, xlab="", ylab="", add=TRUE, ...)
    abline(0,0)
    box(col="grey", bty=bty)
    invisible(w)
}
