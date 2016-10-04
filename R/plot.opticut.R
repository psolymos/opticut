plot.opticut <-
function(x, which = NULL, cut, sort,
las=1, ylab="Relative abundance", xlab="Strata",
show_I=TRUE, show_S=TRUE, hr=TRUE,
theme, mar=c(5, 4, 4, 4) + 0.1, ...)
{
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    if (is.logical(sort)) {
        sort_r <- sort[1L]
        sort_c <- sort[1L]
    } else {
        sort_r <- 1 %in% sort
        sort_c <- 2 %in% sort
    }
    if (!is.null(which)) {
        x$species <- x$species[which]
    }
    ss <- summary(x)
    xx <- ss$summary
    bp <- ss$bestpart
    if (sort_r) {
        bp <- bp[attr(ss$bestpart, "row.order"),,drop=FALSE]
        xx <- xx[attr(ss$bestpart, "row.order"),,drop=FALSE]
    }
    if (sort_c) {
        bp <- bp[,attr(ss$bestpart, "col.order"),drop=FALSE]
    }
    if (!any(xx$logLR >= cut)) {
        warning("All logLR < cut: cut ignored")
    } else {
        bp <- bp[xx$logLR >= cut, , drop=FALSE]
        xx <- xx[xx$logLR >= cut, , drop=FALSE]
    }
    n <- nrow(bp)
    p <- ncol(bp)
    op <- par(las=las, mar=mar)
    on.exit(par(op))
    plot(0, xlim=c(0, p), ylim=c(n, 0),
        type="n", axes=FALSE, ann=FALSE, ...)
    title(ylab=ylab, xlab=xlab)
    axis(1, at=seq_len(p)-0.5,
        labels=colnames(bp), tick=TRUE)
    axis(2, at=seq_len(n)-0.5,
        labels=rownames(bp), tick=TRUE)
    if (show_S)
        axis(3, at=seq_len(ncol(bp))-0.5,
            labels=colSums(bp), tick=FALSE)
    if (show_I)
        axis(4, at=seq_len(n)-0.5,
            labels=format(round(xx$I, 2), nsmall=2), tick=FALSE)
    Cols <- occolors(theme)(100)
    if (hr)
        abline(h=1:n-0.5, col=Cols[1L])
    for (i in seq_len(n)) {
        for (j in seq_len(p)) {
            h <- if (bp[i,j] == 1)
                0.5*xx$I[i]+0.5 else 0.5*(1-xx$I[i])
            ColID <- as.integer(pmin(pmax(1, floor(h * 100) + 1), 100))
            polygon(c(0,1,1,0)+j-1, 0.45*c(-h,-h,h,h)+i-0.5,
                col=Cols[ColID], border=NA)
        }
    }
    box(col="grey")
    invisible(xx)
}
