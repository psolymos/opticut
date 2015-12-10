## generic for plotting model weights
wplot <- function (x, ...)
    UseMethod("wplot")

## plotting model weights, single species
wplot.opticut1 <-
function(x, cut, ylim=c(-1,1), las=1,
ylab="Model weight * Association", xlab="Partitions", theme, ...)
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
    if (!missing(theme) && is.character(theme))
        if (length(theme) == 1L)
            if (theme == "bw")
                #warning("'bw' theme not best suited for wplot")
                theme <- occolors(c("black", "lightgrey", "black"))
    COL <- occolors(theme)(20)
#    COL <- c(colorRampPalette(c("red","yellow"))(10),
#         colorRampPalette(c("yellow","green"))(10))
    br <- seq(-1, 1, 0.1)
    op <- par(las=las)
    on.exit(par(op))
    barplot(rep(0, length(w)), width=1, space=0,
        col=COL[as.integer(base::cut(w, breaks=seq(-1, 1, 0.1)))],
        ylim=ylim, xlab=xlab, ylab=ylab, ...)
    lines(rep(which.max(abs(w))-0.5, 2), c(-1,1), col="grey", lwd=2)
    barplot(w, width=1, space=0, #border=NA,
        col=COL[as.integer(base::cut(w, breaks=seq(-1, 1, 0.1)))],
        ylim=ylim, xlab="", ylab="", add=TRUE, ...)
    abline(0,0)
    box()
    invisible(w)
}

## plotting model weights, multi species
wplot.opticut <-
function(x, which=NULL, cut, sort, las=1,
ylab="Model weight * Association", xlab="Partitions", theme, ...)
{
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    if (!is.null(which) && length(which) == 1L) {
        wplot.opticut1(x$species[[which]],
            cut=cut, las=las, ylab=ylab, xlab=xlab, theme=theme, ...)
    } else {
        if (is.na(x$comb))
            stop("Plot single species (user defined strata matrix, comb=NA):",
                "\nspecify 'which' argument")
        if (x$comb == "rank")
            stop("Plot single species (comb='rank'):"
                ,"\nspecify 'which' argument")
        if (!is.null(which) && length(which) > 1L)
            x$species <- x$species[which]
        if (!missing(theme) && is.character(theme))
            if (length(theme) == 1L)
                if (theme == "bw")
                    #warning("'bw' theme not best suited for wplot")
                    theme <- occolors(c("black", "lightgrey", "black"))
        COL <- occolors(theme)(20)
#        COL <- c(colorRampPalette(c("red","yellow"))(10),
#            colorRampPalette(c("yellow","green"))(10))
        br <- seq(-1, 1, 0.1)
        sss <- summary(x)
        xx <- sss$summary
        if (sort) {
            xx <- xx[attr(sss$bestpart, "row.order"),]
        }
        xx <- xx[xx$logLR >= cut, , drop=FALSE]
        if (nrow(xx) < 2) {
            wplot.opticut1(x$species[[rownames(xx$summary)]],
                cut=cut, las=las, ylab=ylab, xlab=xlab, theme=theme, ...)
        } else {
            nsplit <- xx$nsplit
            nspp <- nrow(xx)
            sppnam <- rownames(xx)
            ## subsetting according to xx (sorted & cut)
            ww <- sapply(x$species[sppnam], "[[", "w")
            ss <- sapply(x$species[sppnam], "[[", "assoc")
            ss[ss==0] <- 1
            ww <- ww * ss
            rownames(ww) <- rownames(x$species[[1]])
            colnames(ww) <- rownames(xx)
            llr <- sapply(x$species[sppnam], "[[", "logLR")
            ww[llr < cut] <- 0
            ww <- ww[rowSums(ww) != 0,,drop=FALSE]
            op <- par(las=las)
            on.exit(par(op))
            plot(0, xlim=c(0, nrow(ww)), ylim=c(ncol(ww),0),
                type="n", axes=FALSE, ann=FALSE, ...)
            title(ylab=ylab, xlab=xlab)
            axis(2, at=1:ncol(ww)-0.5,
                labels=colnames(ww), tick=TRUE)
            axis(1, at=1:nrow(ww)-0.5,
                labels=rownames(ww), tick=TRUE)
            abline(h=1:ncol(ww)-0.5)
            abline(v=0:nrow(ww), col="lightgrey")
            for (i in 1:ncol(ww)) {
                lines(rep(which.max(abs(ww[,i])), 2)-0.5, c(-0.5,0.5)+i-0.5,
                    col="grey", lwd=2)
                for (j in 1:nrow(ww)) {
                    h <- - ww[j,i] * 0.45
                    polygon(c(0,1,1,0)+j-1, c(0,0,h,h)+i-0.5,
                        col=COL[as.integer(base::cut(-h, breaks=seq(-1, 1, 0.1)))])
                }
            }
            box(col="grey")
            invisible(ww)
        }
    }
}

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
    if (!is.null(which)) {
        x$species <- x$species[which]
    }
    ss <- summary(x)
    xx <- ss$summary
    bp <- ss$bestpart
    if (sort) {
        bp <- bp[attr(ss$bestpart, "row.order"),
            attr(ss$bestpart, "col.order"),drop=FALSE]
        xx <- xx[attr(ss$bestpart, "row.order"),,drop=FALSE]
    }
    if (!any(xx$logLR >= cut)) {
        warning("All logLR < cut: cut ignored")
    } else {
        bp <- bp[xx$logLR >= cut, , drop=FALSE]
        xx <- xx[xx$logLR >= cut, , drop=FALSE]
    }
    MinVal <- 0.01
    xx$h0 <- pmin(xx$mu0, xx$mu1) / (xx$mu0 + xx$mu1)
    xx$h0[xx$h0 < MinVal] <- MinVal
    xx$h1 <- pmax(xx$mu0, xx$mu1) / (xx$mu0 + xx$mu1)
    xx$h0[is.na(xx$I)] <- MinVal
    xx$h1[is.na(xx$I)] <- MinVal
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
    Cols <- occolors(theme)(101)
    Br <- seq(0, 1, length=101)
    if (hr)
        abline(h=1:n-0.5, col=Cols[1L])
    for (i in seq_len(n)) {
        for (j in seq_len(p)) {
            h <- if (bp[i,j] == 1)
                xx$h1[i] else xx$h0[i]
            Col <- Cols[which.min(h >= Br)[1L]]
            polygon(c(0,1,1,0)+j-1, 0.45*c(-h,-h,h,h)+i-0.5,
                col=Col, border=NA)
        }
    }
    box(col="grey")
    invisible(xx)
}

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

## frequency in wplot.uncertainty1
## use oComb(attr(x, "est")) to get habitats
## return B x K matrix (0,1)
wplot.uncertainty1 <-
function(x, ylab, ...)
{
    if (missing(ylab))
        ylab <- "Selection"
    cm <- oComb(attr(x, "est"))
    mat <- matrix(0L, nrow(x), nrow(cm))
    colnames(mat) <- rownames(cm)
    for (i in seq_len(nrow(x))) {
        j <- strsplit(as.character(x$best[i]),
            attr(x, "collapse"), fixed=TRUE)[[1L]]
        mat[i,j] <- 1L
    }
    f <- rev(sort(colSums(mat) / nrow(x)))
    barplot(f, col=rev(occolors()(length(f))),
        ylim=c(0,1), ylab=ylab, ...)
    box()
    invisible(mat)
}
