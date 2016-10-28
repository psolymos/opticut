## plotting model weights, multi species
wplot.opticut <-
function(x, which=NULL, cut, sort, las=1,
ylab="Model weight * Association", xlab="Partitions",
theme, mar=c(5, 4, 4, 4) + 0.1, bty="o", ...)
{
    if (missing(cut))
        cut <- getOption("ocoptions")$cut
    if (missing(sort))
        sort <- getOption("ocoptions")$sort
    sort <- if (is.logical(sort))
        sort[1L] else 1 %in% sort
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
                cut=cut, las=las, ylab=ylab, xlab=xlab,
                theme=theme, mar=mar, bty=bty, ...)
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
            op <- par(las=las, mar=mar)
            on.exit(par(op))
            plot(0, xlim=c(0, nrow(ww)), ylim=c(ncol(ww),0),
                type="n", axes=FALSE, ann=FALSE, ...)
            title(ylab=ylab, xlab=xlab)
            axis(2, at=1:ncol(ww)-0.5,
                labels=colnames(ww), tick=TRUE, ...)
            axis(1, at=1:nrow(ww)-0.5,
                labels=rownames(ww), tick=TRUE, ...)
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
            box(col="grey", bty=bty)
            invisible(ww)
        }
    }
}
