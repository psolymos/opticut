## frequency in wplot.uncertainty1
## use oComb(attr(x, "est")) to get habitats
## return B x K matrix (0,1)
wplot.uncertainty1 <-
function(x, ylab, plot=TRUE, ...)
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
    if (plot) {
        sort <- getOption("ocoptions")$sort
        sort <- if (is.logical(sort))
            sort[1L] else 1 %in% sort
        f <- if (sort)
            rev(sort(colSums(mat) / nrow(x))) else colSums(mat) / nrow(x)
        barplot(f, col=rev(occolors()(length(f))),
            ylim=c(0,1), ylab=ylab, ...)
        box()
    }
    invisible(mat)
}
