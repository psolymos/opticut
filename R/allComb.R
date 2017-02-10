## this takes a classification vector
## with at least 2 levels
## and returns a model matrix with binary partitions
allComb <-
function(x, collapse)
{
    if (missing(collapse))
        collapse <-  getOption("ocoptions")$collapse
    f <- droplevels(as.factor(x))
    if (any(grepl(collapse, levels(f), fixed=TRUE)))
        stop("Collapse value found in levels.")
    LEVELS <- levels(f)
    i <- as.integer(f)
    n <- max(i, na.rm=TRUE)
    s <- seq_len(n)
    ac <- kComb(n)
    LABELS <- apply(ac, 2, function(z)
        paste(LEVELS[as.logical(z)], collapse=collapse))
    out <- apply(ac, 2, function(z) z[match(i, s)])
    rownames(out) <- f
    colnames(out) <- LABELS
    attr(out, "collapse") <- collapse
    attr(out, "comb") <- "all"
    #attr(out, "levels") <- LEVELS
    out
}
