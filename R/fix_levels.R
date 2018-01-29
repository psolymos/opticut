fix_levels <-
function(x, sep="_")
{
    if (!is.factor(x))
        stop("x must be a factor")
    sep <- as.character(sep)
    if (length(sep) > 1L)
        stop("sep length must be 1")
    if (identical(getOption("ocoptions")$collapse, sep))
        stop("Nice try, but collapse option and sep argument are identical.")
    levels(x) <- gsub(getOption("ocoptions")$collapse, sep, x, fixed=TRUE)
    x
}
