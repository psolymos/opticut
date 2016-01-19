ocoptions <-
function(...)
{
    opar <- getOption("ocoptions")
    args <- list(...)
    if (length(args)) {
        if (length(args)==1 && is.list(args[[1]])) {
            npar <- args[[1]]
        } else {
            npar <- opar
            npar[match(names(args), names(npar))] <- args
        }
        options("ocoptions"=npar)
    }
    invisible(opar)
}

occolors <-
function(theme)
{
    if (missing(theme))
        theme <-getOption("ocoptions")$theme
    if (is.function(theme))
        return(theme)
    if (length(theme) == 1L) {
        theme <- match.arg(theme, c("br", "gr", "bw"))
        theme <- switch(theme,
            "br" = c("#2c7bb6", "#abd9e9", "#ffffbf",
                "#fdae61", "#d7191c"), # bu-yl-rd
            "gr" = c("#1a9641", "#a6d96a", "#ffffbf",
                "#fdae61", "#d7191c"), # gr-yl-rd
            "bw" = c("#FFFFFF", "#BFBFBF", "#808080",
                "#404040", "#000000")) # bw
    }
    colorRampPalette(theme)
}


