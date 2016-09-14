.onAttach <- function(libname, pkgname){
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields=c("Version", "Date"))
    packageStartupMessage(paste(pkgname, ver[1], "\t", ver[2]))
    if (is.null(getOption("ocoptions")))
        options("ocoptions"=list(
            collapse=" ",
            cut=2,
            sort=TRUE,
            theme="br",
            check_comb=TRUE,
            try_error=FALSE))
    options("pboptions" = list(
        type = if (interactive()) "timer" else "none",
        char = "-",
        txt.width = 50,
        gui.width = 300,
        style = 3,
        initial = 0,
        title = "R progress bar",
        label = "",
        nout = 100L))
    invisible(NULL)
}

.onUnload <- function(libpath){
    options("ocoptions"=NULL)
    invisible(NULL)
}

