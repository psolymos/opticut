.onAttach <- function(libname, pkgname){
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                    fields=c("Version", "Date"))
    packageStartupMessage(paste(pkgname, ver[1], "\t", ver[2]))
    if (is.null(getOption("ocoptions")))
        options("ocoptions"=list(collapse=" ", cut=2, sort=TRUE))
    invisible(NULL)
}

.onUnload <- function(libpath){
    options("ocoptions"=NULL)
    invisible(NULL)
}
