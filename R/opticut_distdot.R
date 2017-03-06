.opticut_dist <- function(x, make_dist=FALSE) {
    List <-c(
        "gaussian",
        "poisson",
        "binomial",
        "negbin",
        "beta",
        "zip",
        "zinb",
        "zip2",
        "zinb2",
        "rsf",
        "rspf")
    if (missing(x))
        return(List)
    if (is.function(x))
        return(FALSE)
    if (is.character(x)) {
        link <- strsplit(x, ":", fixed=TRUE)[[1L]][2L]
        x <- strsplit(x, ":", fixed=TRUE)[[1L]][1L]
        x <- match.arg(x, List)
        full <- if (is.na(link))
            x else paste(x, link, sep=":")
        if (make_dist)
            return(full) else return(x %in% List)
    }
    FALSE
}
