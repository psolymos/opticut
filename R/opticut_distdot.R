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
        x <- strsplit(x, ":", fixed=TRUE)[[1L]][1L]
        x <- match.arg(x, List)
        if (make_dist)
            return(x) else return(x %in% List)
    }
    FALSE
}
