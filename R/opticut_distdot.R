.opticut_dist <- function(x) {
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
    if (is.character(x))
        return(strsplit(x, ":", fixed=TRUE)[[1L]][1L] %in% List)
    FALSE
}
