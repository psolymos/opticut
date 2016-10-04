## higher than kmax is complement,
## e.g. 100 is same as 011 for our purposes
## this returns a 'contrast' matrix corresponding to
## all possible binary partitions of the factor levels n
kComb <-
function(k)
{
    k <- as.integer(k)
    if (k < 2)
        stop("k must be at least 2")
    kmax <- floor(k/2)
    s <- seq_len(k)
    clist <- lapply(seq_len(kmax), function(kk) combn(k, kk))
    ## if kmax is even, take care of cases like
    ## 1100 and 0011
    if (kmax == k/2) {
        COL <- seq_len(ncol(clist[[kmax]])/2)
        clist[[kmax]] <- clist[[kmax]][,COL, drop=FALSE]
    }
    m <- sapply(clist, ncol)
    out <- matrix(0L, k, sum(m))
    z <- 1
    for (i in seq_len(length(clist))) {
        for (j in seq_len(m[i])) {
            out[s %in% clist[[i]][,j],z] <- 1L
            z <- z + 1
        }
    }
    out
}
