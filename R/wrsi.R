## wrsi: weighted relative suitability index
wrsi <-
function(y, x)
{
    Y <- as.integer(ifelse(y > 0, 1L, 0L))
    X <- data.matrix(x)
    n <- length(Y)
    if (nrow(X) != n)
        stop("dimension mismatch: X and Y")
    K <- ncol(X)
    if (is.null(colnames(X)))
        colnames(X) <- paste0("V", seq_len(K))
    ## relative suitability index
    ## # of available units of type k / total # of available units (any type)
    Pavail <- colSums(X) / sum(X)
    ## # of used units of type k / total # of used units (any type)
    Xu <- X * Y
    ## sum(Xu) = sum(Y) except when rowsum != 1
    Pused <- colSums(Xu) / sum(Xu)
    ## crude weighted p-occ
    Pw <- colSums(Xu) / colSums(X)
    ## Weighted Relative Suitability Index
    WRSI <- Pused / Pavail
    #Var <- (1/colSums(Xu)) - (1/sum(Xu)) + (1/colSums(X)) - (1/sum(X))
    res <- data.frame(
        WRSI=WRSI,
        zWRSI=log(WRSI),
        rWRSI=(exp(2 * log(WRSI)) - 1)/(1 + exp(2 * log(WRSI))),
        Pused=Pused,
        Pavail=Pavail,
        Pw=Pw,
        u=colSums(Xu),
        a=colSums(X))
    rownames(res) <- colnames(X)
    class(res) <- c("wrsi", "data.frame")
    res
}
