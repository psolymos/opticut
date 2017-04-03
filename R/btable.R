btable <-
function(table)
{
    if (ncol(table) != nrow(table))
        stop("dimension mismatch")
    if (!all(colnames(table) == rownames(table)))
        stop("dimnames must match")
    f <- function(i, x) {
        c(tp=sum(x[i,i]), fp=sum(x[-i,i]),
        fn=sum(x[i,-i]), tn=sum(x[-i,-i]))
    }
    out <- sapply(seq_len(ncol(table)), f, x=table)
    colnames(out) <- colnames(table)
    out
}
