test_table <-
function(table, FUN=kappa, n=0, type="cohen", method="r", w=NULL, ...)
{
    ref <- etable(table, type=type, w=w)
    D <- c(dim(table), n+1)
    rnd <- if (n > 0)
        c(ref, rtable(n, ref, method=method)) else c(ref)
    dim(rnd) <- D
    k <- sapply(seq_len(n+1),
        function(i, ...) FUN(table, rnd[,,i], ...))
    k
}
