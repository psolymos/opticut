rtable <-
function(n, table, type=c("r", "rc"))
{
    type <- match.arg(type)
    if (type == "rc") {
        out <- array(unlist(r2dtable(n, rowSums(table), colSums(table))),
            c(dim(table), n))
    }
    if (type == "r") {
        r <- rowSums(table)
        K <- length(r)
        f <- function() {
            t(sapply(seq_len(K), function(i)
                rmultinom(1, r[i], table[i,])))
        }
        out <- replicate(n, f())
    }
    #out <- rmultinom(n, sum(table), table) # this does not keep margins
    dimnames(out) <- list(rownames(table), colnames(table), NULL)
    out
}
