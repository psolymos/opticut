untable <-
function(table)
{
    x <- y <- integer(sum(table))
    S <- c(0, cumsum(table))
    i <- row(table)
    j <- col(table)
    for (k in 2:length(S)) {
        if (S[k]-S[k-1] > 0) {
            x[(S[k-1]+1):S[k]] <- i[k-1]
            y[(S[k-1]+1):S[k]] <- j[k-1]
        }
    }
    list(x = factor(rownames(table)[x], levels=rownames(table)),
        y = factor(colnames(table)[y], levels=colnames(table)))
}
