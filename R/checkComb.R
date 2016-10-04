## this checks a design matrix for complementary rows
## e.g. 1100 vs 0011
checkComb <- function(x) {
    n <- NCOL(x)
    if (n < 2)
        return(TRUE)
    mat <- matrix(FALSE, n, n)
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            ## upper.tri
            mat[i,j] <- all(x[,i] == 1-x[,j]) # comp(lementary)
            ## lower.tri
            mat[j,i] <- all(x[,i] == x[,j]) # same
        }
    }
    out <- !any(mat)
    attr(out, "comp") <- cbind(
        i=row(mat)[upper.tri(row(mat))][which(mat[upper.tri(mat)])],
        j=col(mat)[upper.tri(col(mat))][which(mat[upper.tri(mat)])])
    attr(out, "same") <- cbind(
        i=row(mat)[lower.tri(row(mat))][which(mat[lower.tri(mat)])],
        j=col(mat)[lower.tri(col(mat))][which(mat[lower.tri(mat)])])
    out
}
