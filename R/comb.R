## -------------- OptiCut -----------

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

## this takes a classification vector
## with at least 2 levels
## and returns a model matrix with binary partitions
allComb <-
function(x, collapse = getOption("ocoptions")$collapse)
{
    f <- droplevels(as.factor(x))
    LEVELS <- gsub("\\s", "", levels(f))
    i <- as.integer(f)
    n <- max(i)
    s <- seq_len(n)
    ac <- kComb(n)
    LABELS <- apply(ac, 2, function(z)
        paste(LEVELS[as.logical(z)], collapse=collapse))
    out <- apply(ac, 2, function(z) z[match(i, s)])
    rownames(out) <- f
    colnames(out) <- LABELS
    out
}

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


## x is a named vector of ranks, referring to factor levels
## in some classification vector, 1=highest abundance.
oComb <-
function(x, collapse = " ")
{
    if (length(x) < 2L)
        stop("length of x must be >1")
    if (is.null(names(x)))
        names(x) <- seq_len(length(x))
    o <- x[order(x, decreasing = FALSE)]
    out <- diag(1L, length(o))
    out[upper.tri(out)] <- 1L
    out <- out[,-ncol(out)]
    rownames(out) <- names(o)
    colnames(out) <- seq_len(ncol(out))
    for (i in seq_len(ncol(out))) {
        colnames(out)[i] <- paste(names(x)[names(x) %in%
            rownames(out)[out[,i] > 0]],
            collapse = collapse)
    }
    attr(out, "rank") <- o
    out
}

rankComb <-
function(Y, X, Z, dist="gaussian", ...)
{
    if (!is.factor(Z))
        stop("Z must be a factor")
    Z0 <- model.matrix(~Z)
    m <- .opticut1(Y, X, Z1=Z0[,-1L,drop=FALSE],
        linkinv=TRUE, dist=dist, ...)
    lc <- c(m$coef[1], m$coef[1] + m$coef[2:ncol(Z0)])
    names(lc) <- levels(Z)
    x <- rank(-lc)
#    oc <- oComb(x, collapse = getOption("ocoptions")$collapse)
    oc <- oComb(x, collapse = " ")
    out <- oc[match(Z, rownames(oc)),]
    attr(out, "est") <- m$linkinv(lc)
    out
}
