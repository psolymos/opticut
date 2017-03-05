## The Lorenz curve is a continuous piecewise linear function
## representing the distribution of income or wealth.
## p_i=i/n, i=1,...,n
## L_i=sum_{j=1}^{i} x_j / sum_{j=1}^{n} x_j
## p_0=L_0=0
lorenz <-
function(x, n = rep(1, length(x)), na.last=TRUE)
{
    if (any(x < 0))
        stop("x must not be < 0")
    o <- order(x, na.last=na.last)
    xo <- x[o]
#    p <- seq_len(length(xo)) / length(xo)
    p <- cumsum(n[o]) / sum(n)
    L <- cumsum(xo) / sum(xo)
    p <- c(0, p)
    L <- c(0, L)
    J <- p - L

    G <- sum(xo * seq_len(length(xo)))
    G <- 2 * G / (length(xo) * sum(xo))
    G <- G - 1 - (1 / length(xo))

    m1 <- which.max(J)
    out <- cbind(p=p, L=L, x=c(0, xo))
    ## x: habitat suitability cutoff is the back scaled
    ##       L value (original x) from the graph
    ## L: cumulative distribution of the metric of interest
    ## p: cumulative distribution of the available population
    ## S: asymmetry is the sum of x and y coordinates
    ##    at the point of slope 1 (symmetry: S=1)
    ## G: Gini coefficient of 0 mean perfect equality,
    ##    values close to 1 indicate high inequality.
    ## Youden index
    attr(out, "summary") <- c(
        "t" = unname(m1),
        "x(t)" = unname(xo[m1]),
        "p(t)" = unname(p[m1 + 1]),
        "L(t)" = unname(L[m1 + 1]),
        "G" = G,
        "S" = unname(L[m1 + 1] + p[m1 + 1]),
        "J" = max(p - L))
    class(out) <- c("lorenz", "matrix")
    out
}
