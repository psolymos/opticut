## The Lorenz curve is a continuous piecewise linear function
## representing the distribution of income or wealth.
## p_i=1/n, i=1,...,n
## L_i=sum_{j=1}^{i} x_j / sum_{j=1}^{n} x_j
## p_0=L_0=0
lorenz <-
function(x, na.last=TRUE)
{
    o <- order(x, na.last=na.last)
    xo <- x[o]
    p <- seq_len(length(xo)) / sum(length(xo))
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

summary.lorenz <-
function (object, ...) {
    out <- attr(object, "summary")
    class(out) <- "summary.lorenz"
    out
}

print.summary.lorenz <-
function(x, digits, ...)
{
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    cat("Lorenz curve summary\n\n")
    print(x[-1], digits = digits)
    cat("\n")
    invisible(x)
}


quantile.lorenz <-
function(x, probs = seq(0, 1, 0.25), ...)
{
    structure(sapply(probs, function(z) min(x[x[,"L"] >= z, "x"])),
        names=paste0(format(100 * probs, trim=TRUE,
        digits = max(2L, getOption("digits"))), "%"))
}

plot.lorenz <-
function(x, type=c("L", "x"), tangent=NA, h=NA, v=NA, ...)
{
    type <- match.arg(type)
    xx <- if (type == "L")
        x[,c("p", "L")] else x[,c("p", "x")]
    ss <- summary(x)
    zz <- if (type == "L")
        list(x=ss["p(t)"], y=ss["L(t)"]) else list(x=ss["p(t)"], y=ss["x(t)"])
    if (type == "x" && !is.na(tangent)) {
        tangent <- NA
        warning("tangent cannot be plotted when type = 'x'")
    }
    plot(xx, type="l", xaxs = "i", yaxs = "i", ...)
    if (!is.na(tangent))
        abline(zz$y-zz$x, 1, col=tangent)
    if (!is.na(h))
        lines(c(0, zz$x), c(zz$y, zz$y), col=h)
    if (!is.na(v))
        lines(c(zz$x, zz$x), c(0, zz$y), col=v)
    invisible(x)
}

x <- c(rexp(100, 10), rexp(200, 1))
l <- lorenz(x)
head(l)
tail(l)
summary(l)
op <- par(mfrow=c(2,1))
plot(l, lwd=2, tangent=2, h=3, v=4)
abline(0, 1, lty=2, col="grey")
abline(1, -1, lty=2, col="grey")
plot(l, type="x", lwd=2, h=3, v=4)
par(op)


## classifier function: x, g
## - lorenz(x) based cut is used to crosstabulate with g and return prop matrix

## modeling function: combine rankComb and classifier
## summary, plot etc functions for that using freq of g + use/avail coloring

set.seed(1234)
n <- 200
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 %in% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
lam <- exp(0.5 + 0.5*x1 + -0.2*x2)
Y <- rpois(n, lam)
X <- model.matrix(~x2)
Z <- as.factor(x0)

lorenzComb <-
function(Y, X, Z, dist="gaussian", collapse, ...)
{
    if (missing(X))
        X <- matrix(1L, length(Y), 1L)
    if (!is.factor(Z))
        stop("Z must be a factor")
    if (length(Y) != length(Z))
        stop("length(Y) must equal length(Z)")
    if (missing(collapse))
        collapse <-  getOption("ocoptions")$collapse
    Z0 <- model.matrix(~Z)
    m <- .opticut1(Y, X, Z1=Z0[,-1L,drop=FALSE],
        dist=dist, full_model=TRUE, ...)
    f <- fitted(m)
    l <- lorenz(f)
    xt <- attr(l, "summary")["x(t)"]
    h <- ifelse(f >= xt, 1, 0)
    tot <- as.matrix(mefa4::Xtab(f ~ Z + h))
    tot <- tot / sum(f)
    freq <- as.matrix(mefa4::Xtab(~ Z + h))
    out <- list(fitted=f, lorenz=l, total=tot, freq=freq, part=h)
    out
}

## http://www.tannerhelland.com/3643/grayscale-image-algorithm-vb6/
col2gray <- function(col, method="BT.709") {
    method <- match.arg(method, c("BT.709", "BT.601",
        "desaturate", "average", "maximum", "minimum",
        "red", "green", "blue"))
    col <- col2rgb(col) / 255
    if (method == "BT.709")
        out <- 0.2126 * col["red",] + 0.7152 * col["green",] + 0.0722 * col["blue",]
    if (method == "BT.601")
        out <- 0.299 * col["red",] + 0.587 * col["green",] + 0.114 * col["blue",]
    if (method == "desaturate")
        out <- (apply(col, 2, max) + apply(col, 2, min)) / 2
    if (method == "average")
        out <- colMeans(col)
    if (method == "maximum")
        out <- apply(col, 2, max)
    if (method == "minimum")
        out <- apply(col, 2, min)
    if (method == "red")
        out <- col["red",]
    if (method == "green")
        out <- col["green",]
    if (method == "blue")
        out <- col["blue",]
    gray(out)
}

col2gray(occolors()(5))
n <- 25
col <- occolors()(n)
plot(0, type="n", ann=FALSE, axes=FALSE,
    xlim=c(0, n), ylim=c(0, 5))
for (i in 1:n) {
    polygon(c(i-1, i, i, i-1), c(0, 0, 1, 1),
        col=col[i], border=col[i])
    polygon(c(i-1, i, i, i-1), c(0, 0, 1, 1)+1,
        col=col2gray(col, "lightness")[i], border=col2gray(col, "lightness")[i])
    polygon(c(i-1, i, i, i-1), c(0, 0, 1, 1)+2,
        col=col2gray(col, "average")[i], border=col2gray(col, "average")[i])
    polygon(c(i-1, i, i, i-1), c(0, 0, 1, 1)+3,
        col=col2gray(col, "luminosity")[i], border=col2gray(col, "luminosity")[i])
    polygon(c(i-1, i, i, i-1), c(0, 0, 1, 1)+4,
        col=col[i], border=col[i])
}
text(rep(n/2, 4), (1:5)-0.5,
    c("color", "lightness", "average", "luminosity", "color"))
