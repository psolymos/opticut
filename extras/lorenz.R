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
    out <- list(
        data=cbind(p=p, L=L),
        ## xcut: habitat suitability cutoff is the back scaled
        ##       L value (original x) from the graph
        lambda = unname(xo[m1]),
        ## L: cumulative distribution of the metric of interest
        L = unname(L[m1 + 1]),
        ## p: cumulative distribution of the available population
        p = unname(p[m1 + 1]),
        ## S: asymmetry is the sum of x and y coordinates
        ##    at the point of slope 1 (symmetry: S=1)
        S = unname(L[m1 + 1] + p[m1 + 1]),
        ## G: Gini coefficient of 0 mean perfect equality,
        ##    values close to 1 indicate high inequality.
        G = G,
        ## Youden index
        J = max(p - L))
    class(out) <- "lorenz"
    out
}

print.lorenz <-
function(x, digits, ...)
{
    if (missing(digits))
        digits <- max(3L, getOption("digits") - 3L)
    cat("Lorenz curve statistics\n\n",
    "  Threshold = ", format(x$lambda, digits = digits),
    " (L = ", format(x$L, digits = digits),
    ", p = ", format(x$p, digits = digits), ")\n",
    "  Youden index = ", format(x$J, digits = digits), "\n",
    "  Gini index = ", format(x$G, digits = digits), "\n",
    "  Asymmetry = ", format(x$J, digits = digits), "\n\n", sep="")
    invisible(x)
}

## remove diags, add=T/F for adding the line to another one
plot.lorenz <-
function(x, lwd=2, diag=NA, antidiag=NA, tangent=NA, threshold=NA, ...)
{
    plot(x$data, type="l", xaxs = "i", yaxs = "i", lwd=lwd, ...)
    if (!is.na(diag))
        abline(0, 1, col=diag)
    if (!is.na(antidiag))
        abline(1,-1, col=antidiag)
    if (!is.na(threshold)) {
        lines(c(0, x$p), c(x$L, x$L), col=threshold)
        lines(c(x$p, x$p), c(0, x$L), col=threshold)
    }
    if (!is.na(tangent))
        abline(x$L-x$p, 1, col=tangent)
    invisible(x)
}

## classifier function: x, g
## - lorenz(x) based cut is used to crosstabulate with g and return prop matrix

## modeling function: combine rankComb and classifier
## summary, plot etc functions for that using freq of g + use/avail coloring

f_lorenz <- function(j, i, z)
{
    x <- DAT[BB[,j],]
    y <- as.numeric(YY[BB[,j], i])
    off <- if (i %in% colnames(OFF))
        OFF[BB[,j], i] else OFFmean[BB[,j]]
    HABV <- x[,z]
    ## Lorenz-tangent approach for core habitat delineation
    habmod <- glm_skeleton(try(glm(y ~ HABV + ROAD01,
        x,
        family=poisson(),
        offset=off,
        #weights=w,
        x=FALSE, y=FALSE, model=FALSE)))
    ## need to correct for linear effects
    ## so that we estimate potential pop in habitats (and not realized)
    XHSH <- model.matrix(~ HABV + ROAD01, x)
    XHSH[,"ROAD01"] <- 0 # not predicting edge effects
    ## some levels might be dropped (e.g. Marsh)
    XHSH <- XHSH[,names(habmod$coef)]

    ## density based
    lam <- exp(drop(XHSH %*% habmod$coef))
    cv <- Lc_cut(lam, transform=FALSE) # $lam is threshold
    Freq <- table(hab=HABV, lc=ifelse(lam >= cv$lam, 1, 0))
    Prob <- Freq[,"1"] / rowSums(Freq)
    Prob[is.na(Prob)] <- 0
    Hi <- names(Prob)[Prob > 0.5]

    ## probability based
    p <- 1-exp(-lam)
    cv2 <- Lc_cut(lam, transform=FALSE) # $lam is threshold
    Freq2 <- table(hab=HABV, lc=ifelse(p >= cv2$lam, 1, 0))
    Prob2 <- Freq2[,"1"] / rowSums(Freq2)
    Prob2[is.na(Prob2)] <- 0
    Hi2 <- names(Prob2)[Prob2 > 0.5]

    ## final assembly
    out <- list(species=i, iteration=j, habmod=habmod$coef,
        lam=list(
            hi=Hi,
            freq=Freq,
            lc=cv),
        p=list(
            hi=Hi2,
            freq=Freq2,
            lc=cv2))
    out
}

lorenz_breaks <-
function(x, probs=seq(0, 1, 0.1), type=c("lc","pr"),
make_unique=FALSE, digits=6, na.rm=FALSE)
{
    type <- latch.arg(type)
    if (any(probs < 0) || any(probs > 1))
        stop("probs must be in [0,1]")
    probs <- sort(probs)
    if (na.rm)
        x <- x[!is.na(x)]
    o <- order(x)
    x <- cumsum(x[o]) / sum(x)
    if (type=="lc")
        q <- probs
    if (type=="pr")
        q <- quantile(x, probs=probs, na.rm=TRUE)
    xo <- x[o]
    i <- sapply(q, function(z) min(xo[x >= z]))
    names(i) <- format(probs)
    if (make_unique)
        i <- unique(round(i, 6))
    i
}


col2grey <- function(col, method="luminosity") {
    method <- match.arg(method, c("lightness", "average", "luminosity"))
    col <- col2rgb(col) / 255
    if (method == "lightness")
        out <- (apply(col, 2, max) + apply(col, 2, min)) / 2
    if (method == "average")
        out <- colMeans(col)
    if (method == "luminosity")
        out <- 0.21 * col["red",] + 0.72 * col["green",] + 0.07 * col["blue",]
    grey(out)
}

n <- 25
col <- occolors()(n)
plot(0, type="n", ann=FALSE, axes=FALSE,
    xlim=c(0, n), ylim=c(0, 5))
for (i in 1:n) {
    polygon(c(i-1, i, i, i-1), c(0, 0, 1, 1),
        col=col[i], border=col[i])
    polygon(c(i-1, i, i, i-1), c(0, 0, 1, 1)+1,
        col=col2grey(col, "lightness")[i], border=col2grey(col, "lightness")[i])
    polygon(c(i-1, i, i, i-1), c(0, 0, 1, 1)+2,
        col=col2grey(col, "average")[i], border=col2grey(col, "average")[i])
    polygon(c(i-1, i, i, i-1), c(0, 0, 1, 1)+3,
        col=col2grey(col, "luminosity")[i], border=col2grey(col, "luminosity")[i])
    polygon(c(i-1, i, i, i-1), c(0, 0, 1, 1)+4,
        col=col[i], border=col[i])
}
text(rep(n/2, 4), (1:5)-0.5,
    c("color", "lightness", "average", "luminosity", "color"))
