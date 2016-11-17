lorenz <-
function(x)
{
    out <- .lorenz(x)
    out$x <- x
    class(out) <- "lorenz"
    out
}

.lorenz <-
function(x)
{
    o <- order(x)
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
    list(xcut = unname(xo[m1]),
        L = unname(L[m1 + 1]),
        p = unname(p[m1 + 1]),
        S = unname(L[m1 + 1] + p[m1 + 1]),
        G = G,
        J = max(p - L))
}

print.lorenz <-
function(x)
{
    str(x)
}

plot.lorenz <-
function(x)
{
    plot(x)
}

## classifier function: x, g
## - lorenz(x) based cut is used to crosstabulate with g and return prop matrix

## modeling function: combine rankComb and classifier
## summary, plot etc functions for that

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
