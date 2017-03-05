.lorenz_cut1 <-
function(x, g, fix_fitted=FALSE)
{
    g <- as.factor(g)
    sum_by <- function(x, by) {
        X <- t(model.matrix(~g-1))
        rownames(X) <- levels(g)
        cbind(x=as.numeric(X %*% x), by=rowSums(X))
    }
    if (fix_fitted)
        x <- x + abs(min(x))
    ## negative values are problematic for both Lc and rsf
    if (any(x < 0))
        stop("Negative fitted values found.")
    lc <- summary(lorenz(x))
    sb0 <- sum_by(x >= lc["x(t)"], g)
    plc <- sb0[,"x"] / sb0[,"by"]
    bp0 <- ifelse(plc >= 1-lc["p(t)"], 1, 0)
    sb1 <- sum_by(x, g)
    pU <- sb1[,"x"] / sum(sb1[,"x"])
    pA <- sb1[,"by"] / sum(sb1[,"by"])
    bp1 <- ifelse(pU >= pA, 1, 0)
    Mean <- sb1[,"x"] / sb1[,"by"]

    cbind(n=sb0[,"by"], x0=sb0[,"x"], x1=sb1[,"x"], mean=Mean,
         bp_max=ifelse(Mean >= max(Mean), 1, 0),
         p_lc=plc, x_t=lc["x(t)"], p_t=1-lc["p(t)"], bp_lc=bp0,
         p_U=pU, p_A=pA, s=pU/pA, bp_si=bp1)
}
