quantile.lorenz <-
function(x, probs = seq(0, 1, 0.25), ...)
{
    structure(sapply(probs, function(z) min(x[x[,"p"] >= z, "x"])),
        names=paste0(format(100 * probs, trim=TRUE,
        digits = max(2L, getOption("digits"))), "%"))
}
