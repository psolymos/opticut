kappa <-
function(predicted, reference)
{
    a <- sum(diag(predicted)) / sum(predicted)
    a0 <- sum(diag(reference)) / sum(reference)
    k <- (a - a0) / (1 - a0)
    c(a=a, a0=a0, k=k)
}
