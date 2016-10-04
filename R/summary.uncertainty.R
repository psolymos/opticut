## summary:
## - highest % split label
## - selection freq
## - I
## - confint
## plot also uses level argument
summary.uncertainty <-
function(object, level=0.95, ...)
{
    prob <- c((1-level)/2, 1-(1-level)/2)
    ucl <- lapply(object$uncertainty, function(z)
        rev(sort(table(z$best)))[1L] / (object$B + 1))
#    ucq <- sapply(object$uncertainty, function(z)
#        quantile(z$I, prob))
    ucq <- sapply(1:length(ucl), function(i) {
        z <- object$uncertainty[[i]]
        j <- z$best == names(ucl[[i]])
        c(mean(z$I[j]), quantile(z$I[j], prob))
    })
    P <- unname(unlist(ucl))
    I <- ucq[1,]
    q <- t(ucq[-1,])
    colnames(q) <- c("Lower", "Upper")
    object$uctab <- data.frame(split=sapply(ucl, names),
        P=P, I=I, PI=P*I, P*q)
#        I=object$summary$I,
#        lower=ucq[1L,], upper=ucq[2L,])
    object$level <- level
    class(object) <- "summary.uncertainty"
    object
}
