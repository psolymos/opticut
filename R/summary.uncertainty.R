## summary:
## - highest % split label
## - R = selection freq or reliability
## - I = indicator value or contrast (only for highest)
## - confint for I (only for highest)
summary.uncertainty <-
function(object, level=0.95, ...)
{
    prob <- c((1-level)/2, 1-(1-level)/2)
    ucl <- lapply(object$uncertainty, function(z)
        rev(sort(table(z$best)))[1L] / (object$B + 1))
    ucq <- sapply(1:length(ucl), function(i) {
        z <- object$uncertainty[[i]]
        j <- z$best == names(ucl[[i]])
        c(mean(z$I[j]), quantile(z$I[j], prob))
    })
    R <- unname(unlist(ucl))
    I <- ucq[1,]
    q <- t(ucq[-1,])
    colnames(q) <- c("Lower", "Upper")
    object$uctab <- data.frame(split=sapply(ucl, names),
        R=R, I=I, q)
    object$level <- level
    class(object) <- "summary.uncertainty"
    object
}
