auc.roc <- function(x, ...) {
    x$inv_spec <- 1-x$fpr
    dx <- diff(x$inv_spec)
    sum(dx * x$tpr[-1]) / sum(dx)
}
