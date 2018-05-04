roc <- function(labels, scores){
    o <- order(scores, decreasing=TRUE)
    labels <- labels[o]
    scores <- scores[o]
    out <- data.frame(
        tpr=cumsum(labels)/sum(labels),
        fpr=cumsum(!labels)/sum(!labels),
        labels=labels,
        scores=scores)
    class(out) <- c("roc", "data.frame")
    out
}
