etable <-
function(table, type="cohen", w=NULL)
{
    type <- match.arg(type,
        c("majority", "random", "weighted", "cohen"))
    N <- sum(table)
    K <- ncol(table)
    p <- rowSums(table) / N
    q <- colSums(table) / N
    if (type == "majority") {
        ## mcc: majority class classifier, No Information Rate (NIR)
        ecm <- table
        ecm[] <- 0
        ecm[,which.max(p)] <- rowSums(table)
    }
    if (type == "random") {
        ## rgc: random-guess classifier
        ecm <- table
        ecm[] <- N/K * p
    }
    if (type == "weighted") {
        ## rwgc: random-weighted-guess classifier
        ecm <- table
        if (is.null(w)) {
            w <- p
        } else {
            w <- w / sum(w)
        }
        ecm[] <- N * p %*% t(p)
    }
    if (type == "cohen") {
        ## Cohen
        ecm <- table
        ecm[] <- N * p %*% t(q)
    }
    ecm
}
