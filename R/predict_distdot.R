## Z is 0/1 or design matrix without intercept
.predict_dist <- function(cf, dist, link, X, Z, linkinv) {
    K <- NCOL(Z) + 1L
    XX <- cbind(X[,1L,drop=FALSE], Z, X[,-1L,drop=FALSE])
    if (dist %in% c("zip", "zinb")) {
        if (is.na(link))
            link <- "logit"
        p <- binomial(link)$linkinv(cf[length(cf)])
        return((1 - p) * linkinv(drop(XX %*% cf[-length(cf)])))
    }
    if (dist %in% c("zip2", "zinb2")) {
        p <- linkinv(drop(XX[,seq_len(K)] %*% cf[seq_len(K)]))
        lam <- poisson("log")$linkinv(X %*% cf[-seq_len(K)])
        return(drop(p * lam))
    }
    if (dist == "beta") {
        return(linkinv(drop(XX %*% cf[-length(cf)])))
    }
    linkinv(drop(XX %*% cf))
}
