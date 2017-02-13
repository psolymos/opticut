ctable <- function(x, y) {
    if (length(y) != length(x))
        stop("Hey! x and y lengths must match")
    x <- as.factor(x)
    x <- droplevels(x)
    levs <- sort(levels(x))
    x <- factor(x, levels=levs)
    y <- factor(y, levels=levs)
    if (any(is.na(y)))
        stop("NAs, and levels in y that are not in x are not allowed, sorry")
    as.matrix(table(Reference=x, Prediction=y))
}
btable <- function(table) {
    if (ncol(table) != nrow(table))
        stop("dimension mismatch")
    if (!all(colnames(table) == rownames(table)))
        stop("dimnames must match")
    f <- function(i, x) {
        c(tp=sum(x[i,i]), fp=sum(x[-i,i]),
        fn=sum(x[i,-i]), tn=sum(x[-i,-i]))
    }
    out <- sapply(seq_len(ncol(table)), f, x=table)
    colnames(out) <- colnames(table)
    out
}
untable <- function(table) {
    x <- y <- integer(sum(table))
    S <- c(0, cumsum(table))
    i <- row(table)
    j <- col(table)
    for (k in 2:length(S)) {
        if (S[k]-S[k-1] > 0) {
            x[(S[k-1]+1):S[k]] <- i[k-1]
            y[(S[k-1]+1):S[k]] <- j[k-1]
        }
    }
    list(x = factor(rownames(table)[x], levels=rownames(table)),
        y = factor(colnames(table)[y], levels=colnames(table)))
}

multiclass <- function(x, y=NULL, beta=1) {
    if (!is.null(y))
        x <- ctable(x, y)
    N <- sum(x)
    cm <- btable(x)
    Stat <- c(
        Acc = mean((cm["tp",] + cm["tn",]) / N),
        Err = mean((cm["fp",] + cm["fn",]) / N),
        Prec_m = sum(cm["tp",]) / sum(cm["tp",] + cm["fp",]),
        Prec_M = mean(cm["tp",] / (cm["tp",] + cm["fp",])),
        Spec_m = sum(cm["tn",]) / sum(cm["fp",] + cm["tn",]),
        Spec_M = mean(cm["tn",] / (cm["fp",] + cm["tn",])),
        Rec_m = sum(cm["tp",]) / sum(cm["tp",] + cm["fn",]),
        Rec_M = mean(cm["tp",] / (cm["tp",] + cm["fn",])))
    Stat <- c(Stat,
        ## F-score is harmonic mean
        ## G-score is geometric mean sqrt(prec * recall)
        F_m = unname((beta^2 + 1) * Stat["Prec_m"] * Stat["Rec_m"] /
            (beta^2 * Stat["Prec_m"] + Stat["Rec_m"])),
        F_M = unname((beta^2 + 1) * Stat["Prec_M"] * Stat["Rec_M"] /
            (beta^2 * Stat["Prec_M"] + Stat["Rec_M"])))
    Mat <- rbind(Accuracy=(cm["tp",] + cm["tn",]) / N,
        Error=(cm["fp",] + cm["fn",]) / N,
        Precision=cm["tp",] / (cm["tp",] + cm["fp",]),
        Specificity=cm["tn",] / (cm["fp",] + cm["tn",]),
        Recall=cm["tp",] / (cm["tp",] + cm["fn",])) # recall = sensitivity
    #Mat <- rbind(Mat, jouden=Mat["Specificity",]+Mat["Recall",]-1)
    out <- list(ctable=x, btable=cm, average=Stat, beta=beta,
        single=Mat)
    class(out) <- "multiclass"
    out
}

etable <- function(table, type="cohen", w=NULL) {
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
rtable <- function(n, table, type=c("r", "rc")) {
    type <- match.arg(type)
    if (type == "rc") {
        out <- array(unlist(r2dtable(n, rowSums(table), colSums(table))),
            c(dim(table), n))
    }
    if (type == "r") {
        r <- rowSums(table)
        K <- length(r)
        f <- function() {
            t(sapply(seq_len(K), function(i)
                rmultinom(1, r[i], table[i,])))
        }
        out <- replicate(n, f())
    }
    #out <- rmultinom(n, sum(table), table) # this does not keep margins
    dimnames(out) <- list(rownames(table), colnames(table), NULL)
    out
}

kappa <- function(predicted, reference) {
    a <- sum(diag(predicted)) / sum(predicted)
    a0 <- sum(diag(reference)) / sum(reference)
    k <- (a - a0) / (1 - a0)
    c(a=a, a0=a0, k=k)
}

eval_fun <- function(table, n=0, ptype="cohen", rtype="r", w=NULL) {
    ref <- etable(table, type=ptype, w=w)
    D <- c(dim(table), n+1)
    rnd <- if (n > 0)
        c(ref, rtable(n, ref, type=rtype)) else c(ref)
    dim(rnd) <- D
    k <- sapply(seq_len(n+1), function(i) kappa(table, rnd[,,i]))
    k
}
