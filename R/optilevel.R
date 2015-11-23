## figure out optimal way of:
## - merging factor levels, or
## - summing proportional covariates

## todo:
## - adding in Z matrix for confounders

## estimation

.optilevel <-
function(Y, X, alpha=0, dist="gaussian", ...)
{
    if (!is.function(dist)) {
        dist <- match.arg(dist,
            c("gaussian","poisson","binomial","negbin",
            "beta","zip","zinb","ordered", "rspf"))
        if (dist == "rspf")
            stop("dist='rspf' not supported for optilevel")
    }
    if (is.null(colnames(X)))
        colnames(X) <- paste0("V", seq_len(ncol(X)))
    colnames(X) <- gsub("\\s", "", colnames(X))

    m_full <- .opticut1(Y, X, Z1=NULL, dist=dist, ...)
    cf_full <- m_full$coef
    if (dist %in% c("zip","zinb","beta"))
        cf_full <- cf_full[-length(cf_full)]
    rnk_full <- rank(cf_full)
    k <- length(cf_full)
    n <- length(Y)
    AIC_m_full <- -2*m_full$logLik + 2*k
    BIC_m_full <- -2*m_full$logLik + log(n)*k
    IC_full <- (1-alpha)*AIC_m_full + alpha*BIC_m_full
    cfmat <- matrix(NA, k, k)
    cfmat[1,] <- cf_full
    rnkmat <- matrix(NA, k, k)
    rnkmat[1,] <- rnk_full
    colnames(rnkmat) <- colnames(cfmat) <- colnames(X)
    delta <- rep(NA, k)
    delta[1] <- 0
    ICvec <- rep(NA, k)
    IC_best <- ICvec[1] <- IC_full

    j <- 1
    Delta <- -1
    delta_list <- list()
    IC_list <- list()
    cfmat_list <- list()
    rnkmat_list <- list()
    while (Delta < 0) {
        cfmat_list[[j]] <- matrix(NA, k-j, k)
        rnkmat_list[[j]] <- matrix(NA, k-j, k)
        delta_list[[j]] <- numeric(k-j)
        IC_list[[j]] <- numeric(k-j)
        for (i in seq_len(k-j)) {
            rnk <- rnkmat[j,]
            gr <- as.character(rnk)
            l1 <- unique(gr[rnk == i])
            l2 <- unique(gr[rnk == i+1])
            gr[gr %in% c(l1, l2)] <- paste(l1, l2, sep="+")
            XX <- mefa4::groupSums(X, 2, gr)

            m <- .opticut1(Y, XX, Z1=NULL, dist=dist, ...)
            cf <- m$coef
            if (dist %in% c("zip","zinb","beta"))
                cf <- cf[-length(cf_full)]
            rnk <- rank(cf)
            kk <- length(cf)
            AIC_m <- -2*m$logLik + 2*kk
            BIC_m <- -2*m$logLik + log(n)*kk
            IC <- (1-alpha)*AIC_m + alpha*BIC_m
            IC_list[[j]][i] <- IC
            delta_list[[j]][i] <- IC - IC_best
            names(cf) <- names(rnk) <- colnames(XX)
            cfmat_list[[j]][i,] <- cf[gr]
            rnkmat_list[[j]][i,] <- rnk[gr]
        }
        best <- which.min(delta_list[[j]])
        Delta <- delta_list[[j]][best]
        if (Delta < 0) {
            IC_best <- IC_list[[j]][best]
            ICvec[j+1] <- IC_list[[j]][best]
            delta[j+1] <- Delta
            cfmat[j+1,] <- cfmat_list[[j]][best,]
            rnkmat[j+1,] <- rnkmat_list[[j]][best,]
        }
        j <- j + 1
        if (j >= k)
            break
    }
    list(delta=delta, ic=ICvec, coef=cfmat, rank=rnkmat,
        deltalist=delta_list, iclist=IC_list,
        coeflist=cfmat_list, ranklist=rnkmat_list)
}

optilevel <-
function(y, x, alpha=0, dist="gaussian", ...)
{
    if (is.null(dim(x))) {
        if (!is.factor(x))
            x <- as.factor(x)
        if (nlevels(x) != length(unique(x)))
            stop("zombie (empty) levels in x")
        X <- model.matrix(~x-1)
        colnames(X) <- levels(x)
    } else {
        if (any(colSums(abs(x)) == 0))
            stop("zombie (sum=0) columns in x")
        X <- as.matrix(x)
    }
    out <- .optilevel(Y=y, X=X, alpha=alpha, dist=dist, ...)
    levs <- list()
    for (i in seq_len(length(out$ranklist))) {
        levi <- sapply(1:max(out$rank[i,]), function(j)
            paste(colnames(out$coef)[out$rank[i,] == j], collapse=" "))
        levs[[i]] <- levi[out$rank[i,]]
        names(levs[[i]]) <- colnames(out$coef)
    }
    out$levels <- levs
    out
}
