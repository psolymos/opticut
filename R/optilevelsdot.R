## figure out optimal way of:
## - merging factor levels, or
## - summing proportional covariates
.optilevels <-
function(Y, X, Z = NULL, alpha=0, dist="gaussian", ...)
{
    if (!is.function(dist)) {
        dist <- match.arg(dist,
            c("gaussian", "poisson", "binomial", "negbin",
            "beta", "zip", "zinb"))
    }
    if (is.null(colnames(X)))
        colnames(X) <- paste0("X", seq_len(ncol(X)))
    colnames(X) <- gsub("\\s", "", colnames(X))
#    nx <- ncol(X)
    X0 <- X
    if (!is.null(Z)) {
        Z <- data.matrix(Z)
        if (is.null(colnames(Z)))
            colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
        colnames(Z) <- gsub("\\s", "", colnames(Z))
        nz <- ncol(Z)
        X <- cbind(X, Z)
    } else {
        nz <- 0
    }

    m_full <- .opticut1(Y, X, Z1=NULL, dist=dist, ...)
    cf_full <- m_full$coef
    ## drop phi and variance components
    if (dist %in% c("zip","zinb","beta"))
        cf_full <- cf_full[-length(cf_full)]

    cf_full_z <- cf_full[(length(cf_full) - nz + 1):length(cf_full)]
    cf_full <- cf_full[1:(length(cf_full) - nz)]

    rnk_full <- rank(cf_full)
    k <- length(cf_full)
    n <- length(Y)
    AIC_m_full <- -2*m_full$logLik + 2*(k+nz)
    BIC_m_full <- -2*m_full$logLik + log(n)*(k+nz)
    IC_full <- (1-alpha)*AIC_m_full + alpha*BIC_m_full
    cfmat <- matrix(NA, k, k)
    cfmat[1,] <- cf_full
    cfzmat <- matrix(NA, k, nz)
    cfzmat[1,] <- cf_full_z
    rnkmat <- matrix(NA, k, k)
    rnkmat[1,] <- rnk_full
    colnames(rnkmat) <- colnames(cfmat) <- colnames(X)[1:(ncol(X) - nz)]
    colnames(cfzmat) <- colnames(Z)
    delta <- rep(NA, k)
    delta[1] <- 0
    ICvec <- rep(NA, k)
    IC_best <- ICvec[1] <- IC_full

    j <- 1
    Delta <- -1
    delta_list <- list()
    IC_list <- list()
    cfmat_list <- list()
    cfzmat_list <- list()
    rnkmat_list <- list()
    while (Delta < 0) {
        cfmat_list[[j]] <- matrix(NA, k-j, k)
        cfzmat_list[[j]] <- matrix(NA, k-j, nz)
        rnkmat_list[[j]] <- matrix(NA, k-j, k)
        delta_list[[j]] <- numeric(k-j)
        IC_list[[j]] <- numeric(k-j)
        for (i in seq_len(k-j)) {
            rnk <- rnkmat[j,]
            gr <- as.character(rnk)
            l1 <- unique(gr[rnk == i])
            l2 <- unique(gr[rnk == i+1])
            gr[gr %in% c(l1, l2)] <- paste(l1, l2, sep="+")
            XX <- mefa4::groupSums(X[,1:(ncol(X) - nz),drop=FALSE], 2, gr)
            if (!is.null(Z))
                XX <- cbind(XX, Z)

            m <- .opticut1(Y, XX, Z1=NULL, dist=dist, ...)
            cf <- m$coef
            if (dist %in% c("zip","zinb","beta"))
                cf <- cf[-length(cf_full)]

            cf_z <- cf[(length(cf) - nz + 1):length(cf)]
            cf <- cf[1:(length(cf) - nz)]

            rnk <- rank(cf)
            kk <- length(cf)
            AIC_m <- -2*m$logLik + 2*(kk+nz)
            BIC_m <- -2*m$logLik + log(n)*(kk+nz)
            IC <- (1-alpha)*AIC_m + alpha*BIC_m
            IC_list[[j]][i] <- IC
            delta_list[[j]][i] <- IC - IC_best
            names(cf) <- names(rnk) <- colnames(XX)[1:(ncol(XX) - nz)]
            cfmat_list[[j]][i,] <- cf[gr]
            cfzmat_list[[j]][i,] <- cf_z
            rnkmat_list[[j]][i,] <- rnk[gr]
        }
        best <- which.min(delta_list[[j]])
        Delta <- delta_list[[j]][best]
        if (Delta < 0) {
            IC_best <- IC_list[[j]][best]
            ICvec[j+1] <- IC_list[[j]][best]
            delta[j+1] <- Delta
            cfmat[j+1,] <- cfmat_list[[j]][best,]
            cfzmat[j+1,] <- cfzmat_list[[j]][best,]
            rnkmat[j+1,] <- rnkmat_list[[j]][best,]
        }
        j <- j + 1
        if (j >= k)
            break
    }
    list(delta=delta, ic=ICvec,
        coef=cfmat, zcoef=if (is.null(Z)) Z else cfzmat,
        rank=rnkmat,
        deltalist=delta_list, iclist=IC_list,
        coeflist=cfmat_list, zcoeflist=if (is.null(Z)) Z else cfzmat_list,
        ranklist=rnkmat_list,
        alpha=alpha, dist=dist, Y=Y, X=X0, Z=Z)
}
