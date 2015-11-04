source("~/repos/opticut/R/optilevel.R")

## settings
n <- 500 # number of observations
k <- 5 # number of habitat types
#family <- poisson() # type of glm
family <- binomial() # type of glm
prop <- TRUE # factor or proportions
b <- c(-1, -0.2, -0.2, 0.5, 1)
names(b) <- LETTERS[1:k]
#barplot(plogis(b),ylim=c(0,1),ylab="probability", main="True values")
alpha <- 0 # 0=AIC, 1=BIC

## simulation
set.seed(123)
if (prop) {
    x <- replicate(k, exp(rnorm(n)))
    x <- x / rowSums(x) # proportions
    X <- model.matrix(~.-1, data=data.frame(x))
} else {
    x <- as.factor(sample(LETTERS[1:k], n, replace=TRUE))
    X <- model.matrix(~x-1)
}
mu <- drop(X %*% b)
nu <- family$linkinv(mu)
y <- NULL
if (family$family == "poisson")
    y <- rpois(n, nu)
if (family$family == "binomial")
    y <- rbinom(n, 1, nu)
stopifnot(!is.null(y))

## see for yourself

z <- optilevel(y, x, family=family, alpha=0)
## compare last non-NA row to `true` value
family$linkinv(z$coef)
family$linkinv(b)
## optimal classification
z$rank

## with confounding

source("~/repos/opticut/R/optilevel.R")

## settings
n <- 500 # number of observations
k <- 5 # number of habitat types
#family <- poisson() # type of glm
family <- binomial() # type of glm
b <- c(-1, c(-0.2, -0.2, 0.5, 1) - -1)
names(b) <- LETTERS[1:k]
alpha <- 0 # 0=AIC, 1=BIC

## simulation
set.seed(123)
z <- as.factor(sample(LETTERS[1:k], n, replace=TRUE))
Z <- model.matrix(~z)
x <- rnorm(n)
X <- model.matrix(~x)

XX <- cbind(X[,1,drop=FALSE], Z[,-1,drop=FALSE], X[,-1,drop=FALSE])
mu <- drop(XX %*% c(b, -0.5))
nu <- family$linkinv(mu)
y <- NULL
if (family$family == "poisson")
    y <- rpois(n, nu)
if (family$family == "binomial")
    y <- rbinom(n, 1, nu)
stopifnot(!is.null(y))

Y <- y

source("~/repos/opticut/R/opticut.R")
.optilevel <- 
function(Y, X, Z, alpha=0, ...) 
{
    library(mefa4)
    m_full <- .opticut1(Y, X, Z[,-1,drop=FALSE], dist="binomial")
    #m_full <- glm(Y ~ .-1, data=data.frame(X), ...)
    #m_full <- glm(Y ~ .-1, data=data.frame(X), family=family)
    aic <- -2 * m_full$logLik + 2 * length(m_full$coef)
    bic <- -2 * m_full$logLik + log(length(Y)) * length(m_full$coef)
    IC_full <- (1-alpha)*aic + alpha*bic
    cf_full <- c(m_full$coef[1L], m_full$coef[1L] + m_full$coef[-1L])
    rnk_full <- rank(cf_full)
    k <- length(cf_full)

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
    Delta <- 1
    delta_list <- list()
    IC_list <- list()
    cfmat_list <- list() 
    rnkmat_list <- list()
    while (Delta > 0) {
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
            XX <- groupSums(X, 2, gr)
            m <- glm(Y ~ .-1, data=data.frame(XX), ...)
            #m <- glm(Y ~ .-1, data=data.frame(XX), family=family)
            IC <- (1-alpha)*AIC(m) + alpha*BIC(m)
            IC_list[[j]][i] <- IC
            delta_list[[j]][i] <- IC_best - IC
            cf <- coef(m)
            rnk <- rank(cf)
            names(cf) <- names(rnk) <- colnames(XX)
            cfmat_list[[j]][i,] <- cf[gr]
            rnkmat_list[[j]][i,] <- rnk[gr]
        }
        best <- which.max(delta_list[[j]])
        Delta <- delta_list[[j]][best]
        if (Delta > 0) {
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
        ranklist=rnkmat_list, deltalist=delta_list, iclist=IC_list)
}

optilevel <- 
function(y, x, alpha=0, ...)
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
    .optilevel(Y=y, X=X, alpha=alpha, ...)
}


