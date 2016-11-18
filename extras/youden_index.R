## Possoin case, no offset
n <- 200
set.seed(1234)
b <- c(1, -1, 0.5)
x1 <- exp(runif(n))
x2 <- rnorm(n)
lam <- exp(drop(model.matrix(~x1+x2) %*% b))

fun <- function(i) {
    y <- rpois(n, lam)
    y01 <- ifelse(y>0, 1, 0)
    m <- glm(y ~ x1 + x2, family=poisson)
    lam_hat <- exp(drop(model.matrix(~x1+x2) %*% coef(m)))
    cbind(y01,lam_hat)
}


ststfun <- function(x, transform=FALSE) {
    y01 <- x[,1]
    lam_hat <- x[,2]
    rc <- ROC_cut(y01, lam_hat, transform=transform)
    lc <- Lc_cut(lam_hat, transform=transform)
    c(roc=unlist(rc), lc=unlist(lc))
}

library(ineq)
library(pROC)

par(mfrow=c(1,2))
plot(Lc(lam_hat))
abline(h=lc$L, v=lc$p)
abline(lc$L-lc$p, 1)
plot(roc(y01,lam_hat, smooth=TRUE))
abline(h=rc$tpr, v=1-rc$fpr)
abline(rc$tpr+1-rc$fpr, -1)



aa <- data.frame(t(replicate(100, ststfun(fun(), TRUE))))
plot(data.frame(aa))
plot(2*aa$roc.AUC-1, aa$lc.G)
plot(aa$roc.lam, aa$lc.lam)
plot(aa$roc.J, aa$lc.J)



ROC_cut <- function(y01, lam, transform=FALSE) {
    if (transform)
        lam <- 1-exp(-lam)
    df <- data.frame(t(sapply(sort(unique(c(-Inf, lam, Inf))), function(z) {
        CUT <- ifelse(lam >= z, 1, 0)
        c(lam=z, fpr=sum(y01==0 & CUT!=0)/sum(y01==0), tpr=sum(y01==1 & CUT==1)/sum(y01==1))
    })))

    dx <- rev(df$fpr)[-1]-rev(df$fpr)[-length(df$fpr)]
    my <- (rev(df$tpr)[-1]+rev(df$tpr)[-length(df$tpr)])/2
    auc <- sum(dx*my)

    df$qfpr <- qnorm(df$fpr)
    df$qtpr <- qnorm(df$tpr)
    df <- df[is.finite(df$qfpr) & is.finite(df$qtpr),]
    model <- lm(qtpr ~ qfpr, df)
    df$tpr_hat <- pnorm(fitted(model))
    df$J <- df$tpr - df$fpr
    df$J_hat <- df$tpr_hat - df$fpr
    m1 <- which.max(df$J_hat)

    x <- sort(lam)
    G <- sum(x * 1:length(x))
    G <- 2 * G/(length(x) * sum(x))
    G <- G - 1 - (1/length(x))

    list(lam=as.numeric(df$lam[m1]),
        tpr=as.numeric(df$tpr_hat[m1]),
        fpr=as.numeric(df$fpr[m1]),
        G=2*auc-1,
        AUC=auc,
        J=max(df$J_hat))
}

Lc_cut <-
function (lam, transform=FALSE) 
{
    if (transform)
        lam <- 1-exp(-lam)
    o <- order(lam)
    x <- lam[o]
    p <- seq_len(length(x))/sum(length(x))
    L <- cumsum(x)/sum(x)
    p <- c(0, p)
    L <- c(0, L)
    J <- p - L

    G <- sum(x * 1:length(x))
    G <- 2 * G/(length(x) * sum(x))
    G <- G - 1 - (1/length(x))

    m1 <- which.max(J)
    list(lam=unname(x[m1]), L=unname(L[m1+1]), 
        p=unname(p[m1+1]), S=unname(L[m1+1]+p[m1+1]), 
        G=G, J=max(p - L))
}

