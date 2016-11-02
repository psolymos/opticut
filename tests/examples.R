library(opticut)

## run examples, even \dontrun sections

help_pages <- c("opticut-package",
    "dolina",
    "opticut", "optilevels", "sindex",
    "beta2i",
    "allComb", "rankComb",
    "bestmodel", "uncertainty",
    "occolors", "ocoptions")

for (i in help_pages) {
    cat("\n\n---------- opticut example:", i, "----------\n\n")
    eval(parse(text=paste0("example(", i,
        ", package = 'opticut', run.dontrun = TRUE)")))
}

## testing methods

ocoptions(try_error=TRUE)

set.seed(2345)
n <- 50
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 %in% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
x3 <- ifelse(x0 %in% 2:4, 1, 0)
lam1 <- exp(0.5 + 1*x1 + -0.2*x2)
Y1 <- rpois(n, lam1)
lam2 <- exp(1 + 0.5*x3)
Y2 <- rpois(n, lam2)
Y3 <- rpois(n, exp(0))
Y <- cbind(Spp1=Y1, Spp2=Y2, Spp3=Y3)
oc <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="rank")
oca <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="all")
opticut:::.extractOpticut(oc)
uc <- uncertainty(oc, type="asymp", B=999)

as.data.frame(oc)
as.data.frame(summary(oc))
as.data.frame(uc)
as.data.frame(summary(uc))

fun <- function(Y, X, linkinv, ...) {
    mod <- stats::glm(Y ~ .-1, data=X, family="poisson", ...)
    list(coef=coef(mod),
        logLik=logLik(mod),
        linkinv=family(mod)$linkinv)
}
comb <- allComb(x0)
ocfun <- opticut(Y ~ x2, strata=comb, dist=fun)

strata(oc)
strata(oca)
strata(ocfun)

## testing distributions

#### negbin

summary(opticut(Y ~ x2, strata=x0, dist="negbin", comb="rank"))

#### beta

Y2 <- Y / rowSums(Y)
Y2[Y2 == 0] <- 0.01
Y2[Y2 == 1] <- 0.99
summary(opticut(Y2 ~ x2, strata=x0, dist="beta", comb="rank"))

#### zip2, zinb2

summary(opticut(Y ~ x2, strata=x0, dist="zip2", comb="rank"))
summary(opticut(Y ~ x2, strata=x0, dist="zinb2", comb="rank"))

#### ordered

opticut(Y[,3] ~ x2, strata=x0, dist="ordered", comb="rank")$species

#### rsf, rspf

opticut(ifelse(Y[,3] > 0, 1, 0) ~ x2, strata=x0,
    dist="rspf", comb="rank", m=0)$species
