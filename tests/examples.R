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
comb <- allComb(x0)[,1:6]
attr(comb, "collapse") <- NULL
attr(comb, "comb") <- NULL
colnames(comb) <- letters[1:6]
m1 <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="rank")
m2 <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="all")
m3 <- opticut(Y ~ x2, strata=comb, dist="poisson")
fun <- function(Y, X, linkinv, ...) {
    mod <- stats::glm(Y ~ .-1, data=X, family="poisson", ...)
    list(coef=coef(mod),
        logLik=logLik(mod),
        linkinv=family(mod)$linkinv)
}
m4 <- opticut(Y ~ x2, strata=x0, dist=fun)
ocoptions(try_error=FALSE)

str(m1$strata)
str(m2$strata)
str(m3$strata)
str(m4$strata)

opticut:::.extractOpticut(m1)
opticut:::.extractOpticut(m2)
opticut:::.extractOpticut(m3)
opticut:::.extractOpticut(m4)

u1a <- uncertainty(m1, type="asymp", B=99)
u2a <- uncertainty(m2, type="asymp", B=99)
u3a <- uncertainty(m3, type="asymp", B=99)
## asymp needs Hessian: dist=fun cannot provide that
u4a <- try(uncertainty(m4, type="asymp", B=999), silent=TRUE)
stopifnot(inherits(u4a, "try-error"))

u1b <- uncertainty(m1, type="boot", B=9)
u2b <- uncertainty(m2, type="boot", B=9)
u3b <- uncertainty(m3, type="boot", B=9)
u4b <- uncertainty(m4, type="boot", B=9)

u1c <- uncertainty(m1, type="multi", B=9)
## type=multi cannot use object with comb=all
u2c <- try(uncertainty(m2, type="multi", B=9), silent=TRUE)
stopifnot(inherits(u2c, "try-error"))
## type=multi cannot use object with comb=NA (custom partitions)
u3c <- try(uncertainty(m3, type="multi", B=9), silent=TRUE)
stopifnot(inherits(u3c, "try-error"))
u4c <- uncertainty(m4, type="multi", B=9)

strata(m1)
strata(m2)
strata(m3)
strata(m4)

bestmodel(m1)
bestmodel(m2)
bestmodel(m3)
## dist=fun cannot return the best model (--> uncertainty(type=asymm) fails)
bm4 <- try(bestmodel(m4), silent=TRUE) # dist=fun problem
stopifnot(inherits(bm4, "try-error"))

str(bestpart(m1))
str(bestpart(m2))
str(bestpart(m3))
str(bestpart(m4))
bestpart(u1c)

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(u1a)
summary(u2a)
summary(u3a)

print(m1)
print(m2)
print(m3)
print(m4)
print(u1a)
print(u2a)
print(u3a)

as.data.frame(m1)
as.data.frame(summary(m1))
as.data.frame(u1a)
as.data.frame(summary(u1a))

## testing distributions

#### negbin

summary(m_nb <- opticut(Y ~ x2, strata=x0, dist="negbin", comb="rank"))

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
