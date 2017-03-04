#devtools::install_github("psolymos/opticut")
library(opticut)

## --- testing methods ---

ocoptions(try_error=TRUE)
set.seed(2345)
n <- 50
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 %in% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
x3 <- ifelse(x0 %in% 2:4, 1, 0)
mu1 <- 0.5 + 1*x1 + -0.2*x2
mu2 <- 1 + 0.5*x3
mu3 <- rep(0, n)
Y1 <- rpois(n, exp(mu1))
Y2 <- rpois(n, exp(mu2))
Y3 <- rpois(n, exp(mu3))
Y <- cbind(Spp1=Y1, Spp2=Y2, Spp3=Y3)
m1 <- multicut(Y ~ x2, strata=x0, dist="poisson")
fun <- function(Y, X, linkinv, ...) {
    mod <- stats::glm(Y ~ .-1, data=X, family="poisson", ...)
    list(coef=coef(mod),
        logLik=logLik(mod),
        linkinv=family(mod)$linkinv)
}
m4 <- opticut(Y ~ x2, strata=x0, dist=fun)
ocoptions(try_error=FALSE)

subset(m1, c(3,1))
subset(m1, c(TRUE, FALSE, TRUE))
subset(m1, c("Spp1", "Spp3"))

str(m1$strata)

strata(m1)
strata(m4)

bestmodel(m1)
## dist=fun cannot return the best model (--> uncertainty(type=asymm) fails)
bm4 <- try(bestmodel(m4), silent=TRUE) # dist=fun problem
stopifnot(inherits(bm4, "try-error"))

getMLE(m1, 1, vcov=FALSE)
mle4 <- try(getMLE(m4, 1, vcov=FALSE), silent=TRUE)
stopifnot(inherits(mle4, "try-error"))
getMLE(m1, 1, vcov=TRUE)
mle4 <- try(getMLE(m4, 1, vcov=TRUE), silent=TRUE)
stopifnot(inherits(mle4, "try-error"))

summary(m1)
summary(m4)

print(m1)
print(m4)

summary(fitted(m1))
f4 <- try(fitted(m4), silent=TRUE)
stopifnot(inherits(f4, "try-error"))

summary(predict(m1))
p4 <- try(predict(m4), silent=TRUE)
stopifnot(inherits(p4, "try-error"))
summary(predict(m1, gnew=x0, xnew=data.frame(x2=x2)))
pr4 <- try(predict(m4, gnew=x0, xnew=data.frame(x2=x2)), silent=TRUE)
stopifnot(inherits(pr4, "try-error"))

as.data.frame(m1)
as.data.frame(summary(m1))

## --- testing distributions ---

ocoptions(cut=-Inf)

## gaussian

summary(o <- multicut(Y ~ x2, strata=x0, dist="gaussian"))
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=x0, xnew=data.frame(x2=x2)))
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## poisson

summary(o <- multicut(Y ~ x2, strata=x0, dist="poisson"))
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=x0, xnew=data.frame(x2=x2)))
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## binomial

Y1 <- ifelse(Y > 0, 1, 0)
summary(o <- multicut(Y1 ~ x2, strata=x0, dist="binomial"))
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=x0, xnew=data.frame(x2=x2)))
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## negbin

summary(o <- multicut(Y ~ x2, strata=x0, dist="negbin"))
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=x0, xnew=data.frame(x2=x2)))
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## beta

Y2 <- Y / rowSums(Y)
Y2[Y2 == 0] <- 0.01
Y2[Y2 == 1] <- 0.99
summary(o <- multicut(Y2 ~ x2, strata=x0, dist="beta"))
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=x0, xnew=data.frame(x2=x2)))
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## zip

Yzi <- Y
Yzi[1,] <- 0
summary(o <- multicut(Yzi ~ x2, strata=x0, dist="zip"))
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=x0, xnew=data.frame(x2=x2)))
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## zinb

summary(o <- multicut(Yzi ~ x2, strata=x0, dist="zinb"))
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=x0, xnew=data.frame(x2=x2)))
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## zip2

summary(o <- multicut(Yzi ~ x2, strata=x0, dist="zip2"))
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=x0, xnew=data.frame(x2=x2))) #--- FIXME!!!
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## zinb2

summary(o <- multicut(Yzi ~ x2, strata=x0, dist="zinb2"))
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=x0, xnew=data.frame(x2=x2))) # -- FIXME!!!
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## rsf

library(ResourceSelection)
n.used <- 1000
m <- 1
n <- n.used * m
set.seed(1234)
x <- data.frame(x0=as.factor(sample(1:3, n, replace=TRUE)),
    x1=rnorm(n), x2=runif(n))
cfs <- c(1, -0.5, 0.1, -1, 0.5)
dd <- simulateUsedAvail(x, cfs, n.used, m, link="logit")

Y <- dd$status
X <- model.matrix(~ x1 + x2, dd)

## intercept + partition + covariates
o <- multicut(Y ~ x1 + x2, dd, strata=dd$x0, dist="rsf")
o$species
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=o$strata, xnew=o$X))
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## intercept + partition
o <- multicut(Y ~ 1, dd, strata=x0, dist="rsf")
o$species
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=o$strata, xnew=NULL))
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## rspf

o <- multicut(Y ~ x1 + x2, dd, strata=x0, dist="rspf")
o$species
summary(fitted(o))
summary(predict(o))
summary(predict(o, gnew=o$strata, xnew=o$X))
getMLE(o, 1, vcov=FALSE)
getMLE(o, 1, vcov=TRUE)

## --- testing ... passing ---

set.seed(1234)
n <- 50
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 %in% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
lam <- exp(0.5 + 1*x1 + -0.2*x2)
A <- ifelse(x0 %in% c(1,3), 1, 10)
Y <- rpois(n, lam*A)

## no offset: incorrect
no <- multicut(Y ~ x2, strata=x0, dist="poisson")
## with offsets: log Area
wo <- multicut(Y ~ x2, strata=x0, dist="poisson",
    offset=log(A), weights=rep(1,n))
no$species[[1]]
wo$species[[1]]

