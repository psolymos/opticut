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

subset(m1, c(3,1))
subset(m2, c(TRUE, FALSE, TRUE))
subset(m3, c("Spp1", "Spp3"))

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

summary(subset(u1b, c(3,1)))
summary(subset(u2b, c(TRUE, FALSE, TRUE)))
summary(subset(u3b, c("Spp1", "Spp3")))

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
strata(u1b)
strata(u2b)
strata(u3b)
strata(u4b)

bestmodel(m1)
bestmodel(m2)
bestmodel(m2)
bestmodel(m3)
## dist=fun cannot return the best model (--> uncertainty(type=asymm) fails)
bm4 <- try(bestmodel(m4), silent=TRUE) # dist=fun problem
stopifnot(inherits(bm4, "try-error"))

str(bestpart(m1))
str(bestpart(m2))
str(bestpart(m2, pos_only=TRUE))
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

fitted(m1)
fitted(m2)
fitted(m3)
f4 <- try(fitted(m4), silent=TRUE)
stopifnot(inherits(f4, "try-error"))

predict(m1)
predict(m2)
predict(m3)
p4 <- try(predict(m4), silent=TRUE)
stopifnot(inherits(p4, "try-error"))
predict(m1, gnew=x0, xnew=data.frame(x2=x2))
pr2 <- try(predict(m2, gnew=x0, xnew=data.frame(x2=x2)), silent=TRUE)
stopifnot(inherits(pr2, "try-error"))
pr3 <- try(predict(m3, gnew=x0, xnew=data.frame(x2=x2)), silent=TRUE)
stopifnot(inherits(pr3, "try-error"))
pr4 <- try(predict(m4, gnew=x0, xnew=data.frame(x2=x2)), silent=TRUE)
stopifnot(inherits(pr4, "try-error"))

as.data.frame(m1)
as.data.frame(summary(m1))
as.data.frame(u1a)
as.data.frame(summary(u1a))

## --- testing distributions ---

ocoptions(cut=-Inf)

## gaussian

summary(o <- opticut(Y ~ x2, strata=x0, dist="gaussian"))
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=2)
u <- uncertainty(o, type="multi", B=2)

## poisson

summary(o <- opticut(Y ~ x2, strata=x0, dist="poisson"))
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=2)
u <- uncertainty(o, type="multi", B=2)

## binomial

Y1 <- ifelse(Y > 0, 1, 0)
summary(o <- opticut(Y1 ~ x2, strata=x0, dist="binomial"))
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=2)
u <- uncertainty(o, type="multi", B=2)

## negbin

summary(o <- opticut(Y ~ x2, strata=x0, dist="negbin"))
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=2)
u <- uncertainty(o, type="multi", B=2)

## beta

Y2 <- Y / rowSums(Y)
Y2[Y2 == 0] <- 0.01
Y2[Y2 == 1] <- 0.99
summary(o <- opticut(Y2 ~ x2, strata=x0, dist="beta"))
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=2)
u <- uncertainty(o, type="multi", B=2)

## zip

Yzi <- Y
Yzi[1,] <- 0
#B <- replicate(2, sample.int(n, replace=TRUE))
B <- sapply(2:3, function(i) which((1:n) != i)) # jackknife
B[1,] <- 1
summary(o <- opticut(Yzi ~ x2, strata=x0, dist="zip"))
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=B)
u <- uncertainty(o, type="multi", B=B)

## zinb

summary(o <- opticut(Yzi ~ x2, strata=x0, dist="zinb"))
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=B)
u <- uncertainty(o, type="multi", B=B)

## zip2

summary(o <- opticut(Yzi ~ x2, strata=x0, dist="zip2"))
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=B)
u <- uncertainty(o, type="multi", B=B)

## zinb2

summary(o <- opticut(Yzi ~ x2, strata=x0, dist="zinb2"))
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=B)
u <- uncertainty(o, type="multi", B=B)

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
o <- opticut(Y ~ x1 + x2, dd, strata=x0, dist="rsf")
o$species
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=2)
u <- uncertainty(o, type="multi", B=2)
## intercept + partition
o <- opticut(Y ~ 1, dd, strata=x0, dist="rsf")
o$species
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=2)
u <- uncertainty(o, type="multi", B=2)

## rspf

o <- opticut(Y ~ x1 + x2, dd, strata=x0, dist="rspf")
o$species
u <- uncertainty(o, type="asymp", B=9)
u <- uncertainty(o, type="boot", B=2)
u <- uncertainty(o, type="multi", B=2)

## --- ... in uncertainty should produce an error ---

set.seed(1234)
n <- 50
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 %in% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
lam <- exp(0.5 + 1*x1 + -0.2*x2)
A <- ifelse(x0 %in% c(1,3), 1, 10)
Y <- rpois(n, lam*A)

## no offset: incorrect
no <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="rank")
## with offsets: log Area
wo <- opticut(Y ~ x2, strata=x0, dist="poisson",
    offset=log(A), weights=rep(1,n), comb="rank")
no$species[[1]]
wo$species[[1]]

nu <- uncertainty(no, type="multi", B=2)
## passing ... is not enough for resampling, treated as user error
wu <- try(uncertainty(wo, type="multi", B=2), silent=TRUE)
stopifnot(inherits(wu, "try-error"))

## --- zip2 & zinb2 coef inversion

## implementation:
## - MLE returns unmodified coefs (P of 0 in ZI)
## - .opticut1 returns:
##       -1*coef[1:2]
##       linkinv: binomial(link)$linkinv(eta)
## - asymp uncertainty uses MLE, thus have to invert and use linkinv after

## less 0 in g=1 stratum: assoc is 1+ or 2-
yzi <- c(rep(0, 10), rpois(40, 6), rep(0, 30), rpois(20, 4))
g <- rep(1:2, each=50)
table(yzi, g)
o1 <- opticut(yzi, strata=g, dist="zip2")
o2 <- opticut(yzi, strata=g, dist="zinb2")
## assoc must be positive for comb=rank
stopifnot(o1$species[[1]]$assoc == 1)
stopifnot(o2$species[[1]]$assoc == 1)
## MLE is negative (prob of 0)
stopifnot(getMLE(o1, 1)$coef[2] < 0)
stopifnot(getMLE(o2, 1)$coef[2] < 0)
stopifnot(coef(bestmodel(o1, 1)[[1]], "zero")[2] < 0)
stopifnot(coef(bestmodel(o2, 1)[[1]], "zero")[2] < 0)
## uncertainty should show 1+ (mu0 < mu1)
o1$species[[1]]
o2$species[[1]]
u <- uncertainty(o1, type="asymp", B=9)$uncertainty[[1]]
stopifnot(all(u$mu0 < u$mu1))
u <- uncertainty(o2, type="asymp", B=9)$uncertainty[[1]]
stopifnot(all(u$mu0 < u$mu1))
B <- sapply(2:10, function(i) which((1:length(g)) != i)) # jackknife
u <- uncertainty(o1, type="boot", B=B)$uncertainty[[1]]
stopifnot(all(u$mu0 < u$mu1))
u <- uncertainty(o2, type="boot", B=B)$uncertainty[[1]]
stopifnot(all(u$mu0 < u$mu1))
u <- uncertainty(o1, type="multi", B=B)$uncertainty[[1]]
stopifnot(all(u$mu0 < u$mu1))
u <- uncertainty(o2, type="multi", B=B)$uncertainty[[1]]
stopifnot(all(u$mu0 < u$mu1))
