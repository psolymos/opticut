source("~/repos/opticut/R/opticut.R")

set.seed(1234)
n <- 200
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 %in% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
table(x0,x1)
lam1 <- exp(0.5 + 0.5*x1 + -0.2*x2)
boxplot(lam1~x0)
Y1 <- rpois(n, lam1)
lam2 <- exp(0.1 + 0.5*ifelse(x0==4,1,0) + 0.2*x2)
boxplot(lam2~x0)
Y2 <- rpois(n, lam2)
lam3 <- exp(0.1 + -0.2*x2)
boxplot(lam3~x0)
Y3 <- rpois(n, lam3)
Y <- cbind(SPP1=Y1, SPP2=Y2, SPP3=Y3)
X <- model.matrix(~x2)
Z <- modelComb(x0)

opticut1(Y1, X, Z, dist="poisson")
opticut1(Y2, X, Z, dist="poisson")
(mod <- opticut(Y ~ x2, strata=x0, dist="poisson"))
opticut(Y ~ x2, strata=Z[,1:4], dist="poisson")

## dist as function: lme4
library(lme4)
set.seed(1234)
n <- 200
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 %in% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
ee <- rnorm(n/5)
g <- rep(1:5, each=n/5)
lam1 <- exp(0.5 + 0.5*x1 + -0.2*x2 + ee[g])
Y1 <- rpois(n, lam1)

X <- model.matrix(~x2)
Z <- modelComb(x0)

lmefun <- function(Y, X, linkinv, gr, ...) {
    X <- as.matrix(X)
    m <- glmer(Y ~ X-1 + (1|gr), family=poisson("log"), ...)
    list(coef=fixef(m),
        logLik=logLik(m),
        linkinv=family(m)$linkinv)
}
lmefun(Y1, X, gr=g)
.opticut1(Y1, X, Z1=NULL, dist=lmefun, gr=g)
aa <- opticut1(Y1, X, Z, dist=lmefun, gr=g)



## Gaussian
Y <- rnorm(n, log(lam1) + 10, 0.5)
(mod <- opticut(Y ~ x2, strata=x0, dist="gaussian"))

## same with occupancy
set.seed(1234)
n <- 1000
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 %in% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
table(x0,x1)
p1 <- plogis(0.5 + 0.5*x1 + -0.2*x2)
boxplot(p1~x0)
Y1 <- rbinom(n, 1, p1)
p2 <- plogis(0.1 + 0.5*ifelse(x0==4,1,0) + 0.2*x2)
boxplot(p2~x0)
Y2 <- rbinom(n, 1, p2)
Y <- cbind(SPP1=Y1, SPP2=Y2)
X <- model.matrix(~x2)
Z <- modelComb(x0)

summary(opticut(Y ~ x2, strata=x0, dist="binomial"))




### LEGENDRE EXAMPLE

gr <- rep(1:5, each=5)
spp <- cbind(Species1=rep(c(4,6,5,3,2), each=5),
    Species2=c(rep(c(8,4,6), each=5), 4,4,2, rep(0,7)),
    Species3=rep(c(18,2,0,0,0), each=5))
dist="poisson"
Y <- spp[,1]
X <- matrix(1,25,1)
Z <- modelComb(gr)

summary(mod <- opticut(spp ~ 1, strata=gr, dist="poisson"))


library(vegan)
data(mite)
data(mite.env)
mite.env$hab <- with(mite.env, interaction(Shrub, Topo, drop=TRUE))
summary(mod0 <- opticut(as.matrix(mite) ~ 1, strata=mite.env$hab, dist="poisson"))
plot(mod0)

## sequential
system.time(opticut(as.matrix(mite) ~ 1, strata=mite.env$hab, dist="poisson"))
## parallel -- compare system times
library(parallel)
cl <- makeCluster(3)
system.time(opticut(as.matrix(mite) ~ 1, strata=mite.env$hab, dist="poisson", cl=cl))
stopCluster(cl)
## forking -- will not work on Windows
system.time(opticut(as.matrix(mite) ~ 1, strata=mite.env$hab, dist="poisson", cl=3))


summary(mod <- opticut(as.matrix(mite) ~ SubsDens + WatrCont, mite.env,
    strata=mite.env$hab, dist="poisson"))
plot(mod)

s0 <- summary(mod0, cut=-Inf, sort=FALSE)$summary
s <- summary(mod, cut=-Inf, sort=FALSE)$summary
table(s0=s0$split, s=s$split)
table(s0=s0$assoc, s=s$assoc)

## dune data -- cover type info
## http://www.davidzeleny.net/anadat-r/doku.php/en:data:dune
library(vegan)
data(dune)
data(dune.env)

## ordinal regr
## (when nlevels() < 3 use logistic regression instead !!!
Dune <- as.matrix(dune)
#Dune <- Dune[,apply(Dune, 2, function(z) length(unique(z)))>2]
x <- opticut(Dune ~ 1, strata=dune.env$Management, dist="ordered")
summary(x)

X <- matrix(1,20,1)
Z <- modelComb(dune.env$Management)
opticut1(Dune[,1], X, Z=Z, dist="ordered")

library(pbapply)
Dune <- ifelse(as.matrix(dune)>0,1,0)
x <- opticut(Dune ~ 1, strata=dune.env$Management, dist="binomial")
summary(x)

Dune <- as.matrix(dune+0.5) / 10
x <- opticut(Dune ~ 1, strata=dune.env$Management, dist="beta")
summary(x)

## these cutoff values do not work
Dune <- as.matrix(dune)
Dune[dune == 0] <- 0.001
Dune[dune == 1] <- 0.02
Dune[dune == 2] <- 0.1
Dune[dune == 3] <- 2.5
Dune[dune == 4] <- 5
Dune[dune == 5] <- 8.75
Dune[dune == 6] <- 18.75
Dune[dune == 7] <- 37.5
Dune[dune == 8] <- 62.5
Dune[dune == 9] <- 87.5
## need a small constant to get further away from 0
Dune <- (Dune+.05) / 100
#Dune <- Dune[,apply(Dune, 2, function(z) length(unique(z)))>2]
x <- opticut(Dune ~ 1, strata=dune.env$Management, dist="beta")
summary(x)
x <- opticut(Dune[,3] ~ 1, strata=dune.env$Management, dist="beta")
summary(x)

library(parallel)
cl <- makeCluster(3)

X <- 10
Z <- 4
.f <- function(X, Y, Z)
  X*Y*Z
f <- function(X, Y, Z) .f(X, Y, Z)
sapply(1:3, function(k) f(X=X,Y=k,Z=Z))

clusterExport(cl, c("f",".f","X","Z"))
parSapply(cl, 1:3, function(k) f(X=X,Y=k,Z=Z))

rm(X,Z)
ff <- function(X, Y, Z, cl) {
    X <- 10
    Z <- 4
    clusterExport(cl, c("f",".f"))
#    clusterExport(cl, c("X","Z"), envir=parent.frame())
    e <- new.env()
    assign("X", X, envir=e)
    assign("Z", X, envir=e)
    clusterExport(cl, c("X","Z"), envir=e)
    parSapply(cl, 1:3, function(k) f(X=X,Y=k,Z=Z))
}
ff(X, 1:3, Z, cl)


stopCluster(cl)


}


## Lc tangent
if (FALSE) {


g <- rep(1:5, each=5)
spp <- cbind(Species1=rep(c(4,6,5,3,2), each=5),
    Species2=c(rep(c(8,4,6), each=5), 4,4,2, rep(0,7)),
    Species3=rep(c(18,2,0,0,0), each=5))
dist="poisson"
X <- matrix(1,25,1)
#Z <- modelComb(g)
Y <- spp[,2] # ifelse(spp[,2]>0,1,0)
dist="poisson"
#print(opticut1(Y, X, Z, dist=dist), cut=-Inf)

g <- as.factor(g)
Z <- model.matrix(~g)
m <- .opticut1(Y, X, Z1=Z[,-1,drop=FALSE], 
    linkinv=TRUE, dist=dist)
x <- m$linkinv(c(m$coef[1], m$coef[1] + m$coef[2:ncol(Z)]))
Z1 <- lorenzComb(x[as.integer(g)], g)

table(g, Z1)
opticut1(Y, X, Z=modelComb(g), dist=dist)
opticut1(Y, X, Z=Z1, dist=dist)

## Lc works for binomial case, but not perfectly

## same with occupancy
set.seed(1234)
n <- 100
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 %in% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
table(x0,x1)
p1 <- plogis(0.5 + 0.5*x1 + -0.2*x2)
Y1 <- rbinom(n, 1, p1)
p2 <- plogis(0.1 + 0.5*ifelse(x0==4,1,0) + 0.2*x2)
Y2 <- rbinom(n, 1, p2)
X <- model.matrix(~x2)
Z <- modelComb(x0)

Y <- Y1
g <- as.factor(x0)
Z <- model.matrix(~g)
m <- .opticut1(Y, X, Z1=Z[,-1,drop=FALSE], 
    linkinv=TRUE, dist="binomial")
x <- m$linkinv(c(m$coef[1], m$coef[1] + m$coef[2:ncol(Z)]))
Z1 <- lorenzComb(x[as.integer(g)], g)

table(g, Z1)
opticut1(Y, X, Z=modelComb(g), dist="binomial")
opticut1(Y, X, Z=Z1, dist="binomial")

## iterative search based on ranking

Y <- Y2


print(opticut1(Y, X, Z=modelComb(g), dist="binomial"), cut=-Inf)
print(opticut1(Y, X, Z=rankComb(Y, X, g, "binomial"), dist="binomial"), cut=-Inf)

#comb <- sapply(2:20, function(z) ncol(allComb(z)))
comb <- 2^(1:20)-1
M <- 10
plot((2:20)[1:M], comb[1:M], type="l", ylim=c(0,comb[M]), col=2)
lines((2:20)[1:M], (2:20)[1:M]-1, col=4)

cbind(a=comb, b=(2:20)-1)

## presence-only data
## single species model only:
## because the used distr is different for
## each species by definition.

source("c:/Dropbox/pkg/indval/opticut3.R")
library(ResourceSelection)
## settings
n.used <- 1000
m <- 10
n <- n.used * m
set.seed(1234)
x <- data.frame(x0=as.factor(sample(1:3, n, replace=TRUE)), 
    x1=rnorm(n), x2=runif(n))
cfs <- c(1, -0.5, 0.1, -1, 0.5)
## Logistic RSPF model
dd <- simulateUsedAvail(x, cfs, n.used, m, link="logit")

Y <- dd$status
X <- model.matrix(~ x1 + x2, dd)
Z <- modelComb(as.integer(dd$x0))
link="logit"
Z1=NULL

.opticut1(Y, X, Z1=NULL, dist="rspf", m=0, B=0)
.opticut1(Y, X, Z1=Z[,1], dist="rspf", m=0, B=0)
    
(mod <- opticut1(Y, X, Z, dist="rspf", m=0, B=0))




## multi-species opticut

library(rioja)
data(aber)
strat.plot(aber, scale.percent=TRUE, y.rev=TRUE)

Z <- diag(1, nrow(aber), nrow(aber))
Z[lower.tri(Z)] <- 1
Z <- Z[,-1]
Y <- data.matrix(aber)
Y <- Y[,apply(Y, 2, max) > 10]

#mod <- opticut(Y ~ 1, strata=Z, dist="gaussian")

YY <- Y/100
YY[YY==0] <- min(YY[YY>0])
## response cannot be 0
mod <- opticut(YY ~ 1, strata=Z, dist="beta")

xx <- mod$species[[1]]
ll0 <- mean(sapply(mod$species, attr, "logL_null"))
ll <- rowMeans(sapply(mod$species, "[[", "logL"))
w <- exp(ll) / sum(exp(ll))
logLR <- ll - ll0

xx$w <- w
xx$logL <- ll
xx$logLR <- logLR
attr(xx, "logL_null") <- ll0

