## Simulations: K=2 case

## The power of a statistical test is 
## the probability that it correctly rejects the null hypothesis 
## when the null hypothesis is false.

## Here H0 is: no significant diff among hab A and B (pval >= alpha).
## H1: hab A and B are different (pval < alpha).
## in our setup the H0 is false
## so we examine how often we found that H0 is rejected (p < alpha)
## and also that the habitat identified is the right one

library(mefa4)
library(opticut)
library(indicspecies)

B <- 200

oc_sim_K2 <- 
function(b1=0.5, b2=0, b3=0.05, b4=0, mu0=10, n=100, pos=TRUE)
{
    K <- 2
    n <- 2*(n %/% 2)
    i <- 1:n
    h0 <- rep(LETTERS[1:2], each=n/2)
    j <- sample.int(n/2, round((n/2)*(b4/2)))
    h1 <- h0
    h1[j] <- h0[j+n/2]
    h1[j+n/2] <- h0[j]
    mu <- c(mu0*(1-b1), mu0)
    Mean <- rep(mu, each=n/2)
    ## this can remove group effects
    ## same pattern in each group
#    Conf <- runif(n, -0.5*diff(mu), 0.5*diff(mu))
    ## counteracts the habitat effect
    Conf <- runif(n,0,1) * rep(c(diff(mu), -diff(mu)), each=n/2)
    Noise <- rnorm(n, 0, 0.1*mu0)
    Y <- Mean + b2*Conf + b3*Noise    
    if (pos)
        Y[Y < 0] <- 0
    out <- data.frame(Y=Y, h0=h0, h1=h1, i=i, x=i-mean(i),
        Mean=Mean, Conf=Conf, Noise=Noise)
    attr(out, "settings") <- list(b1=b1, b2=b2, b3=b3, b4=b4, 
                                  K=K, n=n, mu0=mu0, pos=pos)
    out
}

est_fun1 <- function(..., R=999, level=0.95) {
    dat <- oc_sim_K2(...)
    m0 <- opticut(Y ~ 1, dat, strata=h1)$species[[1]]
    m1 <- opticut(Y ~ Conf, dat, strata=h1)$species[[1]]
    ## IndVal
    iv <- multipatt(data.frame(spp1=dat$Y), dat$h1, func = "IndVal.g", 
        duleg=TRUE, control = how(nperm=R))
    ## Phi coef
    rg <- multipatt(data.frame(spp1=dat$Y), dat$h1, func = "r.g", 
        duleg=TRUE, control = how(nperm=R))
    ## F-ratio
    fv <- anova(lm(Y ~ h1, dat))
    
    out <- matrix(NA, 5, 2)
    colnames(out) <- c("stat", "pass")
    rownames(out) <- c("I0", "IX", "IV", "PH", "FR")
    out[1,1] <- m0$I[1]
    out[1,2] <- ifelse(m0$logLR >= 2, 1, 0)
    out[2,1] <- m1$I[1]
    out[2,2] <- ifelse(m1$logLR >= 2, 1, 0)
    out[3,1] <- iv$sign[1,"stat"]
    out[3,2] <- ifelse(iv$sign[1,"p.value"] <= 1-level, 1, 0)
    out[4,1] <- rg$sign[1,"stat"]
    out[4,2] <- ifelse(rg$sign[1,"p.value"] <= 1-level, 1, 0)
    out[5,1] <- fv[1,"F value"]
    out[5,2] <- ifelse(fv[1,"Pr(>F)"] <= 1-level, 1, 0)
    if (rownames(m0) != "B")
        out[1,2] <- -out[1,2]
    if (rownames(m1) != "B")
        out[2,2] <- -out[2,2]
    if (iv$sign[1,"s.B"] != 1)
        out[3,2] <- -out[3,2]
    if (rg$sign[1,"s.B"] != 1)
        out[4,2] <- -out[4,2]
    ## F-ratio cannot tell which is low/high
    t(out)    
}

vals <- expand.grid(
    b2=seq(0, 1, by=0.1),
    b4=seq(0, 1, by=0.1))

vals2 <- expand.grid(
    b2=seq(0, 1, by=0.1),
    b4=seq(0, 1, by=0.1),
    b1=c(0.1, 0.5, 0.9),
    b3=c(0.1, 0.5, 1))

library(parallel)
cl <- makeCluster(6)
clusterExport(cl, c("est_fun", "est_fun1", "oc_sim_K2", "B","vals","vals2"))
clusterEvalQ(cl, library(opticut))
clusterEvalQ(cl, library(indicspecies))

## contrast (I)
res1 <- parLapply(cl, seq(0.1, 0.9, by=0.1), function(z) 
    replicate(B, est_fun1(b1=z)))

## confounding
res2 <- parLapply(cl, seq(0, 1, by=0.1), function(z) 
    replicate(B, est_fun1(b2=z)))

## noise
res3 <- parLapply(cl, seq(0.1, 1, by=0.1), function(z) 
    replicate(B, est_fun1(b3=z)))

## mixing
res4 <- parLapply(cl, seq(0, 1, by=0.1), function(z) 
    replicate(B, est_fun1(b4=z)))

## conf & mixing
res5 <- parLapply(cl, 1:nrow(vals), function(z) 
    replicate(B, est_fun1(b2=vals[z,1], b4=vals[z,2])))

## all
res6 <- parLapply(cl, 1:nrow(vals2), function(z) 
    replicate(B, est_fun1(b2=vals2[z,1], b4=vals2[z,2],
        b1=vals2[z,3], b3=vals2[z,4])))

stopCluster(cl)

#save(vals, vals2, res1, res2, res3, res4, res5, res6, B
#    file="~/Dropbox/collaborations/opticut/R/opticut-simuls.Rdata")

load("~/Dropbox/collaborations/opticut/opticut-simuls.Rdata")

r1 <- sapply(res1, function(z) rowMeans(abs(z[2,,])))
r2 <- sapply(res2, function(z) rowMeans(abs(z[2,,])))
r3 <- sapply(res3, function(z) rowMeans(abs(z[2,,])))
r4 <- sapply(res4, function(z) rowMeans(abs(z[2,,])))
r5 <- sapply(res5, function(z) rowMeans(abs(z[2,,])))
r6 <- sapply(res6, function(z) rowMeans(abs(z[2,,])))

r1v <- sapply(res1, function(z) rowMeans(z[2,,]>0))
r2v <- sapply(res2, function(z) rowMeans(z[2,,]>0))
r3v <- sapply(res3, function(z) rowMeans(z[2,,]>0))
r4v <- sapply(res4, function(z) rowMeans(z[2,,]>0))
r5v <- sapply(res5, function(z) rowMeans(z[2,,]>0))
r6v <- sapply(res6, function(z) rowMeans(z[2,,]>0))
summary(t(r1-r1v))

op <- par(mfrow=c(2,5))
rr <- r6v[,vals2$b1==0.1 & vals2$b3==0.1]
for (i in 1:5)
image(unique(vals$b2), unique(vals$b4), 
    matrix(-rr[i,], length(unique(vals$b2)), length(unique(vals$b4))),
    main=rownames(rr)[i],
    xlab="Confounding", ylab="Misclassification")

rr <- r6v[,vals2$b1==0.9 & vals2$b3==1]
for (i in 1:5)
image(unique(vals$b2), unique(vals$b4), 
    matrix(-rr[i,], length(unique(vals$b2)), length(unique(vals$b4))),
    main=rownames(rr)[i],
    xlab="Confounding", ylab="Misclassification")
par(op)

r6c <- sapply(res6, function(z) rowSums(z[2,,]>0))
xx <- data.frame(success=as.numeric(t(r6c)), 
    failure=B-as.numeric(t(r6c)),
    method=rep(rownames(r6c), nrow(vals2)), vals2)
xx$method <- relevel(xx$method, "IX")
mm <- glm(cbind(success, failure) ~ method + (b1 + b2 + b3 + b4), xx, family=binomial)
summary(mm)
av <- anova(mm)
av$Perc <- round(100 * anova(mm)$Deviance / 1041142, 2)
sum(av$Perc, na.rm=TRUE)
av


## Gaussian example

library(opticut)
Y <- c(0, 0, 3, 0, 2, 3, 0, 5, 5, 6, 3, 4)
z <- as.factor(rep(LETTERS[1:3], each=4))

opticut1(Y, Z=z)
opticut1(ifelse(Y>0,1,0), Z=z, dist="binomial")

print(opticut1(Y, Z=allComb(z)), cut=-Inf)
print(opticut1(ifelse(Y>0,1,0), Z=allComb(z), dist="binomial"), cut=-Inf)
