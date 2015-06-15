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

