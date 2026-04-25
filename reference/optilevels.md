# Optimal Number of Factor Levels

Finds the optimal number of factor levels given the data and a model
using a likelihood-based agglomerative algorithm.

## Usage

``` r
optilevels(y, x, z = NULL, alpha = 0, dist = "gaussian", ...)

# S3 method for class 'optilevels'
bestmodel(object, ...)
```

## Arguments

- y:

  vector of observations.

- x:

  a factor or a matrix of proportions (i.e. the values 0 and 1 should
  have consistent meaning across the columns, often through a unit sum
  constraint). It is the user's responsibility to ensure that values
  supplied for `x` are sensible. `x` is not expected to include an
  intercept.

- z:

  a design matrix with predictor variables besides the one(s) defined
  via the argument `x`. It is the user's responsibility to ensure that
  values supplied for `z` are sensible and it also makes sense to bind
  `x` and `z` together. Variables in `z` should be centered (mean 0)
  (and possibly normalized by SD), because the design matrix from `x` is
  not expected to include an intercept.

- alpha:

  numeric \[0-1\], weighting factor for calculating information criteria
  for model selection (i.e. IC = (1-alpha)\*AIC + alpha\*BIC, also
  referred to as CAIC: consistent AIC).

- dist:

  character, distribution argument passed to underlying functions, see
  listed on the help page of [`opticut`](opticut.md) (except for
  `dist = "zip2"`, `dist = "zinb2"` `dist = "rsf"`, and
  `dist = "rspf"`).

- object:

  fitted object.

- ...:

  other arguments passed to the underlying functions, see
  [`opticut1`](opticut.md).

## Value

An object of class 'optilevels' that is a list with the following
elements:

- `"delta"`:

  delta IC values along the selection path considering best models.

- `"ic"`:

  IC values along the selection path considering best models.

- `"coef"`:

  matrix of coefficients (linear predictor scale) corresponding to
  argument `x` along the selection path considering best models.

- `"zcoef"`:

  matrix of coefficients (linear predictor scale) corresponding to
  argument `z` when not `NULL` along the selection path considering best
  models, or `NULL`.

- `"rank"`:

  matrix ranks based on the coefficients along the selection path
  considering best models. Ranking uses the default
  `ties.method = "average"` in
  [`rank`](https://rdrr.io/r/base/rank.html).

- `"deltalist"`:

  delta IC values along the selection path considering all competing
  models.

- `"iclist"`:

  IC values along the selection path considering all competing models.

- `"coeflist"`:

  matrix of coefficients (linear predictor scale) corresponding to
  argument `x` along the selection path considering all competing
  models.

- `"zcoeflist"`:

  matrix of coefficients (linear predictor scale) corresponding to
  argument `z` when not `NULL` along the selection path considering all
  competing models, or `NULL`.

- `"ranklist"`:

  matrix ranks based on the coefficients along the selection path
  considering all competing models.

- `"levels"`:

  list of (merged) factor levels along the selection path considering
  best models.

- `"Y"`:

  vector of observations (argument `y`).

- `"X"`:

  design matrix component corresponding to argument `x`.

- `"Z"`:

  design matrix component corresponding to argument `z`.

- `"alpha"`:

  weighting argument.

- `"dist"`:

  distribution argument.

- `"factor"`:

  logical, indicating if argument `x` is a factor (`TRUE`) or a matrix
  (`FALSE`).

`bestmodel` returns the best supported model for further manipulation
(e.g. prediction).

## Author

Peter Solymos \<psolymos@gmail.com\>

## See also

[`opticut`](opticut.md) and [`multicut`](multicut.md) for fitting best
binary and multi-level response models.

## Examples

``` r
## --- Factor levels with Gaussian distribution
## simple example from Legendre 2013
## Indicator Species: Computation, in
## Encyclopedia of Biodiversity, Volume 4
## https://dx.doi.org/10.1016/B978-0-12-384719-5.00430-5
gr <- as.factor(paste0("X", rep(1:5, each=5)))
spp <- cbind(Species1=rep(c(4,6,5,3,2), each=5),
    Species2=c(rep(c(8,4,6), each=5), 4,4,2, rep(0,7)),
    Species3=rep(c(18,2,0,0,0), each=5))
rownames(spp) <- gr
## must add some noise to avoid perfect fit
spp[6, "Species1"] <- 7
spp[1, "Species3"] <- 17
spp
#>    Species1 Species2 Species3
#> X1        4        8       17
#> X1        4        8       18
#> X1        4        8       18
#> X1        4        8       18
#> X1        4        8       18
#> X2        7        4        2
#> X2        6        4        2
#> X2        6        4        2
#> X2        6        4        2
#> X2        6        4        2
#> X3        5        6        0
#> X3        5        6        0
#> X3        5        6        0
#> X3        5        6        0
#> X3        5        6        0
#> X4        3        4        0
#> X4        3        4        0
#> X4        3        2        0
#> X4        3        0        0
#> X4        3        0        0
#> X5        2        0        0
#> X5        2        0        0
#> X5        2        0        0
#> X5        2        0        0
#> X5        2        0        0

ol <- optilevels(spp[,"Species3"], gr)
ol[c("delta", "coef", "rank", "levels")]
#> $delta
#> [1]  0 -4 NA NA NA
#> 
#> $coef
#>        X1 X2 X3 X4 X5
#> [1,] 17.8  2  0  0  0
#> [2,] 17.8  2  0  0  0
#> [3,]   NA NA NA NA NA
#> [4,]   NA NA NA NA NA
#> [5,]   NA NA NA NA NA
#> 
#> $rank
#>      X1 X2 X3 X4 X5
#> [1,]  5  4  2  2  2
#> [2,]  3  2  1  1  1
#> [3,] NA NA NA NA NA
#> [4,] NA NA NA NA NA
#> [5,] NA NA NA NA NA
#> 
#> $levels
#> $levels[[1]]
#>         X1         X2         X3         X4         X5 
#>       "X1"       "X2" "X3+X4+X5" "X3+X4+X5" "X3+X4+X5" 
#> 
#> $levels[[2]]
#>         X1         X2         X3         X4         X5 
#>       "X1"       "X2" "X3+X4+X5" "X3+X4+X5" "X3+X4+X5" 
#> 
#> 

## get the final factor level
gr1 <- gr
levels(gr1) <- ol$level[[length(ol$level)]]
table(gr, gr1)
#>     gr1
#> gr   X1 X2 X3+X4+X5
#>   X1  5  0        0
#>   X2  0  5        0
#>   X3  0  0        5
#>   X4  0  0        5
#>   X5  0  0        5

## compare the models
o0 <- lm(spp[,"Species3"] ~ gr - 1)
o1 <- lm(spp[,"Species3"] ~ gr1 - 1)
data.frame(AIC(o0, o1), delta=AIC(o0, o1)$AIC - AIC(o0))
#>    df       AIC delta
#> o0  6 -3.103558     0
#> o1  4 -7.103558    -4
ol$delta # should be identical
#> [1]  0 -4 NA NA NA

## --- Proportions with Poisson distribution
## simulation
set.seed(123)
n <- 500 # number of observations
k <- 5 # number of habitat types
b <- c(-1, -0.2, -0.2, 0.5, 1)
names(b) <- LETTERS[1:k]
x <- replicate(k, exp(rnorm(n)))
x <- x / rowSums(x) # proportions
X <- model.matrix(~.-1, data=data.frame(x))
lam <- exp(drop(crossprod(t(X), b)))
y <- rpois(n, lam)

z <- optilevels(y, x, dist="poisson")

## best model refit
bestmodel(z)
#> 
#> Call:  stats::glm(formula = Y ~ . - 1, family = Family, data = XX)
#> 
#> Coefficients:
#>      X1  `X2+X3`       X4       X5  
#> -0.8941  -0.3098   0.5044   1.0391  
#> 
#> Degrees of Freedom: 500 Total (i.e. Null);  496 Residual
#> Null Deviance:       586.2 
#> Residual Deviance: 548.4     AIC: 1297

## estimates
plogis(z$coef)
#>             X1        X2        X3        X4       X5
#> [1,] 0.2880184 0.3761074 0.4665556 0.6289233 0.737896
#> [2,] 0.2902740 0.4231689 0.4231689 0.6234846 0.738678
#> [3,]        NA        NA        NA        NA       NA
#> [4,]        NA        NA        NA        NA       NA
#> [5,]        NA        NA        NA        NA       NA
plogis(b)
#>         A         B         C         D         E 
#> 0.2689414 0.4501660 0.4501660 0.6224593 0.7310586 
## optimal classification
z$rank
#>      X1 X2 X3 X4 X5
#> [1,]  1  2  3  4  5
#> [2,]  1  2  2  3  4
#> [3,] NA NA NA NA NA
#> [4,] NA NA NA NA NA
#> [5,] NA NA NA NA NA

## get the final matrix
x1 <- mefa4::groupSums(x, 2, z$levels[[length(z$levels)]])
head(x)
#>           [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,] 0.2258750 0.21671095 0.14615367 0.17407228 0.23718812
#> [2,] 0.2256225 0.10514521 0.10039234 0.20888497 0.35995493
#> [3,] 0.4995209 0.29345552 0.10323013 0.04264214 0.06115135
#> [4,] 0.1150706 0.22726090 0.09395935 0.20075919 0.36294998
#> [5,] 0.1998835 0.03883339 0.01372342 0.53850776 0.20905194
#> [6,] 0.3048449 0.04987853 0.15529267 0.46033369 0.02965023
head(x1)
#>             X1      X2+X3         X4         X5
#> [1,] 0.2258750 0.36286462 0.17407228 0.23718812
#> [2,] 0.2256225 0.20553755 0.20888497 0.35995493
#> [3,] 0.4995209 0.39668566 0.04264214 0.06115135
#> [4,] 0.1150706 0.32122025 0.20075919 0.36294998
#> [5,] 0.1998835 0.05255681 0.53850776 0.20905194
#> [6,] 0.3048449 0.20517120 0.46033369 0.02965023

## compare the models
m0 <- glm(y ~ x - 1, family="poisson")
m1 <- glm(y ~ x1 - 1, family="poisson")
data.frame(AIC(m0, m1), delta=AIC(m0, m1)$AIC - AIC(m0))
#>    df      AIC     delta
#> m0  5 1298.452  0.000000
#> m1  4 1297.381 -1.071501
z$delta # should be identical
#> [1]  0.000000 -1.071501        NA        NA        NA

if (FALSE) { # \dontrun{
## dolina example with factor
data(dolina)
dolina$samp$stratum <- as.integer(dolina$samp$stratum)
y <- dolina$xtab[dolina$samp$method == "Q", "ppyg"]
x <- dolina$samp$mhab[dolina$samp$method == "Q"]
z <- scale(model.matrix(~ stratum + lmoist - 1,
    dolina$samp[dolina$samp$method == "Q",]))

## without additional covariates
dol1 <- optilevels(y, x, z=NULL, dist="poisson")
dol1$rank
summary(bestmodel(dol1))

## with additional covariates
dol2 <- optilevels(y, x, z, dist="poisson")
dol2$rank
summary(bestmodel(dol2))

## compare the two models
AIC(bestmodel(dol1), bestmodel(dol2))
} # }
```
