# Options for the opticut Package

A convenient way of handling options related to the opticut package.

## Usage

``` r
ocoptions(...)
```

## Arguments

- ...:

  arguments in `tag = value` form, or a list of tagged values. The tags
  must come from the parameters described below.

## Value

When parameters are set by `ocoptions`, their former values are returned
in an invisible named list. Such a list can be passed as an argument to
`ocoptions` to restore the parameter values. Tags are the following:

- collapse:

  character value to be used when merging factor levels, the default is
  `"+"`.

- cut:

  log likelihood ratio value, model/species with lower values are
  excluded from summaries and plots, the default is `2`.

- sort:

  logical value indicating if species/partitions should be meaningfully
  sorted, the default is `TRUE`. It can take numeric value when only
  species (`1`) or partitions (`2`) are to be sorted (`1:2` is
  equivalent to `TRUE`).

- theme:

  the color theme to be used based on [`occolors`](occolors.md), the
  default is `"br"`.

- check_comb:

  check the design matrices for complementary partitions using
  [`checkComb`](allComb.md), the default is `TRUE`.

- try_error:

  if [`opticut`](opticut.md) and [`multicut`](multicut.md) should
  [`try`](https://rdrr.io/r/base/try.html) to exclude species where the
  models failed (`TRUE`), the default is to stop when an error is
  encountered (`FALSE`).

- scale:

  the scaling factor used to calculate indicator potential (I) based on
  the estimated contrast (x): I = abs(tanh(x\*scale)), the default is
  0.5.

- fix_fitted:

  [`bestpart.multicut`](multicut.md) uses [`lorenz`](lorenz.md) which
  requires nonnegative fitted values, however models with identity link
  can lead to negative expected values. When `TRUE` the fitted
  values (x) are adjusted as x' = x + abs(min(x)) to ensure
  nonnegativity. The default is `FALSE`.

- robust_loglik:

  if ill-defined models resulting in perfect fit (infinite log
  likelihood, or `NA`, `NaN`) should be allowed. The default `TRUE`
  makes such ill-defined log likelihoods a very small real number
  `-(.Machine$double.xmax^(1/3))`. `FALSE` is equivalent to allowing
  every model to safeguard against such cases or not.

## Author

Peter Solymos \<psolymos@gmail.com\>

## Examples

``` r
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

## current settings
print(unlist(ocoptions())) # these give identical answers
#>      collapse           cut          sort         theme    check_comb 
#>           "+"           "2"        "TRUE"          "br"        "TRUE" 
#>     try_error         scale    fix_fitted robust_loglik 
#>       "FALSE"         "0.5"       "FALSE"        "TRUE" 
unlist(getOption("ocoptions"))
#>      collapse           cut          sort         theme    check_comb 
#>           "+"           "2"        "TRUE"          "br"        "TRUE" 
#>     try_error         scale    fix_fitted robust_loglik 
#>       "FALSE"         "0.5"       "FALSE"        "TRUE" 
summary(ocall <- opticut(spp ~ 1, strata=gr, dist="gaussian", comb="all"))
#> Multivariate opticut results, comb = all, dist = gaussian
#> 
#> Call:
#> opticut.formula(formula = spp ~ 1, strata = gr, dist = "gaussian", 
#>     comb = "all")
#> 
#> Best supported models with logLR >= 2:
#>          split assoc      I mu0  mu1 logLR      w
#> Species2 X1+X3   +++ 0.9866 2.0  7.0 14.82 0.4995
#> Species1 X2+X3   +++ 0.8617 3.0  5.6 16.74 0.7035
#> Species3    X1   +++ 1.0000 0.5 17.8 54.26 1.0000
#> 15 binary splits
#> 

## resetting pboptions and checking new settings
ocop <- ocoptions(collapse="&", sort=FALSE)
unlist(getOption("ocoptions"))
#>      collapse           cut          sort         theme    check_comb 
#>           "&"           "2"       "FALSE"          "br"        "TRUE" 
#>     try_error         scale    fix_fitted robust_loglik 
#>       "FALSE"         "0.5"       "FALSE"        "TRUE" 
## running again with new settings
summary(ocall <- opticut(spp ~ 1, strata=gr, dist="gaussian", comb="all"))
#> Multivariate opticut results, comb = all, dist = gaussian
#> 
#> Call:
#> opticut.formula(formula = spp ~ 1, strata = gr, dist = "gaussian", 
#>     comb = "all")
#> 
#> Best supported models with logLR >= 2:
#>          split assoc      I mu0  mu1 logLR      w
#> Species1 X2&X3   +++ 0.8617 3.0  5.6 16.74 0.7035
#> Species2 X1&X3   +++ 0.9866 2.0  7.0 14.82 0.4995
#> Species3    X1   +++ 1.0000 0.5 17.8 54.26 1.0000
#> 15 binary splits
#> 

## resetting original
ocoptions(ocop)
unlist(getOption("ocoptions"))
#>      collapse           cut          sort         theme    check_comb 
#>           "+"           "2"        "TRUE"          "br"        "TRUE" 
#>     try_error         scale    fix_fitted robust_loglik 
#>       "FALSE"         "0.5"       "FALSE"        "TRUE" 
```
