# Ranking Based Binary Partitions

Blindly fitting a model to all possible partitions is wasteful use of
resources. Instead, one can rank the K levels (strata) based on expected
response values to explore only K-1 binary partitions along the gradient
defined by the ranks of the expected values.

## Usage

``` r
oComb(x, collapse)
rankComb(Y, X, Z, dist = "gaussian", collapse, ...)
```

## Arguments

- Y:

  numeric, vector of observations.

- X:

  numeric, design matrix.

- Z:

  factor, must have at least 2 unique levels.

- dist:

  character, distribution argument passed to underlying functions, see
  listed on the help page of [`opticut`](opticut.md).

- x:

  and a numeric vector.

- collapse:

  character, what to paste between levels. Defaults to
  `getOption("ocoptions")$collapse`.

- ...:

  other arguments passed to the underlying functions, see
  [`opticut`](opticut.md).

## Value

`oComb` returns the 'contrast' matrix based on the rank vector as input.
Ranked from lowest to highest expected value among the partitions.

The function `rankComb` fits the model with multiple (K \> 2) factor
levels to find out the ranking, and returns a binary classification
matrix as returned by `oComb` corresponding to the ranking.

## Author

Peter Solymos \<psolymos@gmail.com\>

## See also

[`allComb`](allComb.md) for alternative partitioning algorithm.

[`opticut`](opticut.md) for the user interface.

## Examples

``` r
## simulate some data
set.seed(1234)
n <- 200
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 %in% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
lam <- exp(0.5 + 0.5*x1 + -0.2*x2)
Y <- rpois(n, lam)

## binary partitions
head(rc <- rankComb(Y, model.matrix(~x2), as.factor(x0), dist="poisson"))
#>   2 1+2 1+2+4
#> 4 0   0     1
#> 4 0   0     1
#> 2 1   1     1
#> 2 1   1     1
#> 1 0   1     1
#> 4 0   0     1
attr(rc, "est") # expected values in factor levels
#>        1        2        3        4 
#> 2.767840 2.792267 1.563285 1.757870 
aggregate(exp(0.5 + 0.5*x1), list(x0=x0), mean) # true values
#>   x0        x
#> 1  1 2.718282
#> 2  2 2.718282
#> 3  3 1.648721
#> 4  4 1.648721

## simple example
oComb(1:4, "+")
#>   1 1+2 1+2+3
#> 1 1   1     1
#> 2 0   1     1
#> 3 0   0     1
#> 4 0   0     0
## using estimates
oComb(attr(rc, "est"))
#>   3 3+4 1+3+4
#> 1 0   0     1
#> 2 0   0     0
#> 3 1   1     1
#> 4 0   1     1
```
