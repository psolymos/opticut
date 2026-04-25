# Finding All Possible Binary Partitions

These functions are used to find all possible binary partitions. Finding
all combinations require a classification vector with K \> 1 strata.

## Usage

``` r
allComb(x, collapse)
kComb(k)
checkComb(x)
```

## Arguments

- x:

  a vector for `allComb` (can be of any type but treated as factor, must
  have at least 2 unique values); and a numeric matrix for `checkComb`.

- collapse:

  character, what to paste between levels. Defaults to
  `getOption("ocoptions")$collapse`.

- k:

  numeric, number of levels (strata) in a given classification (K \> 1).

## Value

`kComb` returns a contrast matrix corresponding to all possible binary
partitions of the factor with K levels. Complements are not counted
twice, i.e. (0,0,1,1) is equivalent to (1,1,0,0). The number of such
possible combinations is M = 2^(K - 1) - 1.

`allComb` takes a classification vector with at least 2 levels and
returns a model matrix with binary partitions.

`checkComb` checks if combinations are unique and non-complementary
(misfits are returned as attributes). Returns a logical value.

## Author

Peter Solymos \<psolymos@gmail.com\>

## See also

[`opticut`](opticut.md) for the user interface.

[`rankComb`](rankComb.md) and [`lorenz`](lorenz.md) for alternative
partitioning algorithms.

## Examples

``` r
kComb(k = 2)
#>      [,1]
#> [1,]    1
#> [2,]    0
kComb(k = 3)
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
kComb(k = 4)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,]    1    0    0    0    1    1    1
#> [2,]    0    1    0    0    1    0    0
#> [3,]    0    0    1    0    0    1    0
#> [4,]    0    0    0    1    0    0    1

## finding all combinations
(f <- rep(LETTERS[1:4], each=2))
#> [1] "A" "A" "B" "B" "C" "C" "D" "D"
(mc <- allComb(f, collapse = "_"))
#>   A B C D A_B A_C A_D
#> A 1 0 0 0   1   1   1
#> A 1 0 0 0   1   1   1
#> B 0 1 0 0   1   0   0
#> B 0 1 0 0   1   0   0
#> C 0 0 1 0   0   1   0
#> C 0 0 1 0   0   1   0
#> D 0 0 0 1   0   0   1
#> D 0 0 0 1   0   0   1
#> attr(,"collapse")
#> [1] "_"
#> attr(,"comb")
#> [1] "all"
## checking for complementary entries
checkComb(mc) # TRUE
#> [1] TRUE
#> attr(,"comp")
#>      i j
#> attr(,"same")
#>      i j
## adding complementary entries to the matrix
mc2 <- cbind(z = 1 - mc[,1], mc[,c(1:ncol(mc), 1)])
colnames(mc2) <- 1:ncol(mc2)
mc2
#>   1 2 3 4 5 6 7 8 9
#> A 0 1 0 0 0 1 1 1 1
#> A 0 1 0 0 0 1 1 1 1
#> B 1 0 1 0 0 1 0 0 0
#> B 1 0 1 0 0 1 0 0 0
#> C 1 0 0 1 0 0 1 0 0
#> C 1 0 0 1 0 0 1 0 0
#> D 1 0 0 0 1 0 0 1 0
#> D 1 0 0 0 1 0 0 1 0
checkComb(mc2) # FALSE
#> [1] FALSE
#> attr(,"comp")
#>      i j
#> [1,] 1 2
#> [2,] 1 9
#> attr(,"same")
#>      i j
#> [1,] 9 2
```
