# Best model, Partition, and MLE

Generic functions for accessing best model, best partition, and Maximum
Likelihood Estimate from fitted objects.

## Usage

``` r
bestmodel(object, ...)
bestpart(object, ...)
getMLE(object, ...)
```

## Arguments

- object:

  fitted model object.

- ...:

  other arguments passed to the underlying functions.

## Value

`bestmodel` returns the best supported model for further manipulation
(e.g. prediction).

`bestpart` returns a matrix with the best supported partitions for each
species (species as columns).

`getMLE` returns a named list corresponding to the best supported model.
The list has the following elements: `coef` is the Maximum Likelihood
Estimate (MLE), `vcov` is the variance-covariance matrix for the MLE,
`dist` is the distribution inherited from input `object`.

## Author

Peter Solymos \<psolymos@gmail.com\>

## See also

[`opticut`](opticut.md), [`multicut`](multicut.md),
[`uncertainty`](uncertainty.md).
