# Scaling for the Indicator Potential

Transformation of estimated contrasts to indicator potential.

## Usage

``` r
beta2i(x, scale = 1)
```

## Arguments

- x:

  numeric, real valued coefficients.

- scale:

  numeric, scaling constant.

## Value

Returns a numeric vector (I = abs(tanh(x \* scale))).

## Author

Peter Solymos \<psolymos@gmail.com\>

## See also

[`opticut`](opticut.md) and [`multicut`](multicut.md) use the scaled I
values as indicator potential.

[`ocoptions`](ocoptions.md) for setting value for the default scaling
factor.

## Examples

``` r
x <- seq(-5, 5, 0.1)
Col <- occolors(c("red", "blue"))(10)
plot(x, beta2i(x), type = "n")
s <- seq(1, 0.1, -0.1)
for (i in 1:10) {
    lines(x, beta2i(x, scale = s[i]), col = Col[i])
    text(1.5 - 0.2, beta2i(1.5, scale = s[i]), s[i], col = Col[i])
}
```
