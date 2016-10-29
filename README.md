# opticut: Likelihood Based Optimal Partitioning for Indicator Species Analysis

Likelihood based optimal partitioning for indicator species analysis and more.
Finding the best binary partition for each species based on
model selection, possibly controlling for modifying/confounding
variables as described in Kemencei et al. (2014).

![](https://github.com/psolymos/opticut/raw/master/extras/oc-logo.png)

## Versions

[![Build Status](https://travis-ci.org/psolymos/opticut.svg?branch=master)](https://travis-ci.org/psolymos/opticut) [![codecov](https://codecov.io/gh/psolymos/opticut/branch/master/graph/badge.svg)](https://codecov.io/gh/psolymos/opticut)

User visible changes in the package are listed in the [NEWS](https://github.com/psolymos/opticut/blob/master/NEWS.md) file.

Install development version from GitHub:

```
library(devtools)
install_github("psolymos/opticut")
```

## Report a problem

Use the [issue tracker](https://github.com/psolymos/opticut/issues)
to report a problem.

## License

[GPL-2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## Typical workflow

```
library(opticut)

## --- community data ---
y <- cbind(
    Sp1 = c(4,6,3,5, 5,6,3,4, 4,1,3,2),
    Sp2 = c(0,0,0,0, 1,0,0,1, 4,2,3,4),
    Sp3 = c(0,0,3,0, 2,3,0,5, 5,6,3,4))

## --- stratification ---
g <-      c(1,1,1,1, 2,2,2,2, 3,3,3,3)

## --- find optimal partitions for each species ---
oc <- opticut(y, strata = g, dist = "poisson")
summary(oc)

##  Multivariate opticut results, comb = rank, dist = poisson
##
##  Call:
##  opticut.default(Y = y, strata = g, dist = "poisson")
##
##  Best supported models with logLR >= 2:
##      split assoc      I  mu0  mu1 logLR      w
##  Sp3   2 3    ++ 0.6471 0.75 3.50 4.793 0.6962
##  Sp2     3   +++ 0.8571 0.25 3.25 9.203 0.9577
##  2 binary splits
##  1 species not shown

## --- visualize the results ---
plot(oc, cut = -Inf)

## --- quantify uncertainty ---
uc <- uncertainty(oc, type = "asymp", B = 999)
summary(uc)

##  Multivariate opticut uncertainty results
##  type = asymp, B = 999, level = 0.95
##
##      split R      I   Lower  Upper
##  Sp1   1 2 1 0.2860 0.02341 0.5668
##  Sp3   2 3 1 0.6218 0.21456 0.8813
##  Sp2     3 1 0.8274 0.51229 0.9680
```

## References

Kemencei, Z., Farkas, R., Pall-Gergely, B., Vilisics, F., Nagy, A., Hornung,
E. & Solymos, P. (2014): Microhabitat associations of land snails in
forested dolinas: implications for coarse filter conservation.
_Community Ecology_ **15**:180--186. [[PDF](https://drive.google.com/file/d/0B-q59n6LIwYPWnBjLUxvcXJVUXc/view?usp=sharing)]
