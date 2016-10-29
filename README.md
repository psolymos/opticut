# opticut: Likelihood Based Optimal Partitioning for Indicator Species Analysis

Likelihood based optimal partitioning for indicator species analysis and more.
Finding the best binary partition for each species based on
model selection, possibly controlling for modifying/confounding
variables as described in Kemencei et al. (2014).

![](https://github.com/psolymos/opticut/raw/master/extras/oc-logo.png)

User visible changes in the package are listed in the [NEWS](https://github.com/psolymos/opticut/blob/master/NEWS.md) file.

## Versions

[![Build Status](https://travis-ci.org/psolymos/opticut.svg?branch=master)](https://travis-ci.org/psolymos/opticut) [![codecov](https://codecov.io/gh/psolymos/opticut/branch/master/graph/badge.svg)](https://codecov.io/gh/psolymos/opticut)

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
## community data
y <- cbind(
    Sp1=c(4,6,3,5, 5,6,3,4, 4,1,3,2),
    Sp2=c(0,0,0,0, 1,0,0,1, 4,2,3,4),
    Sp3=c(0,0,3,0, 2,3,0,5, 5,6,3,4))

## stratification
g <-    c(1,1,1,1, 2,2,2,2, 3,3,3,3)

## find optimal partitions for each species
oc <- opticut(formula = y ~ 1, strata = g, dist = "poisson")
summary(oc)

## visualize the results
plot(oc, cut = -Inf)

## quantify uncertainty
uc <- uncertainty(oc, type = "asymp", B = 999)
```

## References

Kemencei, Z., Farkas, R., Pall-Gergely, B., Vilisics, F., Nagy, A., Hornung,
E. & Solymos, P. (2014): Microhabitat associations of land snails in
forested dolinas: implications for coarse filter conservation.
_Community Ecology_ **15**:180--186. [[PDF](https://drive.google.com/file/d/0B-q59n6LIwYPWnBjLUxvcXJVUXc/view?usp=sharing)]
