library(opticut)

## run examples, even \dontrun sections

help_pages <- c("opticut-package",
    "dolina",
    "opticut", "optilevels", "sindex",
    "allComb", "rankComb",
    "bestmodel", "uncertainty",
    "occolors", "ocoptions")

for (i in help_pages) {
    cat("\n\n---------- opticut example:", i, "----------\n\n")
    eval(parse(text=paste0("example(", i,
        ", package = 'opticut', run.dontrun = TRUE)")))
}

## testing coercion methods

set.seed(2345)
n <- 50
x0 <- sample(1:4, n, TRUE)
x1 <- ifelse(x0 %in% 1:2, 1, 0)
x2 <- rnorm(n, 0.5, 1)
x3 <- ifelse(x0 %in% 2:4, 1, 0)
lam1 <- exp(0.5 + 1*x1 + -0.2*x2)
Y1 <- rpois(n, lam1)
lam2 <- exp(1 + 0.5*x3)
Y2 <- rpois(n, lam2)
Y3 <- rpois(n, exp(0))
Y <- cbind(Spp1=Y1, Spp2=Y2, Spp3=Y3)
oc <- opticut(Y ~ x2, strata=x0, dist="poisson", comb="rank")
uc <- uncertainty(oc, type="asymp", B=999)

as.data.frame(oc)
as.data.frame(summary(oc))

as.matrix(uc)
as.matrix(summary(uc))
as.data.frame(uc)
as.data.frame(summary(uc))
