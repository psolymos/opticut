#devtools::install_github("psolymos/opticut", ref="multiclass")
library(opticut)
library(mefa4)
library(parallel)
source("~/repos/opticut/extras/ip/ipredict.R")
source("~/repos/opticut/extras/ip/ipredict.multicut.R")
source("~/repos/opticut/extras/ip/ipredict.opticut.R")
source("~/repos/opticut/extras/ip/loo.R")
load("~/Dropbox/collaborations/opticut/R/abmi-data/abmi-data.Rdata")
opticut.formula <- opticut:::opticut.formula
multicut.formula <- opticut:::multicut.formula
.opticut_dist <- opticut:::.opticut_dist

## number of workers
NCL <- 14
## min number of detections
NMIN <- 20
## test is small snail data
TEST <- FALSE

VALS <- expand.grid(
    TAXON=names(ABMI$detections),
    METHOD=c("opticut", "multicut"),
    DIST=c("poisson", "binomial"),
    SCALE=c("sites", "sites_pc"))
VALS <- VALS[!(VALS$TAXON %in% c("vascular_plants", "bryophytes", "lichens") &
    VALS$DIST=="poisson"),]
VALS <- VALS[!(VALS$TAXON != "birds" & VALS$SCALE=="sites_pc"),]
VALS <- VALS[order(VALS$TAXON, VALS$METHOD, VALS$DIST, VALS$SCALE),]

if (TEST) {
    VALS <- droplevels(VALS[VALS$TAXON=="mites",])
}

rn <- rownames(ABMI$sites[ABMI$sites$Year >= 2009 &
    ABMI$sites$NRNAME %in% c("Boreal", "Foothills"),])
for (i in names(ABMI$detections))
    rn <- sort(intersect(rn, rownames(ABMI$detections[[i]])))
#str(rn)

#v <- 1
for (v in 1:nrow(VALS)) {

## "birds" "mites" "vascular_plants" "bryophytes" "lichens"
TAXON <- as.character(VALS$TAXON[v]) # "birds"
## "sites" "sites_pc"
SCALE <- as.character(VALS$SCALE[v]) # "sites"
## binomial, poisson
DIST <- as.character(VALS$DIST[v]) # "binomial"
## opticut, multicut
METHOD <- as.character(VALS$METHOD[v]) # "opticut"

X <- ABMI[[SCALE]][rn,]

table(Succ=cut(X$Succ, c(-1, 0.2, 0.8, 2)),
    Alien=cut(X$Alien, c(-1, 0.2, 0.8, 2)))
X$S <- cut(X$Succ, c(-1, 0.2, 0.8, 2))
levels(X$S) <- c("S0","S1","S2")
X$A <- cut(X$Alien, c(-1, 0.2, 0.8, 2))
levels(X$A) <- c("A0","A1","A2")
X$SA <- interaction(X$S, X$A, drop=TRUE, sep="")
table(X$SA, useNA="a")
g0 <- X$SA

Y <- ABMI$detections[[TAXON]][rn,]
Y <- Y[,colSums(Y>0) >= 5]

if (TAXON=="birds")
    Y <- Y[,colnames(Y) %in%
        rownames(ABMI$species$birds)[ABMI$species$birds$singing]]

if (TEST) {
    data(dolina)
    g0 <- dolina$samp$mhab[dolina$samp$method=="T"]
    Y <- dolina$xtab[dolina$samp$method=="T",]
    Y <- Y[,colSums(Y>0) >= NMIN]
    TAXON <- "dolina"
}
f <- paste("ocip", TAXON, METHOD, DIST,
    ifelse(SCALE=="sites", "ha", "pc"), sep="-")
cat(rep("-", getOption("width")), "\n", f, " ", as.character(Sys.time()),
    "\n", rep("-", getOption("width")), "\n", sep="")
flush.console()


if (DIST == "binomial")
    Y <- ifelse(Y > 0, 1, 0)


cl <- if (NCL > 1)
    makeCluster(NCL) else NULL

if (METHOD=="opticut") {
    o <- opticut(Y ~ 1, strata=g0, dist=DIST, cl=cl)
} else {
    o <- multicut(Y ~ 1, strata=g0, dist=DIST, cl=cl)
}
print(o)
print(o$dist)
print(range(Y))

if (NCL > 1) {
    clusterEvalQ(cl, library(opticut))
    clusterEvalQ(cl, library(pbapply))
    ## these are future/unexported functions for opticut
    clusterExport(cl, c("opticut.formula", "multicut.formula", ".opticut_dist",
        "ipredict", "ipredict.opticut", "ipredict.multicut"))
    ## these are needed because of the update, which is fine
    clusterExport(cl, c("Y", "g0", "DIST"))
}
t0 <- proc.time()
ip <- lotso(o, cl=cl)
t1 <- proc.time()-t0
ip2 <- loto(o, cl=cl)

ip$settings <- list(NMIN=NMIN, TAXON=TAXON, SCALE=SCALE,
    DIST=DIST, METHOD=METHOD, proc_time=t1)

if (NCL > 1)
    stopCluster(cl)

fn <- paste0("~/Dropbox/collaborations/opticut/R/abmi-data/", f, ".Rdata")
save(o, ip, ip2, file=fn)
print(ip$kappa[,1])
}


## peek at the results

fl <- c(
    "ocip-birds-multicut-binomial-ha.Rdata",
    "ocip-birds-multicut-binomial-pc.Rdata",
    "ocip-birds-multicut-poisson-ha.Rdata",
    "ocip-birds-multicut-poisson-pc.Rdata",
    "ocip-birds-opticut-binomial-ha.Rdata",
    "ocip-birds-opticut-binomial-pc.Rdata",
    "ocip-birds-opticut-poisson-ha.Rdata",
    "ocip-birds-opticut-poisson-pc.Rdata",
    "ocip-bryophytes-multicut-binomial-ha.Rdata",
    "ocip-bryophytes-opticut-binomial-ha.Rdata",
    "ocip-lichens-multicut-binomial-ha.Rdata",
    "ocip-lichens-opticut-binomial-ha.Rdata",
    "ocip-mites-multicut-binomial-ha.Rdata",
    "ocip-mites-multicut-poisson-ha.Rdata",
    "ocip-mites-opticut-binomial-ha.Rdata",
    "ocip-mites-opticut-poisson-ha.Rdata",
    "ocip-vascular_plants-multicut-binomial-ha.Rdata",
    "ocip-vascular_plants-opticut-binomial-ha.Rdata")

fl <- c("ocip-dolina-multicut-binomial-ha_clNULL.Rdata",
    "ocip-dolina-multicut-binomial-ha.Rdata",
    "ocip-dolina-multicut-poisson-ha_clNULL.Rdata",
    "ocip-dolina-multicut-poisson-ha.Rdata",
    "ocip-dolina-opticut-binomial-ha_clNULL.Rdata",
    "ocip-dolina-opticut-binomial-ha.Rdata",
    "ocip-dolina-opticut-poisson-ha_clNULL.Rdata",
    "ocip-dolina-opticut-poisson-ha.Rdata")


xx <- list()
oo <- list()
for (f in fl) {
    fn <- paste0("~/Dropbox/collaborations/opticut/R/abmi-data/", f)
    e <- new.env()
    load(fn, envir=e)
    xx[[f]] <- e$ip
    oo[[f]] <- e$o
}

tmp <- t(sapply(xx, function(z) z$kappa))
V2 <- data.frame(VALS, a=tmp[,1])
## why is a0 not the same? -- cohen consides col props, use random?
## opti is better than multi: binary is more robust for extrapolation?
## bin/pois not very different
## generally good accuracy
## is it more related to sample size?



ct <- xx[["ocip-vascular_plants-opticut-binomial-ha.Rdata"]]$multiclass$ctable
addmargins(ct)
