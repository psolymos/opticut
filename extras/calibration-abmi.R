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

VALS <- expand.grid(
    TAXON=names(ABMI$detections),
    METHOD=c("opticut", "multicut"),
    DIST=c("poisson", "binomial"),
    SCALE=c("sites", "sites_pc"))
VALS <- VALS[!(VALS$TAXON %in% c("vascular_plants", "bryophytes", "lichens") &
    VALS$DIST=="poisson"),]
VALS <- VALS[!(VALS$TAXON != "birds" & VALS$SCALE=="sites_pc"),]
VALS <- VALS[order(VALS$TAXON, VALS$METHOD, VALS$DIST, VALS$SCALE),]

rn <- rownames(ABMI$sites[ABMI$sites$Year >= 2009 &
    ABMI$sites$NRNAME %in% c("Boreal", "Foothills"),])
for (i in names(ABMI$detections))
    rn <- sort(intersect(rn, rownames(ABMI$detections[[i]])))
#str(rn)

## number of workers
NCL <- 2
## min number of detections
NMIN <- 20
## test is small snail data
TEST <- FALSE

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

f <- paste("ocip", TAXON, METHOD, DIST,
    ifelse(SCALE=="sites", "ha", "pc"), sep="-")
cat(rep("-", 50), "\n", f, " ", as.character(Sys.time()),
    "\n", rep("-", 50), "\n", sep="")
flush.console()

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
}

if (DIST == "binomial")
    Y <- ifelse(Y > 0, 1, 0)


cl <- makeCluster(NCL)
if (METHOD=="opticut") {
    o <- opticut(Y ~ 1, strata=g0, dist=DIST, cl=cl)
} else {
    o <- multicut(Y ~ 1, strata=g0, dist=DIST, cl=cl)
}

clusterEvalQ(cl, library(opticut))
clusterEvalQ(cl, library(pbapply))
## these are future/unexported functions for opticut
clusterExport(cl, c("opticut.formula", "multicut.formula", ".opticut_dist",
    "ipredict", "ipredict.opticut", "ipredict.multicut"))
## these are needed because of the update, which is fine
clusterExport(cl, c("Y", "g0", "DIST"))

ip <- lotso(o, cl=cl)
ip$settings <- list(NMIN=NMIN, TAXON=TAXON, SCALE=SCALE,
    DIST=DIST, METHOD=METHOD)
stopCluster(cl)

fn <- paste0("~/Dropbox/collaborations/opticut/R/abmi-data/", f, ".Rdata")
if (!TEST)
    save(o, ip, file=fn)
if (TEST)
   break

}
