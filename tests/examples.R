#devtools::install_github("psolymos/opticut")
library(opticut)

## --- run examples with \dontrun sections ---

help_pages <- c("opticut-package",
    "dolina", "birdrec",
    "opticut", "optilevels", "multicut",
    "sindex", "lorenz",
    "beta2i",
    "allComb", "rankComb",
    #"bestmodel",
    "uncertainty",
    "occolors", "ocoptions")

for (i in help_pages) {
    cat("\n\n---------- opticut example:", i, "----------\n\n")
    eval(parse(text=paste0("example(", i,
        ", package = 'opticut', run.dontrun = TRUE)")))
}
