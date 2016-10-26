library(opticut)

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
