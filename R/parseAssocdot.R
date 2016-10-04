.parseAssoc <- function(x) {
    LRc <- rep(1L, nrow(x))
    LRc[x$logLR > 2] <- 2L
    LRc[x$logLR > 8] <- 3L
    Sign <- c("-","0","+")[x$assoc + 2L]
    Assoc <- character(nrow(x))
    for (i in 1:length(Assoc))
        Assoc[i] <- paste0(rep(Sign[i], LRc[i]), collapse="")
    Assoc[x$assoc == 0] <- "0"
    factor(Assoc, levels=c("---","--","-","0","+","++","+++"))
}
