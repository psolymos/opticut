multiclass <-
function(x, y=NULL, beta=1)
{
    if (!is.null(y))
        x <- ctable(x, y)
    N <- sum(x)
    cm <- btable(x)
    Stat <- c(
        Acc = mean((cm["tp",] + cm["tn",]) / N),
        Err = mean((cm["fp",] + cm["fn",]) / N),
        Prec_m = sum(cm["tp",]) / sum(cm["tp",] + cm["fp",]),
        Prec_M = mean(cm["tp",] / (cm["tp",] + cm["fp",])),
        Spec_m = sum(cm["tn",]) / sum(cm["fp",] + cm["tn",]),
        Spec_M = mean(cm["tn",] / (cm["fp",] + cm["tn",])),
        Rec_m = sum(cm["tp",]) / sum(cm["tp",] + cm["fn",]),
        Rec_M = mean(cm["tp",] / (cm["tp",] + cm["fn",])))
    Stat <- c(Stat,
        ## F-score is harmonic mean
        ## G-score is geometric mean sqrt(prec * recall)
        F_m = unname((beta^2 + 1) * Stat["Prec_m"] * Stat["Rec_m"] /
            (beta^2 * Stat["Prec_m"] + Stat["Rec_m"])),
        F_M = unname((beta^2 + 1) * Stat["Prec_M"] * Stat["Rec_M"] /
            (beta^2 * Stat["Prec_M"] + Stat["Rec_M"])))
    Mat <- rbind(Accuracy=(cm["tp",] + cm["tn",]) / N,
        Error=(cm["fp",] + cm["fn",]) / N,
        Precision=cm["tp",] / (cm["tp",] + cm["fp",]),
        Specificity=cm["tn",] / (cm["fp",] + cm["tn",]),
        Recall=cm["tp",] / (cm["tp",] + cm["fn",])) # recall = sensitivity
    #Mat <- rbind(Mat, Jouden=Mat["Specificity",]+Mat["Recall",]-1)
    Mat <- rbind(Mat, AUC=0.5*(Mat["Specificity",]+Mat["Recall",]))
    out <- list(ctable=x, btable=cm, average=Stat, beta=beta,
        single=Mat)
    class(out) <- "multiclass"
    out
}
