bestmodel.optilevels <-
function (object, ...)
{
    ## need to combine X and Z for best model
    ## factor vs matrix: all the same
    xg <- mefa4::groupSums(object$X, 2,
        object$levels[[length(object$levels)]])
    object$X <- cbind(xg, object$Z)
    CALL <- object$call
    CALL[[1]] <- as.name(".opticut1")
    CALL$y <- NULL
    CALL$x <- NULL
    CALL$z <- NULL
    CALL$alpha <- NULL
#    CALL$dist <- NULL
    CALL$Y <- as.name("Y")
    CALL$X <- as.name("X")
    CALL$Z1 <- as.name("Z")
    CALL$full_model <- TRUE
    out <- eval(CALL, envir=object)
    out
}
