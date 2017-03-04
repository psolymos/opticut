bestmodel.multicut <-
function (object, which=NULL, ...)
{
        Z <- model.matrix(~Z, data.frame(Z=object$strata))[,-1L]
        .extractOpticut(object, which,
            boot=FALSE,
            internal=TRUE,
            full_model=TRUE,
            best=TRUE,
            Z=Z, ...)
}
