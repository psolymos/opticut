bestmodel.multicut <-
function (object, which=NULL, ...)
{
        .extractOpticut(object, which,
            boot=FALSE,
            internal=TRUE,
            full_model=TRUE,
            best=TRUE,
            Z=object$strata, ...)
}
