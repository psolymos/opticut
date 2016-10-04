## this returns a #habclasses x #spp matrix
bestpart.uncertainty <-
function (object, ...)
{
    sapply(lapply(object$uncertainty, bestpart.uncertainty1), rowMeans)
}
