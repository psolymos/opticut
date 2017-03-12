strata.uncertainty <-
function (object, ...) {
    if (inherits(object, "uncertainty_opti"))
        strata.opticut(object, ...)
    if (inherits(object, "uncertainty_multi"))
        strata.multicut(object, ...)
}
