summary.lorenz <-
function (object, ...) {
    out <- attr(object, "summary")
    class(out) <- "summary.lorenz"
    out
}
