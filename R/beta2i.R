beta2i <- function(x, scale=1) {
    abs(tanh(x * scale))
}