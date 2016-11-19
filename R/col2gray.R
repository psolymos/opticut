## http://www.tannerhelland.com/3643/grayscale-image-algorithm-vb6/
col2gray <-
function(col, method="BT.709")
{
    method <- match.arg(method, c("BT.709", "BT.601",
        "desaturate", "average", "maximum", "minimum",
        "red", "green", "blue"))
    col <- col2rgb(col) / 255
    if (method == "BT.709")
        out <- 0.2126*col["red",] + 0.7152*col["green",] + 0.0722*col["blue",]
    if (method == "BT.601")
        out <- 0.299*col["red",] + 0.587*col["green",] + 0.114*col["blue",]
    if (method == "desaturate")
        out <- (apply(col, 2, max) + apply(col, 2, min)) / 2
    if (method == "average")
        out <- colMeans(col)
    if (method == "maximum")
        out <- apply(col, 2, max)
    if (method == "minimum")
        out <- apply(col, 2, min)
    if (method == "red")
        out <- col["red",]
    if (method == "green")
        out <- col["green",]
    if (method == "blue")
        out <- col["blue",]
    gray(out)
}
