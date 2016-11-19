\name{occolors}
\alias{occolors}
\alias{col2gray}
\title{
Color palettes for the opticut package
}
\description{
A convenient way of setting color palettes for the opticut package.
}
\usage{
occolors(theme)
col2gray(col, method="BT.709")
}
\arguments{
  \item{theme}{
character value, character vector, or a function used to interpolate the colors.
The built-in values are \code{"br"} (blue-red divergent palette, colorblind safe),
\code{"gr"} (green-red divergent palette), \code{"bw"} (black and white: grayscale converted \code{"br"} settings).
See \code{\link{colorRampPalette}}, \code{\link{gray}} and the Examples.
}
  \item{col}{
vector of color specification as described on the help page
for the \code{\link{col2rgb}} function. This is converted to grayscale.
}
  \item{method}{
character, the method used for grayscale conversion.
}
}
\description{
Grayscale conversion methods calculate gray levels based on
red (R), green (G), and blue (B) color channels:
\item{"BT.709"}{0.2126 * R + 0.7152 * G + 0.0722 * B, luma correction following the  ITU-R BT.709 recommendation}
\item{"BT.601"}{0.299 * R + 0.587 * G + 0.114 * B, luma correction following the  ITU-R BT.601 recommendation}
\item{"desaturate"}{(max(R, G, B) + min(R, G, B)) / 2, also called lightness}
\item{"average"}{(R + G + B) / 3}
\item{"maximum"}{max(R, G, B)}
\item{"minimum"}{min(R, G, B)}
\item{"red"}{R}
\item{"green"}{G}
\item{"blue"}{B}
}
\value{
\code{occolors} returns a function, see \code{\link{colorRampPalette}}.

\code{col2gray} returns a vector of grey colors based on
the conversion \code{method} and \code{\link{gray}}.
}
\author{
Peter Solymos <solymos@ualberta.ca>

Hexadecimal values for the built-in palettes are taken from
\url{http://colorbrewer2.org/}.
}
\seealso{
\code{\link{colorRampPalette}} for a general description of palettes.

\code{\link{ocoptions}} for setting the color theme options in the
opticut package.
}
\examples{
## using palettes
plot(1:100, rep(2, 100), pch = 15,
    ylim = c(0, 21), axes = FALSE, ann = FALSE,
    col = occolors()(100)) # default 'bg'
text(50, 1, "theme = 'br'")
points(1:100, rep(5, 100), pch = 15,
    col=occolors("gr")(100))
text(50, 4, "theme = 'gr'")
points(1:100, rep(8, 100), pch = 15,
    col=occolors("bw")(100))
text(50, 7, "theme = 'bw'")
points(1:100, rep(11, 100), pch = 15,
    col=occolors(terrain.colors)(100))
text(50, 10, "theme = terrain.colors")
points(1:100, rep(14, 100), pch = 15,
    col=occolors(c("purple", "pink", "orange"))(100))
text(50, 13, "theme = c('purple', 'pink', 'orange')")
points(1:100, rep(17, 100), pch = 15,
    col=occolors(c("#a6611a", "#ffffbf", "#018571"))(100))
text(50, 16, "theme = c('#a6611a', '#ffffbf', '#018571')")
points(1:100, rep(20, 100), pch = 15,
    col=occolors(c("#7b3294", "#ffffbf", "#008837"))(100))
text(50, 19, "theme = c('#7b3294', '#ffffbf', '#008837')")

## grayscale conversions
n <- 25
col <- occolors("br")(n)
method <- c("BT.709", "BT.601",
    "desaturate", "average", "maximum", "minimum",
    "red", "green", "blue")
plot(0, type="n", ann=FALSE, axes=FALSE,
    xlim=c(0, n), ylim=c(3*length(method), 0))
for (j in 1:length(method)) {
    for (i in 1:n) {
        polygon(c(i-1, i, i, i-1), c(0, 0, 1, 1)+((j-1)*3),
            col=col[i], border=col[i])
        polygon(c(i-1, i, i, i-1), c(1, 1, 2, 2)+((j-1)*3),
            col=col2gray(col[i], method=method[j]),
            border=col2gray(col[i], method=method[j]))
        text(n/2, 1+((j-1)*3), method[j])
    }
}
}
\keyword{ utilities }
