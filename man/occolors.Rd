\name{occolors}
\alias{occolors}
\title{
Color palettes for the opticut package
}
\description{
A convenient way of settings color palettes for the opticut package.
}
\usage{
occolors(theme)
}
\arguments{
  \item{theme}{
character value, character vector, or a function used to interpolate the colors,
The built-in values are \code{"br"} (blue-red divergent palette, colorblind safe),
\code{"br"} (green-red divergent palette), \code{"bw"} (black and white).
See \code{\link{colorRampPalette}} and the Examples.
}
}
\value{
Returns a function, see \code{\link{colorRampPalette}}.
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
}
\keyword{ utilities }
