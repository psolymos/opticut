\name{beta2i}
\alias{beta2i}
\title{
Scaling for the Indicator Potential
}
\description{
Transformation of estimated contrasts to indicator potential.
}
\usage{
beta2i(x, scale = 1)
}
\arguments{
  \item{x}{
numeric, real valued coefficients.
}
  \item{scale}{
numeric, scaling constant.
}
}
\value{
Returns a numeric vector (I = abs(tanh(x * scale))).
}
\author{
Peter Solymos <solymos@ualberta.ca>
}
\seealso{
\code{\link{opticut}} and \code{\link{multicut}}
use the scaled I values as indicator potential.

\code{\link{ocoptions}} for setting value for the default scaling factor.
}
\examples{
x <- seq(-5, 5, 0.1)
Col <- occolors(c("red", "blue"))(10)
plot(x, beta2i(x), type = "n")
s <- seq(1, 0.1, -0.1)
for (i in 1:10) {
    lines(x, beta2i(x, scale = s[i]), col = Col[i])
    text(1.5 - 0.2, beta2i(1.5, scale = s[i]), s[i], col = Col[i])
}
}
\keyword{ utilities }
