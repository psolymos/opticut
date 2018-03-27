\name{warblers}
\alias{warblers}
\docType{data}
\title{
Warblers Data Set
}
\description{
Five species of warblers were studies to determine the factors
controlling the species abundances and competition (MacArthur 1958).
}
\usage{data("warblers")}
\format{
A list with 3 elements: \code{xtab} is a list of sample x species matrices
(\code{sec_prc}: percentages of total number of seconds of observations,
\code{num_prc}: percentages of total number of observations in seconds,
\code{sec_cnt}: counts based on percentages and totals of seconds,
\code{num_cnt}: counts based on percentages and totals),
\code{samp} is a data frame with sample level attributes such as height
(6 is the base, 1 is the top of trees, most trees were 50--60 feet tall)
and depth of branches in the canopy
(B: bare or lichen- covered base, M: middle zone of old needles,
T: terminal zone of new, less than 1.5 years old, needles or buds) ,
\code{taxa} is a data frame with scientific and common names for the species.
}
\source{
MacArthur 1958.
}
\references{
MacArthur, R. H., 1958.
Population ecology of some warblers of northeastern coniferous forests.
Ecology 39:599--619.
<doi:10.2307/1931600>
}
\examples{
data(warblers)
str(warblers)

warbh <- mefa4::groupSums(warblers$xtab$num_cnt, 1, warblers$samp$height)
warbd <- mefa4::groupSums(warblers$xtab$num_cnt, 1, warblers$samp$depth)
op <- par(mfrow=c(1,2))
matplot(rownames(warbh), warbh, type="b", xlab="height")
matplot(warbd, type="b", axes=FALSE, xlab="depth")
box()
axis(2)
axis(1, 1:3, rownames(warbd))
par(op)
}
\keyword{datasets}
