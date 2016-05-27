\name{dolina}
\alias{dolina}
\docType{data}
\title{
A comprehensive and micro-scale land snail data set from 16 dolines of the
Aggtelek Karst Area, Hungary.
}
\description{
Data set containing land snail counts as described in Kemecei et al. 2014.
}
\usage{data("dolina")}
\format{
A list with 3 elements: \code{xtab} is a sample x species matrix with counts,
\code{samp} is a data frame with sample level attributes,
\code{taxa} is a data frame with scientific names and families for the species.

Land snails were sampled during daylight hours between
16 and 18 of August, 2007. Samples were taken from four
microhabitat types (\code{dolina$samp$microhab}, \code{dolina$samp$mhab}):
litter (LI), trunks of live trees (TL), dead
wood (also known as coarse woody debris; DW), and rock (RO).
In each of the 16 dolina (\code{dolina$samp$dolina}),
seven samples were collected in the litter microhabitat
along a north-south transect. In the case of the other three
microhabitat types, samples were collected from three random
locations per microhabitat type in each dolina.
A total of 256 samples (\code{dolina$samp$sample}) were
collected, each consisting 2 sub-samples collected
by 2 sampling methods (\code{dolina$samp$method}):
litter samples (Q) and timed search (T).

One liter of litter samples including topsoil were collected
to be examined later in the laboratory. Litter samples were
collected adjacent to live wood, dead wood and rocks, and
not from the wood or rocks themselves. Litter samples in the
litter microhabitat were not collected near wood or rocks
(minimum distance of 2 meters). During 5 minutes per site of
time-restricted direct search we investigated microhabitats in
a 1 meter radius circle around the litter sample location, but
also including tree or rock surfaces.

The vertical zone (\code{dolina$samp$stratum},
bottom, middle or edge of the dolinas),
aspect of these sample locations (\code{dolina$samp$aspect}),
along with litter depth (\code{dolina$samp$lthick}, cm),
and litter moisture (\code{dolina$samp$lmoist},
scored on an ordinal scale: 1=dry,
2=fresh, 3=moist) were also recorded.

Distinction of live animals versus fresh
empty shells was not feasible due to the
method of sorting dry material and the delay in litter sample
processing, so these were combined and constituted the
'fresh' group. Whitened, disintegrating and broken shells
constituted the 'broken' group. This 'broken' group was
excluded from the data set.
}
\source{
Solymos et al. 2016 and Kemencei et al. 2014.
}
\references{
Kemencei, Z., Farkas, R., Pall-Gergely, B., Vilisics, F., Nagy, A., Hornung,
E. & Solymos, P. (2014): Microhabitat associations of land snails in
forested dolinas: implications for coarse filter conservation.
Community Ecology 15:180--186.

Solymos, P., Kemencei, Z. Pall-Gergely, B., Farkas, R., Vilisics, F., Nagy, A., Kisfali, M. & Hornung, E. (2016): Public data from the dolina project. Version 1.0. Zenodo, DOI: 10.5281/zenodo.53080
}
\examples{
data(dolina)
str(dolina)

## species richness by microhabitat
Richness <- rowSums(dolina$xtab > 0)
boxplot(Richness ~ mhab, dolina$samp, ylab="Species richness")
}
\keyword{datasets}
