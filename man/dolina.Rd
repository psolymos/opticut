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
  The format is:
List of 3
 $ xtab: num [1:512, 1:42] 0 0 0 0 0 0 0 0 0 0 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:512] "10A1Q" "10A1T" "10A2Q" "10A2T" ...
  .. ..$ : chr [1:42] "aacu" "amin" "apol" "apur" ...
 $ samp:'data.frame':	512 obs. of  9 variables:
  ..$ sample  : Factor w/ 256 levels "10A1","10A2",..: 1 1 2 2 3 3 4 4 5 5 ...
  ..$ dolina  : int [1:512] 10 10 10 10 10 10 10 10 10 10 ...
  ..$ microhab: Factor w/ 4 levels "dead.wood","litter",..: 2 2 2 2 2 2 2 2 2 2 ...
  ..$ mhab    : Factor w/ 4 levels "A","H","L","R": 1 1 1 1 1 1 1 1 1 1 ...
  ..$ method  : Factor w/ 2 levels "Q","T": 1 2 1 2 1 2 1 2 1 2 ...
  ..$ aspect  : Factor w/ 5 levels "eastern","flat",..: 4 4 4 4 4 4 2 2 3 3 ...
  ..$ stratum : Factor w/ 4 levels "1bottom","2middle",..: 4 4 3 3 2 2 1 1 2 2 ...
  ..$ lmoist  : num [1:512] 1 1 1 1 1 1 1.5 1.5 1 1 ...
  ..$ lthick  : num [1:512] 2 2 2.5 2.5 3 3 0.5 0.5 1.5 1.5 ...
 $ taxa:'data.frame':	42 obs. of  2 variables:
  ..$ scientific.name: Factor w/ 51 levels "Acanthinula aculeata",..: 1 3 39 4 5 6 13 10 16 14 ...
  ..$ author         : Factor w/ 28 levels "(Alder, 1830)",..: 12 19 20 1 11 7 16 25 4 11 ...
}
\details{
A list with 3 elements: \code{xtab} is a sample x species matrix with counts,
\code{samp} is a data frame with sample level attributes,
\code{taxa} is a data frame with scientifc names and authorities for the species.

Land snails were sampled during daylight hours between
16 and 18 of August, 2007. Samples were taken from four
microhabitat types (\code{dolina$samp$microhab}, \code{dolina$samp$mhab}):
litter (A), trunks of live trees (L), dead
wood (also known as coarse woody debris; H), and rock (R).
In each of the 16 dolina (\code{dolina$samp$dolina}),
seven samples were collected in the litter microhabitat
along a north-south transect. In the case of the other three
microhabitat types, samples were collected from three ran-
dom locations per microhabitat type in each dolina.
A total of 256 samples (\code{dolina$samp$sample}) were
collected, each consisting 2 subsamples collected
by 2 sampling methods (\code{dolina$samp$method}):
litter samples (Q) and timed search (T).

One liter of litter samples including topsoil were collected
to be examined later in the laboratory. Litter samples were
collected adjacent to live wood, dead wood and rocks, and
not from the wood or rocks themselves. Litter samples in the
litter microhabitat were not collected near wood or rocks (mi-
nimum distance of 2 meters). During 5 minutes per site of
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
Kemencei et al. 2014.
}
\references{
Kemencei, Z., Farkas, R., Pall-Gergely, B., Vilisics, F., Nagy, A., Hornung, E. & Solymos, P. (2014): Microhabitat associations of land snails in forested dolinas: implications for coarse filter conservation. Community Ecology 15:180--186.
}
\examples{
data(dolina)
str(dolina)
}
\keyword{datasets}
