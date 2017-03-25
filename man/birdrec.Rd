\name{birdrec}
\alias{birdrec}
\docType{data}
\title{
Bird Species Detections
}
\description{
Data set listing 156 species (mostly birds, few amphibians and mammals) detected
at 127 sites (367 point locations) in Alberta, Canada in 2015,
using autonomous reconding technology (ARU; Wildlife Acoustic Song Meter)
for sound recodings.
}
\usage{data("birdrec")}
\format{
A list with 3 elements with matching ordering:
\code{xtab} is a sample x species matrix with number of detections,
\code{samp} is a data frame with sample level attributes.
\code{taxa} is a data frame with species level attributes.

Multiple random recordings at each location were selected
according to a stratified random design
(based on combination of \code{TOY} and \code{TOD}).
These recordings were listened to by trained analysts
and species were identified based on auditory cues.

This data set lists detections from the first 1-minute segment of each
recording. Dates for the 3967 1-minute segments
range between 2015-03-31 and 2015-07-29.
Variables in \code{birdrec$samp} are the following:

\code{PKEY}: primary key for location/time combinations.

\code{POINT}: unique spatial location IDs, each point had
its own ARU unit.

\code{SITE}: site ID (1-4 ARU units deployed per site).

\code{YEAR}: year, 2015.

\code{MONTH}: month from 3 (March) to 7 (July).

\code{MDAY}: day of month, 1-31.

\code{HOUR}: 24-hour of day, values between 0-12.

\code{MINUTE}: minute, 0-59.

\code{YDAY}: ordinal day of the year, 89-209.

\code{RAIN}, \code{WIND}, \code{INDUSTRY}, \code{NOISE}:
level of rain, wind, indistrial noise, and background noise.
0 = no; 1 = light; 2 = moderate; 3 = heavy.

\code{MICROPHONE}:
Every recording contains a certain level of background static due to the
pre-amplifiers; however, problems, such as, electrostatic discharge on the
microphones, faulty wiring, poorly installed microphones and/or
missing microphones can occur causing excess static or dead channels.
0 = no microphone related issues;
1 = left microphone cuts out intermittently;
2 = right microphone cuts out intermittently;
3 = both microphones cut out intermittently ;
4 = left channel failed;
5 = right channel failed;
6 = both channels failed (no cases in the data set);
7 = left side extra static;
8 = right side extra static;
9 = both sides extra static;
10 = other issues;
11 = unbalanced channels.

\code{TOY}: time of year intervals used for stratified random
selection of dates. 8 intervals divied into 3 major units
(early, mid, and late breeding season; \code{YDAY} 140 and 180 were used
as threshold between the major units).

\code{TOD}: time of day, midnight (\code{HOUR} = 0) or
morning (\code{HOUR} > 0).

Variables in \code{birdrec$taxa} are the following:
FIXME

Methodology and metadata is described in ABMI (2016), and Lankau et al. (2015).
}
\source{
Alberta Biodiversity Monitoring Institute (ABMI, www.abmi.ca)
}
\references{
Alberta Biodiversity Monitoring Institute (ABMI), 2016.
Terrestrial field data collection protocols (abridged version) 2016-05-18.
Alberta Biodiversity Monitoring Institute; Edmonton, Alberta, Canada.

Lankau, H.E., MacPhail, A., Knaggs, M. & Bayne, E., 2015.
Acoustic recording analysis protocol.
Bioacoustic Unit, University of Alberta, and
Alberta Biodiversity Monitoring Institute; Edmonton, Alberta, Canada.
}
\examples{
data(birdrec)
str(birdrec)

aggregate(rowSums(birdrec$y),
    list(TOY=birdrec$x$TOY, TOD=birdrec$x$TOD), mean)
boxplot(rowSums(birdrec$y) ~ TOD + TOY, birdrec$x,
    col=c("gold", "tomato"), ylab="# detections")

#FIXME -- add more stuff here
}
\keyword{datasets}
