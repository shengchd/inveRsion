\name{gDat}
\alias{gDat}
\docType{data}
\title{
Sample data of class genoDat
}
\description{
Data set used to illustrate local haplotype coding performed with \code{codeHaplo} 
}
\usage{data(gDat)}
\format{
  The format is:
Formal class 'GenoDat' [package "inveRsion"] with 4 slots
  ..@ genoDat    : int [1:9, 1:10] 0 1 0 1 0 0 1 0 0 1 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:9] "V1" "V2" "V3" "V4" ...
  .. .. ..$ : NULL
  ..@ lociPos    : Named int [1:10] 959 1268 1393 1467 1531 1761 1847 1987 2006 2030
  .. ..- attr(*, "names")= chr [1:10] "V1" "V2" "V3" "V4" ...
  ..@ alleleSum  : num [1:10, 1] 3 1 8 3 9 1 1 9 1 1
  ..@ noMissCount: num [1:10, 1] 9 9 9 9 9 9 9 9 9 9
}

\examples{
data(gDat)
gDat
plot(gDat)
}
\keyword{datasets}
