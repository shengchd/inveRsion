\name{hapCode}
\alias{hapCode}
\docType{data}
\title{
Sample data set of class HaploCode
}
\description{
Illustrative data set, with local haplotypes encoded, to be used as input of \code{scanInv}.
}
\usage{data(hapCode)}
\format{
  The format is:
Formal class 'HaploCode' [package "inveRsion"] with 3 slots
  ..@ haploCode: num [1:2000, 1:583] 8 8 37 37 8 37 8 8 37 37 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:583] "0.602976-0.604061" "0.604061-0.605972" "0.60602-0.608417" "0.608668-0.608855" ...
  ..@ blockSize: num 3
  ..@ minAllele: num 0.3
}
\examples{
data(hapCode)
hapCode
}
\keyword{datasets}
