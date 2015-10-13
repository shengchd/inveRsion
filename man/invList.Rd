\name{invList}
\alias{invList}
\docType{data}
\title{
Sample data of class inversionList}
\description{
list of objects of class \code{inversion}, each of which contains the information of overlapping segments that cover candidate inversions within a chromosome. 
}
\usage{data(invList)}
\format{
  The format is:
Formal class 'inversionList' [package "inveRsion"] with 1 slots
  ..@ results:List of 1
  .. ..$ :Formal class 'inversion' [package "inveRsion"] with 8 slots
 }
\details{
The object is constructed with the function \code{listInv}
}
\examples{
data(invList)
invList
}
\keyword{datasets}
