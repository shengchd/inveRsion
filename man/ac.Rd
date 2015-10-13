\name{ac}
\alias{ac}
\docType{data}
\title{
Sample data set of class accuracy
}
\description{
Illustrative output of function \code{accBic}. It stores the accuracy of classification 
of each chromosome into the inverted population and the frequency of the inversion as a 
for a range of Bic thresholds.
}
\usage{data(hapCode)}
\format{
  The format is:
Formal class '"accuracy"' [package "inveRsion"] with 1 slots
  ..@ out: num 
         bicInt   prob     ac
 [1,]    0.0000 0.3180 0.9045
 [2,]  142.4209 0.3745 0.9610
 [3,]  284.8419 0.3945 0.9810
 [4,]  427.2628 0.4130 0.9995
 [5,]  569.6838 0.4130 0.9995
}
\examples{
data(ac)
ac
}
\keyword{datasets}
