\name{scanRes}
\alias{scanRes}
\docType{data}
\title{
Sample data set of class scan
}
\description{
Sample set that illustrates the result of scanning a chromosome with segments of fix window size. It fits the inversion model for each segment and stores its results.  
}
\usage{data(scanRes)}
\format{
  The format is:
Formal class 'scan' [package "inveRsion"] with 10 slots
  ..@ leftBP  : num [1:189, 1] 0.603 0.604 0.606 0.609 0.612 ...
  ..@ rightBP : num [1:189, 1] 1.15 1.15 1.15 1.15 1.15 ...
  ..@ leftBP2 : num [1:189, 1] 0.604 0.606 0.608 0.609 0.612 ...
  ..@ rightBP2: num [1:189, 1] 1.15 1.15 1.15 1.15 1.15 ...
  ..@ LogLike : num [1:189, 1] 27.2 45.4 89.8 108.2 113.3 ...
  ..@ prob    : num [1:189, 1] 0.996 0.987 0.947 0.936 0.946 ...
  ..@ ent     : num [1:189, 1] 2.41 2.33 2.78 2.87 2.92 ...
  ..@ entTh   : num [1:189, 1] 29 23 26 30 24 28 20 21 24 20 ...
  ..@ bic     : num [1:189, 1] -193.2 -129.4 -107.8 -119.8 -69.1 ...
  ..@ window  : num 0.5
}
\details{
The object is constructed with the function \code{scanInv} and used as an input for \code{listInv}.
}
\examples{
data(scanRes)
plot(scanRes) 
}
\keyword{datasets}
