\name{getROIs-methods}
\docType{methods}
\alias{getROIs-methods}
\alias{getROIs,scan-method}
\alias{getROIs}
\title{Extracts regions of possible inversion events from "scan"}
\description{
lists the regions of interest in a matrix, given by the overlapping of significant trial segments (of fixed window seize) that may be part of an inversion event. 
}
\section{Methods}{
\describe{
\item{\code{signature(object = "scan")}}{ list of regions of interest }
}}
\usage{
getROIs(object, thBic)
}
\arguments{
  \item{object}{\code{scan}. Cromosome scanned for inversions with trial segments of fixed window size.
}
  \item{thBic}{
BIC threshold above which overlapping segments are considered for the definition of each ROI 
}
}
\details{
ROIs are defined as overlapping trial segments with BIC greater than \code{thBIC}. The output is a matrix for which each row is one ROI. The first two columns give intervals defining the left brake-points and the two subsequent columns are the intervals for the right brake-points.   
ROIs are given in mega-basis.

}
\value{
\code{matrix}
}
\author{
Alejandro Caceres \email{acaceres@creal.cat}}

\seealso{
\code{\linkS4class{scan}}
}
\examples{
data(scanRes)
ROI <- getROIs(scanRes, thBic = 0)
ROI
}





\keyword{methods}
