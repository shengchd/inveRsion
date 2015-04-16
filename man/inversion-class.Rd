\name{inversion-class}
\Rdversion{1.1}
\docType{class}
\alias{inversion-class}

\title{Class "inversion" }
\description{Internal class (not to be called by the user) that retrieves the output on the inversion model run in an ROI. \code{inversionList} is a list of objects of class \code{inversion}}
\section{Objects from the Class}{
Objects are created by calls to \code{listInv(object, hapCode, geno, ROI, saveRes, thBic, all)}}
\section{Slots}{
  \describe{
    \item{\code{classification}:}{Object of class \code{"vector"} overall classification (majority vote) of each chromosome for all the trail segments in the ROI }
    \item{\code{leftBP}:}{Object of class \code{"vector"} left brake-points for each of the trial segments}
    \item{\code{rightBP}:}{Object of class \code{"vector"} right brake-points for each of the trial segments}
    \item{\code{bic}:}{Object of class \code{"vector"} BIC for the models on each trail segment }
    \item{\code{intLeftBP}:}{Object of class \code{"vector"}  interval for the left brake-point, in the ROI}
    \item{\code{intRightBP}:}{Object of class \code{"vector"} interval for the right-brake point in the ROI}
    \item{\code{invFreq}:}{Object of class \code{"numeric"} overall inversion frequency }
    \item{\code{RR}:}{Object of class \code{"list"} classification of given by each of the trail segments. }
  }
}
\author{ Alejandro Caceres \email{acaceres@creal.cat}}

\seealso{ \code{\link{listInv}}, \code{\linkS4class{inversionList}}} 
\keyword{classes}
