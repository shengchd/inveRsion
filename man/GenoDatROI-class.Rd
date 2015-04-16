\name{GenoDatROI-class}
\Rdversion{1.1}
\docType{class}
\alias{GenoDatROI-class}
\alias{initialize,GenoDatROI-method}
\alias{show,GenoDatROI-method}
\alias{ROIGenoDat}

\title{Class "GenoDatROI"}
\description{ \code{GenoDat} defined within an region of interest}
\section{Objects from the Class}{
object are created with calls to \code{ROIGenoDat(objectGenoDat, ROI)}
}

\section{Slots}{
  \describe{
    \item{\code{genoDat}:}{\code{"matrix"} . Genotypes }
    \item{\code{lociPos}:}{ \code{"numeric"}. Genomic coordinates }
    \item{\code{alleleSum}:}{ \code{"matrix"}. Total number of variant alleles in the population per SNP}
    \item{\code{noMissCount}:}{\code{"matrix"}. Total number of subjects with no-missing values}
    \item{\code{ROI}:}{\code{"numeric"}. Region of interest. }
}
}
\section{Extends}{
Class \code{"\linkS4class{GenoDat}"}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "GenoDatROI")}}
    \item{show}{\code{signature(object = "GenoDatROI")}}
	 }
}

\usage{
ROIGenoDat(objectGenoDat, ROI)
}
\arguments{
  \item{objectGenoDat}{ \code{GenoDat} 
}
  \item{ROI}{ \code{ numeric}. Region of interest. 2-component vector that defines the limits of a chromosome segments where an inversions is thought to occur. 

}
}
\details{
\code{ROIGenoDat} is the constructor of the class.
}


\author{ Alejandro Caceres \email{acaceres@creal.cat}}


\seealso{
	\code{\linkS4class{GenoDat}} 
	}
\examples{
data(gDat)
gDatROI<-ROIGenoDat(gDat,ROI=c(1268,1847))
gDatROI
}
\keyword{classes}
