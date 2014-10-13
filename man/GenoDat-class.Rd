\name{GenoDat-class}
\Rdversion{1.1}
\docType{class}
\alias{GenoDat-class}
\alias{plot,GenoDat-method}
\alias{show,GenoDat-method}

\title{Class "GenoDat"}
\description{Object that handles genotype data}

\section{Objects from the Class}{
Objects can be created by calls of the form 
\code{setUpGenoDatFile} or \code{setUpGenoDatSNPmat} .
}

\section{Slots}{
  \describe{
    \item{\code{genoDat}:}{\code{"matrix"} . Genotypes }
    \item{\code{lociPos}:}{ \code{"numeric"}. Genomic coordinates }
    \item{\code{alleleSum}:}{ \code{"matrix"}. Total number of variant alleles in the population per SNP}
    \item{\code{noMissCount}:}{\code{"matrix"}. Total number of subjects with no-missing values}
}
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "GenoDat")}: Plots minor allele frequency }
    \item{show}{\code{signature(object = "GenoDat")}: shows data summary}
	 }
}

 
\author{ Alejandro Caceres \email{acaceres@creal.cat} }

\seealso{
 \code{\link{setUpGenoDatFile}},  \code{\link{setUpGenoDatSNPmat}} 

}

\examples{
data(gDat)
gDat
}
\keyword{classes}
