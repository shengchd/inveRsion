\name{inversionList-class}
\Rdversion{1.1}
\docType{class}
\alias{inversionList-class}
\alias{plot,inversionList-method}
\alias{show,inversionList-method}

\title{Class "inversionList"}
\description{Lists output of the inversion model for each region of interest specified.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{listInv(object, hapCode, geno, ROI, saveRes, thBic, all)} and are lists of \code{inversion} objects. 
}
\section{Slots}{
  \describe{
    \item{\code{results}:}{Object of class \code{"list"} list of \code{inversion} objects}
  }
}
\section{Methods}{
  \describe{
    \item{accBic}{\code{signature(object = "inversionList")}: Computes accuracy for chromosome classification when known }
    \item{getClassif}{\code{signature(object = "inversionList")}: Extracts classification for each chromosome }
    \item{plot}{\code{signature(x = "inversionList")}: Plots BIC values for trial segments}
    \item{show}{\code{signature(object = "inversionList")}: shows object}
	 }
}
\author{Alejandro Caceres \email{acaceres@creal.cat} }

\seealso{ \code{\link{listInv}}, \code{\linkS4class{HaploCode}} 
}
\examples{
data(invList)
invList
}
\keyword{classes}
