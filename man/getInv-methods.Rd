\name{getInv-methods}
\docType{methods}
\alias{getInv-methods}
\alias{getInv,scan-method}
\alias{getInv}
\title{gets "scan" into a matrix}
\description{
Each row of the matrix represents a trial segment of fixed window size, for which the inversion model has been fit. It lists the left and right brake-points and output of the fitting: Log-likelihood ratio, probability of inversion, entropy, BIC (Bayes Information Criterion) and number of haplotypes.  
}
\section{Methods}{
\describe{
\item{\code{signature(object = "scan")}}{returns matrix with output of inversion model for each trial segment}
}}

\usage{
getInv(object, thBic, rnd, Like)
}

\arguments{
  \item{object}{ \code{scan}. Cromosome scanned for inversions with trial segments of fixed window size.
  
}
  \item{thBic}{\code{numeric}. BIC threshold above which data is retrieved.

}
  \item{rnd}{\code{logic}. Whether round matrix elements.

}
  \item{Like}{\code{numeric}. Log-likelihood ratio threshold above which data is retrieved.

}
}
\details{
Matrix with output of \code{scanInv}. Each row corresponds to a trail segments with given brake points and significance measures for the inversion model. 
}
\value{
\code{matrix}. 
}
\examples{
data(scanRes)
a<-getInv(scanRes,thBic=2500)
a
}

\keyword{methods}
