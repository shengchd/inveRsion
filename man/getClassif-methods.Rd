\name{getClassif-methods}
\docType{methods}
\alias{getClassif-methods}
\alias{getClassif,inversionList-method}
\alias{getClassif}

\title{ Overall classification }
\description{
Classifies into inverted or non-inverted populations each chromosome in the sample. }
\section{Methods}{
\describe{
\item{\code{signature(object = "inversionList")}}{ for each of the inversions of the list, it returns the classification of each chromosome.}
}}

\usage{
getClassif(object, thBic, wROI, bin,geno,id)
}

\arguments{
  \item{object}{\code{inversionList}. List of inversions obtained from a chromosome scan. 
}
  \item{thBic}{ \code{numeric}. BIC threshold above which significant segments are chosen for the final classification
}
  \item{wROI}{\code{numeric}. ROI number from the list to select classification }
  \item{bin}{\code{logic}. Whether binary or continous classification is retrieved}
  \item{geno}{\code{logic}. Whether inversion genotype is retrieved}
  \item{id}{\code{character}. Vector of subject IDs}
}
\details{
The overall classification of chromosomes into inverted and non-inverted populations is given by the majority vote of the classifications obtained for each trial segment in the ROI, with BIC greater than \code{thBic}. 
}
\value{
\code{numeric}. Vector with values between 0 and 1 representing membership to the non-inverted and inverted population respectively.
}
\author{
Alejandro Caceres \email{acaceres@creal.cat}
}

\seealso{
\code{\linkS4class{inversionList}}
}
\examples{
data(invList)
r<-getClassif(invList)
head(r)

}

\keyword{methods}
