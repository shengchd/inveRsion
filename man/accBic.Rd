\name{accBic}
\alias{accBic}
\alias{accBic-methods}
\alias{accBic,inversionList-method}

\title{accBic computes "accuracy" from "inversionList"}
\description{
\code{accBic} computes the accuracy of the classification of chromosomes into previously known inverted and non-inverted populations. The classification is obtained from a majority vote of the classifications produced by the trial segment models whose BIC is greater than a given threshold. 

}
\usage{
accBic(object, mem, classFile, nsub, npoints, geno, wROI)
}
\arguments{
  \item{object}{ of class \code{inversionList}
}
  \item{mem}{
vector with the numbering of chromosomes known to have the inversion 
}
  \item{classFile}{an alternative to \code{mem}, it passes the file name containing the numbering of chromosomes known to have the inversion. 
}
  \item{nsub}{
total number of subjects (= 2* total number of chromosomes)
}
  \item{npoints}{
number of BIC threshold between 0 and max (BIC) for which the accuracy is to be computed
}
  \item{geno}{ 
  whether the accuracy is assessed for inversion genotype or inversion allele (phased data).  
  
}
  \item{wROI}{
integer indicating the ROI number to be used. The total number of ROIs are the total number of components in the \code{object} list.
}
}
\value{
\item{\code{accuracy}}{ object of class \code{accuracy}}
}
\author{
Alejandro Caceres \email{acaceres@creal.cat}}
\seealso{ \code{\linkS4class{inversionList}}, \code{\linkS4class{accuracy}} }

\examples{
data(invList)
memFile <- system.file("extdata", "mem.txt", package = "inveRsion")
ac <- accBic(invList, classFile = memFile, nsub = 1000, npoints = 10)
plot(ac)
}

\keyword{methods}
