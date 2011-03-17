\name{setUpGenoDatFile}
\alias{setUpGenoDatFile}
\alias{setUpGenoDatSNPmat}

\title{
Loads genotype data onto R 
}
\description{
Loads onto an R session genotype data from text files or PLINK files. 

}
\usage{
setUpGenoDatFile(file = "GenoDat.txt", saveRes = FALSE, sortMinor = TRUE)
setUpGenoDatSNPmat(Chr, Geno, Annot, saveRes = FALSE, saveGeno = FALSE)
}
\arguments{
  \item{file}{ \code{character}. File path with genotype information 
}
  \item{saveRes}{ \code{logical}. Wether results should be saved into file \code{ gDat.RData}

}
  \item{sortMinor}{ \code{logical}. Whether genotypes should be sorted by minor allele frequency.
}
 
 \item{Chr}{ \code{numeric}. Chromosome number}
 \item{Geno}{\code{snpMatrix}. Matrix of \code{raw} with genotype data}
 \item{Annot}{ \code{numeric}. Annotation information read from an \code{.bim} file}
 \item{saveGeno}{\code{logical}. Wether \code{.txt} file should be saved with genotype information}  

}
\value{
\code{GenoDat} object
}
\author{
Alejandro Caceres \email{acaceres@creal.cat}}


\seealso{
 \code{\linkS4class{GenoDat}}

}
\examples{
gen <- system.file("extdata", "genotypes.txt", package = "inveRsion")
gDat <-setUpGenoDatFile(file=gen,sortMinor=TRUE,saveRes=FALSE)
gDat
}

\keyword{ class}
