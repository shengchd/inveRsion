\name{scanInv}
\alias{scanInv}
\title{
Inversion scan 
}
\description{
This function scans a whole chromosome in search for inversion events. The scan is done by fitting an inversion model to all segments in the chromosome with fixed length size.
}
\usage{
scanInv(objectHaploCode, window, maxSteps = 30, geno = FALSE, saveRes = TRUE, saveBlocks=TRUE)
}

\arguments{
  \item{objectHaploCode}{Object of class \code{HaploCode} produced by the \code{codeHaplo} function.
}
  \item{window}{\code{numeric}, size of the window in mega-basis.
}
  \item{maxSteps}{\code{numeric}, maximum number of iteration in the EM algorithm for the inversion model
}
  \item{geno}{\code{logical}. Whether the original data is genotypes or phases haplotypes. 
}
  \item{saveRes}{\code{logical}. Whether results should be saves into file \code{invRes.RData} 
}
  \item{saveBlocks}{\code{logical}. Whether save blocks for each candidate break point.  
}

}
\details{
The function processes the haplotypes coded in \code{objectHaploCode}. If subsequent re-runs are requires for different window sizes, this object can be omitted. The function will thus search the local directory for previous results to speed up further scans. 
}
\value{
object of class \code{scan}
}
\author{
Alejandro Caceres \email{acaceres@creal.cat}}

\seealso{
\code{\linkS4class{HaploCode}} , \code{\linkS4class{scan}}   
}

\examples{
data(hapCode)
window<-0.5
scanRes<-scanInv(hapCode,window=window,saveRes=FALSE,geno=FALSE,saveBlocks=FALSE)
scanRes
plot(scanRes)
}
\keyword{ constructor}
