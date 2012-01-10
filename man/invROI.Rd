\name{invROI}
\alias{invROI}
\title{
Inversion scan 
}
\description{
This function scans a region of interest to detect inversion events. The scan is done by fitting an inversion model to all possible segments within the ROI. 
}
\usage{
invROI(hapCode,geno,ROI,saveRes,thBic)
}

\arguments{
  \item{hapCode}{Object of class \code{HaploCode} produced by the \code{codeHaplo} function.
}
 \item{geno}{\code{logical}. Whether the original data is genotypes or phases haplotypes. 
}
 \item{ROI}{\code{vector}.  4D vector with the limits for the left and break points respectively.
}
  \item{saveRes}{\code{logical}. Whether results should be saves into file \code{invRes.RData} 
}
  \item{thBic}{\code{numeric}. BIC threshold above which candidate segments are chosen to infer inversion status of each chromosome in the sample. 
  }

}
\details{
The function processes the haplotypes coded in \code{objectHaploCode}. It is an intensive search of all possible segments, therefore the number of computations is large if the limits for the break points are not tight.}
\value{
object of class \code{inversionList}
}
\author{
Alejandro Caceres \email{acaceres@creal.cat}}

\seealso{
\code{\linkS4class{HaploCode}} , \code{\linkS4class{inversionList}}   
}

\examples{
data(hapCode)
ROI<-c(0.746,0.754,1.245,1.255)
iROI<-invROI(hapCode,ROI=ROI,saveRes=FALSE,geno=FALSE)
iROI
}
\keyword{ constructor}
