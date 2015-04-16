\name{listInv-methods}
\docType{methods}
\alias{listInv-methods}
\alias{listInv,scan-method}
\alias{listInv}

\title{ Constructor of class inversionList}
\description{
Determines the full inversion sequence from overlapping trial segments of fixed window sized. It computes the limits of the inversion, the population frequency and retrieves the maximum BIC across the trial segments.   
}
\section{Methods}{
\describe{

\item{\code{signature(object = "scan")}}{ }
}}

\usage{
listInv(object, hapCode, geno, ROI, saveRes, thBic, all, saveROI)
}
\arguments{
  \item{object}{ \code{scan}. Results of scanning the genome for an inversion with a fixed window size; output of \code{scanInv}.
}
  \item{hapCode}{\code{HaploCode}. Object with the result of coding haplotypes for each candidate brake-point; output of \code{codeHaplo}
}
  \item{geno}{\code{logical}. Whether original data is genotypes of phased haplotypes.  
}
  \item{ROI}{\code{numeric}. 2-vector specification passes the chromosome segment to be analyzed. 4-vector specification passes the region of interest for the left brake-point (ROI[1] and ROI[2]) and the right brake-point (ROI[3] and ROI[4]). ROI should be specified in mega-basis units.   
}
  \item{saveRes}{\code{logical}.Whether results should be saved into file \code{list Inv.RData}
}
  \item{thBic}{\code{numeric}. BIC threshold above which trial segments are selected.
  }
  \item{all}{\code{logical}. Whether recomputing within the ROI should be done for all possible segment sizes or the window size in \code{scan}.
}
 \item{saveROI}{\code{logical}. Whether saving the blocks for the candidate break-points for all the ROIs.
}

}
\details{
\code{listInv} is both a method for class \code{scan} and constructor of \code{inversionList}. It re-runs the inversion model within the ROIs found in the previous scan. However, it is also possible to explicitly pass the ROIs defined by the user. The re-run is done with the same window size of the scan, which is convenient if enough significant trial segments where found within the inversion segment. If a more detailed re-run is needed set \code{all=TRUE}. This computes the model for all possible trial segments of any length within the ROI. Note that for high SNP density this can be computational intensive.

}

\value{
\code{inversionList}
}

\author{
Alejandro Caceres \email{acaceres@creal.cat}
}

\seealso{
\code{\linkS4class{HaploCode}},\code{\linkS4class{scan}},\code{\linkS4class{inversionList}}
}

\examples{
data(scanRes)
data(hapCode)
invList<-listInv(scanRes,hapCode=hapCode,geno=FALSE,saveRes=FALSE,all=FALSE,thBic=0, saveROI=FALSE)
invList
}

\keyword{methods}

