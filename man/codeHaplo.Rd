\name{codeHaplo}
\alias{codeHaplo}
\title{
Codes haplotypes into decimal integers
}
\description{
The function labels the haplotypes of size \code{blockSize} around each candidate brake-point. For labeling genotype data, the function takes objects of class \code{genoDat} as main argument. For phased data, this argument should be ignored and a file name passed instead. 

}
\usage{
codeHaplo(objectGenoDat, blockSize, minAllele, saveRes = TRUE, file = NULL, ROI,intSNP=FALSE,phasing="inversion&BP")
}

\arguments{
  \item{objectGenoDat}{ \code{Genodat} object; if phased data then provide \code{file} instead}
  \item{blockSize}{\code{numeric}. number of SNPs flanking each side of each candidate brake-point. Default value 3}
  \item{minAllele}{ \code{numeric}. minimum allele frequency for each probe to be considered as a candidate brake-point. Default value 0.1}
  \item{saveRes}{
\code{logical}. Whether results should be saved in file \code{hapCode.RData}}
  \item{file}{\code{character}. File name with phased data }
  \item{ROI}{ \code{numeric}. 2-vector specification passes a chromosome segment to be encoded. 4-vector specification passes the region of interest for the left brake-point (ROI[1] and ROI[2]) and the right brake-point (ROI[3] and ROI[4])}
  \item{intSNP}{\code{logical}. whether build flanking blocks of uniform SNP density across the genome.} 
  \item{phasing}{\code{character}. Either c("ivenrsion&BP","forward", "forward&inverted"), defines the phasing strategy between all four flanking blocks of each break poin. ivenrsion&BP phases the internal flanking block first, then the external blocks to the internal bloks and match them. "forward" uses the blocks in the sequence assumming the forward population; in the inversion model, the backward population is obtained inverting the internal blocks. "forward&inverted" phases the forward and inverted population independently.}

}
\details{
When \code{setUpGenodat} is passed, the coding first computes the local haplotypes for each candidate-brake point form the genotype data and then encodes each haplotype into a decimal integer. The local haplotypes are computed with \code{haplo.em} form \code{haplo.stats} and assigns those with highest posterior probability to each chromosome. In the case of phased data, passed through \code{file}, no local haplotyping is needed and only the labeling is performed.
}
\value{
Object of class \code{HaploCode}   
}

\author{
Alejandro Caceres \email{acaceres@creal.cat}
}
\references{
http://mayoresearch.mayo.edu/mayo/research/schaid_lab/software.cfm
}
\seealso{
 \code{\linkS4class{GenoDat}},  \code{\linkS4class{HaploCode}}, 

}
\examples{
data(gDat)
hapCode <-codeHaplo(gDat,blockSize=3,minAllele=0.3,saveRes=FALSE)
hapCode
}
\keyword{constructor}
