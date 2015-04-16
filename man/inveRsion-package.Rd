\name{inveRsion-package}
\alias{inveRsion-package}
\alias{inveRsion}
\docType{package}
\title{
Detection of Genetic inversions using SNP-array data}
\description{
This package scans the whole genome in search of inversion events. Input data can be genotypes or phased haplotypes. It computes regions where inversions are probed by trial segments of fixed length. An inversion model is fit at each trial segment and significance measures, like Bayes Information Criterion, give evidence of segments belonging to the inversion event. Methods are implemented to identify the complete inversion segment and to classify the chromosomes in the sample as inverted or not. 
}
\details{
\tabular{ll}{
Package: \tab inveRsion\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-11-12\cr
License: \tab GPL version 2 or newer \cr
LazyLoad: \tab yes\cr
Depends: \tab methods, haplo.stats\cr
}
The package is designed as a stream analysis of a sequence of procedures: \code{setUpGenoDatFile} loads the genotype data onto R; \code{HaploCode} performs local haplotyping around the candidate brake points of the inversion; \code{sanInv} takes trial segments of fixed window size and fits the inversion model, sweeping the whole genome; and \code{listInv} summarizes the inversion events. 
}
\author{
Alejandro Caceres
Maintainer: Alejandro Caceres <acaceres@creal.cat>
}
\references{
A Caceres el al. Detection of genetic inversions with SNP-array data, manuscript in preparation. 

SS Sindi and BJ Raphael, Identification and frequency estimation of inversion polymorphisms from haplotype data, J Comput Biol. 2010 Mar;17(3):517-31.

PF O'Reilly et al., invertFREGENE: software for simulating inversions in population genetic data, Bioinformatics. 2010 Mar 15;26(6):838-40. Epub 2010 Jan 26.

}
\keyword{ package }
\examples{
#vignette("inveRsionMan")
}
