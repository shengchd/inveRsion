\name{inveRsion-internal}
\title{Internal inveRsion objects}
\alias{.First.lib}
\alias{.Last.lib}
\alias{getb1b2}
\alias{int2bit}
\alias{iterateInversionModel}
\alias{writeFilegenoDat}
\alias{callEncode}
\alias{datChr}
\alias{flipGeno}
\alias{encodeGeno}
\alias{encodeGenoAcross}
\alias{encodeHaplo}
\alias{manual}
\alias{getblockLocB12}


\usage{
getb1b2(BrakePoint, objectHaploCode)
int2bit(Num, BlockSize)
iterateInversionModel(b1b2,window,BlockSize,ls,maxSteps,geno,Coor,Coor2,candidatePoints,getInvclass) 
writeFileGenoDat(gg, allesum = FALSE, caco)
callEncode(object, BlockSize, MinAlTh, file = NULL, ROI,intSNP,phasing)
datChr(x, Geno, Annot)
flipGeno(datLeft,datRight, maxSteps)
encodeGeno(BrakePoint, objectGenoDat, BlockSize,lp,db)
encodeGenoAcross(BrakePoint, objectGenoDat, BlockSize,lp,db)
encodeHaplo(BrakePoint, haploDat, BlockSize,lp,db)
getblockLocB12(BrakePoint,BlockSize,lp,db)


}

\description{Internal inveRsion objects.}
\details{These are not to be called by the user.}
\keyword{internal}
