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
\alias{encodeHaplo}
\alias{manual}



\usage{
getb1b2(BrakePoint, objectHaploCode)
int2bit(Num, BlockSize)
iterateInversionModel(b1b2,window,BlockSize,ls,maxSteps,geno,Coor,Coor2,candidatePoints,getInvclass) 
writeFileGenoDat(gg, allesum = FALSE, caco)
callEncode(object, BlockSize, MinAlTh, file = NULL, ROI)
datChr(x, Geno, Annot)
flipGeno(datOrg, maxSteps)
encodeGeno(BrakePoint, objectGenoDat, BlockSize)
encodeHaplo(BrakePoint, haploDat, BlockSize)

}

\description{Internal inveRsion objects.}
\details{These are not to be called by the user.}
\keyword{internal}