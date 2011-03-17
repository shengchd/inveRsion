##Define Class##
#GenoDat#
setClass(Class="GenoDat",
	representation=representation(
      genoDat="matrix",
      lociPos="numeric",
      alleleSum="matrix",
      noMissCount="matrix"
	),
	validity=function(object)
	{
	   cat("-Set up object of class: GenoDat- \n")
	   if (NCOL(object@genoDat)!= length(object@lociPos))
	   {
	       stop("number of probes in genoDat and lociPos are not equal \n")
	   }	   	   
	   
	    if (NCOL(object@genoDat)!= NROW(object@alleleSum))
	   {
	       stop("number of probes in genoDat and alleleSum are not equal \n")
	   }	 
	   
	    if (NCOL(object@genoDat)!= NROW(object@noMissCount))
	   {
	       stop("number of probes in genoDat and noMissCount are not equal \n")
	   }	 
	   return(TRUE)
	}
	
)

#GenoDatROI#
setClass(Class="GenoDatROI",
	     representation=representation(ROI="vector"),
         contains="GenoDat",
         validity=function(object)
         {
             cat("-Set up object class: GenoDatROI- \n")
             if(length(object@ROI)!=0) 
               if(object@ROI[1] >= object@ROI[2])
               {
                  print(object@ROI)
                  stop("left ROI limit grater than right ROI limit")
               }
          }
         )
         
setMethod("initialize", signature="GenoDatROI",
           definition=function(.Object,genoDat=matrix(nrow=0,ncol=0),lociPos=numeric(0),
                                alleleSum=matrix(nrow=0,ncol=0),noMissCount=matrix(nrow=0,ncol=0),
                                ROI=vector("numeric"))
                      {
                      
                        if(length(ROI[1])!=0) 
                        {
                          if(length(lociPos)!=0)
                          {
                            if(ROI[1]<=lociPos[1])
                                ROI[1]<-lociPos[1]
                      
                             if(ROI[2]>=lociPos[length(lociPos)])
                               ROI[2]<-lociPos[length(lociPos)]
                          }
                        }  
                          
                        .Object@ROI<-ROI

                        .Object@genoDat<-genoDat
                        .Object@lociPos<-lociPos
                        .Object@alleleSum<-alleleSum
                        .Object@noMissCount<-noMissCount
                        
                        validObject(.Object)
                        return(.Object)   
                    }  
                      
           )
           
                         
##Class constructors##
#GenoDat from file#
setUpGenoDatFile <-
function(file="GenoDat.txt",saveRes=FALSE, sortMinor=TRUE)
{
   
    message("-processing data from file:",file,"-") 
    
    message("  reading file ...")
    genoDat<-read.table(file)
    
    lociPos<-unlist(genoDat[1,])
    genoDat<-genoDat[-1,]

    if(sortMinor)
    {
    
      nsub<-NROW(genoDat)
      caco<-as.factor(rep(1,nsub))
    
      genoDat<-as.matrix(genoDat)
      genoDat<-matrix(as.raw(genoDat+1), nrow=NROW(genoDat), ncol=NCOL(genoDat))
    
      message("  setting: common allele=0, variant allele=1")
      sumAllele<-writeFileGenoDat(genoDat,allesum=TRUE,caco)
    
      genoDat<-read.table("GenoDatTemp.txt")
      unlink("GenoDatTemp.txt")
      genoDat<-t(as.matrix(genoDat))
      
    }else{
    
      genoDat<-as.matrix(genoDat)

    }
    
    message("  setting up GenoDat")
    #get rid of probes with more than 10% with missing values
    nageno<-is.na(genoDat)
    sel<-apply(nageno,2,sum)<NROW(genoDat)*0.1
    
    #set up GenoDatList    
    ans<-new(Class="GenoDat",
         genoDat=genoDat[,sel],
         lociPos=lociPos[sel],
         alleleSum=matrix(unlist(sumAllele[[2]])[sel],ncol=1,nrow=sum(sel)),
         noMissCount=matrix(unlist(sumAllele[[3]])[sel],ncol=1,nrow=sum(sel))
    )
    gDat<-ans
    if(saveRes)
      save(gDat,file="gDat.RData")
    ans

}

#-auxiliary function. write temporary file GenoDatTemp.txt
writeFileGenoDat <-
function(gg,allesum=FALSE,caco)
{
    ggin<-as.vector(gg)
    
    lev<-c(0,1,2,3)
    nprobes<-NCOL(gg)
    nsub<-NROW(gg)
    nlev<-length(lev)

    if(allesum)
    {
        if(missing(caco))
            stop("caco variable has no entries \n")
	
        if(is.null(levels(caco)))
            stop("caco variable must be factor \n")
        
        levcaco<-levels(caco)
        numlevcaco<-length(levcaco)
        levcacoInt<-1:numlevcaco
        cacoInt<-sapply(caco,function(x) levcacoInt[x==levcaco])
    }
    else 
    {
        #irrelevant integers
        cacoInt<-0
        levcaco<-0
        numlevcaco<-0
    }
      
    out<-.C("writeGenoDat",
             as.integer(ggin),
             as.integer(nprobes),
             as.integer(nsub),
             as.integer(lev),
             as.integer(nlev),
             as.integer(allesum),
             as.integer(cacoInt),
             as.integer(levcacoInt),
             as.integer(numlevcaco),
             outInv=as.double(rep(0,nprobes)),
             outAleleSum=as.double(rep(0,numlevcaco*nprobes)),
             outNoMissCount=as.double(rep(0,numlevcaco*nprobes)),PACKAGE="inveRsion");
    ans<-list(out$outInv, t(matrix(out$outAleleSum, nrow=numlevcaco,ncol=nprobes)), t(matrix(out$outNoMissCount,nrow=numlevcaco,ncol=nprobes)));
  
    ans;
}

 

#GenoDat from data matrix of type raw (SNPmatrix)#
setUpGenoDatSNPmat <-
function(Chr,Geno,Annot,saveRes=FALSE,saveGeno=FALSE)
{
    if(dim(Geno)[2]!=dim(Annot)[1])
        stop("Error: number of probes in genotype and annotation data do not coincide. \n")    
	if(missing(Chr))
        stop("Specify Chr (chromosome number). \n")

    message("-processing snpMatrix-")          
    genosel<-datChr(Chr,Geno=Geno,Annot=Annot)

    nsub<-NROW(genosel)
    caco<-as.factor(rep(1,nsub))
    lociPos<-as.vector(Annot[[4]][Annot[[1]]==Chr])
    
    message("  setting: common allele=0, variant allele=1")
    sumAllele<-writeFileGenoDat(genosel,allesum=TRUE,caco)
    
    #write to file to deal with NA
    genoDat<-read.table("GenoDatTemp.txt")
    unlink("GenoDatTemp.txt")
    genoDat<-t(as.matrix(genoDat))
    
    if(saveGeno)
    {
      message("  writing GenoDat.txt")
      write.table(rbind(lociPos,genoDat),file("GenoDat.txt"),col.names=FALSE,row.names=FALSE)
    }
    
    
    message("  setting up GenoDat")
    #get rid of probes with more than 10% with missing values
    nageno<-is.na(genoDat)
    sel<-apply(nageno,2,sum)<NROW(genoDat)*0.1
    
    #set up GenoDatList
    ans<-new(Class="GenoDat",
         genoDat=genoDat[,sel],
         lociPos=lociPos[sel],
         alleleSum=matrix(unlist(sumAllele[[2]])[sel],ncol=1,nrow=sum(sel)),
         noMissCount=matrix(unlist(sumAllele[[3]])[sel],ncol=1,nrow=sum(sel))
    )
    
    gDat<-ans
    if(saveRes)
      save(gDat,file="gDat.RData")
    ans
}

#-auxilary function. gets chromosome data 
datChr<-function(x,Geno,Annot)
{
   redAnnot<-Annot[Annot[[1]]==x,]
   redGeno<-Geno[,Annot[[1]]==x]
   selNames<-order(redAnnot[[1]],redAnnot[[4]])
   ans<-redGeno[,selNames]
}



#class contructor for ROIGenoDat# 
ROIGenoDat <-
function(objectGenoDat,ROI)
{
    if(!is(objectGenoDat,"GenoDat"))
        stop("Error: fisrt argument smust be class GenoDat \n")    
	if(missing(ROI))
       { 
            warning("GenoDatROI object defined with zero ROIs \n")
            ROI=list()
        }
        
    LP<-objectGenoDat@lociPos
    
    if(ROI[1]<=objectGenoDat@lociPos[1])
        ROI[1]<-objectGenoDat@lociPos[1]
                      
    if(ROI[2]>=objectGenoDat@lociPos[length(objectGenoDat@lociPos)])
        ROI[2]<-objectGenoDat@lociPos[length(objectGenoDat@lociPos)]
     
    if(length(ROI)==2)
       sel<-LP>ROI[1] & LP<ROI[2]
  
    if(length(ROI)==4)
       sel<-(LP>ROI[1] & LP<ROI[2]) | (LP>ROI[3] & LP<ROI[4])

    #set up GenoDat
    new(Class="GenoDatROI",
       genoDat=objectGenoDat@genoDat[,sel],
       lociPos=objectGenoDat@lociPos[sel],
       alleleSum=matrix(objectGenoDat@alleleSum[sel,],ncol=NCOL(objectGenoDat@alleleSum),nrow=sum(sel)),
       noMissCount=matrix(objectGenoDat@noMissCount[sel,],ncol=NCOL(objectGenoDat@alleleSum),nrow=sum(sel)),
       ROI=ROI
    )
}


##Plot GenoDat##

setMethod(f="plot",signature=c(x="GenoDat"),
    definition=function(x,y,...)
	{
	       minal<-x@alleleSum/x@noMissCount
		   matplot(x@lociPos/10^6,minal/2,type="p",pch=".",xlab="Probe MB",ylab="Minor Allele Frequency")
	}
	
)

##Show GenoDat##

setMethod("show","GenoDat",
    function(object)
    {
        cat("-Showing object of class: GenoDat- \n")
        
        cat("\n")
        cat("@genoDat: Genotype Data\n")
        
        if(length(object@genoDat)!=0)
        {
             cat("   *", class(object@genoDat[1,1]),  class(object@genoDat),"~", NROW(object@genoDat), "subjects by ", NCOL(object@genoDat), "probes\n")
             if(NCOL(object@genoDat)>=10 & NROW(object@genoDat)>=10)
             {
                cat("   * first 10x10 elements = \n")
                print(formatC(object@genoDat[1:10,1:10]),quote=FALSE)
              }else{
                cat("   * first elements = \n")
                print(formatC(object@genoDat),quote=FALSE)              
              
              }
        }else{print(matrix("numeric",nrow=0,ncol=0))}
        
        cat("\n")
        cat("@lociPos: Probe Corrdinates\n")
        
        if(length(object@lociPos)!=0)
        {
            cat("   *", class(object@lociPos), "~", length(object@lociPos), "probes \n")
            cat("   * first 10 elements = \n")
            ln<-min(length(object@lociPos),10)
            print(formatC(object@lociPos[1:ln]),quote=FALSE)
        }else{print(vector("numeric"))}
        
        cat("\n")
        cat("Other slots... \n")
        
        cat("\n")
        cat("@allaleSum: allele sum at each probe \n")
        cat("@noMissCount: number of no missings values \n")
        
        cat("\n-end showing GenoDat- \n")
    }
)

##Show GenoDatROI##

setMethod("show","GenoDatROI",
    function(object)
    {
        show(as(object,"GenoDat"))
        if(length(object@ROI)!=0)
        {
            cat("\n-The object is defined in the ROI- \n")
            cat(" @ROI\n")
            print(object@ROI)
        }else{
            cat("\n-The object contains no ROI- \n")
        }
            cat("-end showing ROIs- \n")

    }
)