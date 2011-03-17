##Define Class##           
           
setClass(Class="HaploCode",
	     representation=representation(haploCode="matrix",blockSize="numeric",minAllele="numeric"))
           
setMethod("initialize", signature="HaploCode",
           definition=function(.Object,haploCode=matrix(nrow=0,ncol=0),blockSize=3,minAllele=0.1)
                      {
                        .Object@haploCode<-haploCode
                        .Object@blockSize<-blockSize
                        .Object@minAllele<-minAllele
                        return(.Object)   
                    }  
                      
           )

          
##Class constructor##
#User's call#
codeHaplo <-
function(objectGenoDat,blockSize,minAllele,saveRes=TRUE,file=NULL,ROI)
{
     message("-Code haplotypes-")
    if(missing(objectGenoDat))
    {
       if(length(file)==1)
       {
          message("  reading data from haplotype file:")
          message("  ",file," \n")
       }else{
         stop("Error: no GenoDat found") 
       }     
    }else{
      if(!is(objectGenoDat,"GenoDat"))
          stop("Error: first argument must be class GenoDat \n") 
    }   
    
	if(missing(blockSize))
       { 
            warning("\n Default block size of 3 SNPS used \n")
            blockSize=3
        }

    if(missing(minAllele))
       { 
            warning("\n Default minor allele frequency of 0.1 used \n")
            minAllele=0.1
        }
        
    if(!missing(ROI))
      {
         if(length(ROI)!=2)
           if(length(ROI)!=4)
            stop("Error: ROI must have 2 or 4 components") 

         if(ROI[1]>ROI[2])
            stop("Error: invalid ROI. It must satisfy ROI[1]<ROI[2]") 
         
         if(length(ROI)==4)
           if(ROI[3]>ROI[4] | ROI[2]>ROI[3]) 
              stop("Error: invalid ROI. It must satisfy ROI[1]<ROI[2]<ROI[3]<ROI[4]") 

          message("  coding haplotypes between:")
          message("  ", ROI,"\n")
       }
              
    haploCode<-callEncode(objectGenoDat,blockSize,minAllele,file,ROI)
    
    #remove misssings
 #   sel<-is.na(haploCode)
 #   sel<-apply(sel,2,sum)==0
 #   haploCode<-haploCode[,sel]

    ans<-new(Class="HaploCode",
         haploCode=haploCode,
         blockSize=blockSize,
         minAllele=minAllele
    )
    
    hapCode<-ans
    if(saveRes)
      save(hapCode,file="hapCode.RData")
      
     message("") 
     message("-computation done-")
    ans

}

#-calls contructor using either data from a file or from a genoDat object#
#-Defines ROI for each case; for genoDat does it through ROIGenoDat#
#-calls for genoDat objects perform the local haplotyping (haplo.stats) and
#then encode the haplotypes into decimal integers#  
#-calls for haplotype files then encode the haplotypes into decimal integers#  
 
callEncode <-
function(object,BlockSize,MinAlTh,file=NULL,ROI)
{ 
  #select SNPs with minimum allele frequency
  message("  preparing data")
  
  #build haplotypes with SNP distance jump between them
  jump<-1
  
  if(!missing(object))
  {
    if(!missing(ROI))
    {
       roi<-ROI*10^6
       objectGenoDat<-ROIGenoDat(object,roi)
    }else{
       objectGenoDat<-object
    }
       
    alleleSum<-objectGenoDat@alleleSum
    noMissCount<-objectGenoDat@noMissCount
  
    selMinorAl<- alleleSum>MinAlTh*noMissCount
    selMinorAl<-as.vector(selMinorAl)
  
    lociPos<-objectGenoDat@lociPos
    lociSel<-(1:length(lociPos))[selMinorAl] 

    #slecte SNPs with a distance af at least blockSize from the borders
    lociNum<-length(objectGenoDat@lociPos)
    selBLoci<-lociSel>=BlockSize*jump-jump+1 & lociSel<=(lociNum-BlockSize*jump)
    lociSel<-lociSel[selBLoci]

    message("  Number of brakepoints: ", length(lociSel))
    #run hapCodeInt for selected SNPs
    ans<-sapply(lociSel,encodeGeno,objectGenoDat,BlockSize)

    colnames(ans)<-paste(objectGenoDat@lociPos[lociSel]/10^6,objectGenoDat@lociPos[lociSel+1]/10^6,sep="-")
    
  }else{
  
    haploDat<-read.table(file=file,header=FALSE)
    subNum<-dim(haploDat)[1]-1
    
    lociPos<-as.vector(haploDat[1,])
    haploDat<-haploDat[-1,]
    
    if(!missing(ROI))
    {
        roi<-ROI*10^6
        LP<-lociPos
        if(roi[1]<=lociPos[1])
           roi[1]<-lociPos[1]
                      
        if(roi[2]>=lociPos[length(lociPos)])
           roi[2]<-lociPos[length(lociPos)]
        
        if(length(roi)==2)
          sel<-LP>roi[1] & LP<roi[2]
  
        if(length(roi)==4)
          sel<-(LP>roi[1] & LP<roi[2]) | (LP>roi[3] & LP<roi[4])
                    
        haploDat=haploDat[,sel]
        lociPos=lociPos[sel]

     }
     
    lociNum<-dim(haploDat)[2]
    
    #select of SNPs with more than 10% of allele frequency    
    pair<-2*(1:(NROW(haploDat)/2))
    odd<-pair-1

    genoDat<-haploDat[pair,]+haploDat[odd,]

    alleleSum<-colSums(genoDat)
    noMissCount<-rep(NROW(genoDat),NCOL(genoDat))
  
    selMinorAl<- alleleSum>MinAlTh*noMissCount
    selMinorAl<-as.vector(selMinorAl)
  
    lociSel<-(1:length(lociPos))[selMinorAl] 

    #slecte SNPs with a distance af at least blockSize from the borders
    selBLoci<-lociSel>=BlockSize*jump-jump+1 & lociSel<=(lociNum-BlockSize*jump)
    lociSel<-lociSel[selBLoci]
    
    message("  Number of brakepoints: ", length(lociSel))
    #run hapCodeInt for selected SNPs
    ans<-sapply(lociSel,encodeHaplo,haploDat,BlockSize)

    colnames(ans)<-paste(lociPos[lociSel]/10^6,lociPos[lociSel+1]/10^6,sep="-")
     
  }  
  attr(ans,"lociSel")<-lociSel  
  ans
}

#encode haplotypes into decimal integers (binaires with 2*blocksize digits)-phased data#
encodeHaplo <-
function(BrakePoint,haploDat,BlockSize)
{

  cat(".")
  
  #flanking blocks
  blockLocB12<-(BrakePoint-BlockSize+1):(BrakePoint+BlockSize)
      if(blockLocB12[1]<=0) stop("Brake Point must be a distance BlockSize from borders")
  
  blockB12<-haploDat[,blockLocB12]
  
  #define binary base for coding the haplotypes as decimal integers
  base2<-sapply((NCOL(blockB12)-1):0,function(i) 2^(i))
  base2Mat<-as.matrix(base2)[,rep(1,NROW(blockB12))]
  ans<-as.vector(rowSums(blockB12*t(base2Mat)))
  ans
  
}

#get local haplotyping of genotypes and then oncode them into decimal integers#
encodeGeno <-
function(BrakePoint,objectGenoDat,BlockSize)
{

  cat(".")
  #set up input data
  GenoDat<-objectGenoDat@genoDat
  LociPos<-objectGenoDat@lociPos

  #flanking blocks
  #blockLocB12<-seq(from=(BrakePoint-jump*BlockSize)+jump, to=(BrakePoint+jump*BlockSize),by=jump)
  blockLocB12<-(BrakePoint-BlockSize+1):(BrakePoint+BlockSize)
      if(blockLocB12[1]<=0) stop("Brake Point must be a distance BlockSize from borders")

  #prepare data for haplo.stats, binary data indimessageing precese of minor allele 
  blockB12<-GenoDat[,blockLocB12]
  
  blockB12Labels<-LociPos[blockLocB12]
  genoDat1<-(blockB12!=0)*1
  genoDat2<-(blockB12==2)*1

  blockB12bin<-cbind(genoDat1,genoDat2)
  sel<-as.vector(rbind(1:(2*BlockSize),(2*BlockSize+1):(4*BlockSize)))
  blockB12bin<-data.frame(blockB12bin[,sel])+1

  geno<-setupGeno(blockB12bin, miss.val=c(0,NA),locus.label=blockB12Labels)
  #get haplotypes
  hap.f<-haplo.em(geno, locus.label=blockB12Labels,control=haplo.em.control(insert.batch.size =3*BlockSize,n.try=5))

  #get haplotypes (Binary) for each subject strand
  hap1<-(hap.f$haplotype[hap.f$hap1code,]==2)*1
  hap2<-(hap.f$haplotype[hap.f$hap2code,]==2)*1

  #define binary base for coding the haplotypes as decimal integers
  base2<-sapply((NCOL(hap1)-1):0,function(i) 2^(i))
  base2Mat<-as.matrix(base2)[,rep(1,NROW(hap1))]


  #select most frequent haplotypes
  subIdFac<-as.factor(hap.f$subj.id)
  locMaxProb<-unlist(tapply(-hap.f$post,subIdFac,order))==1

  hapCode1<-as.vector(rowSums(hap1*t(base2Mat)))
  hapCode1<-hapCode1[locMaxProb]

  hapCode2<-as.vector(rowSums(hap2*t(base2Mat)))
  hapCode2<-hapCode2[locMaxProb]
  
  ans<-as.vector(rbind(hapCode1,hapCode2))
  ans
  
}

##Show HaploCode##
setMethod("show","HaploCode",
    function(object)
    {
    
        cat("-Showing object of class: HaploCode- \n")
        cat("\n")
        cat("@HaploCode: Binary code for haplotypes of brake-points flanked by SNP blocks\n")

        if(length(object@haploCode)!=0)
        {
             cat("   *", class(object@haploCode[1,1]),  class(object@haploCode),"~", NROW(object@haploCode), "chromosomes by ", NCOL(object@haploCode), "brake points\n")

             if(NCOL(object@haploCode)>=10 & NROW(object@haploCode)>=10)
             {
                cat("   * first 10x10 elements = \n")
                print(formatC(object@haploCode[1:10,1:10]),quote=FALSE)
              }else{
                cat("   * first elements = \n")
                print(formatC(object@haploCode),quote=FALSE)              
              
              }

        }else{print(matrix("numeric",nrow=0,ncol=0))}           
         
        cat("\n")
        if(length(object@blockSize)!=0)
          cat("@blockSize: Block Size of ", object@blockSize, "SNPs \n")
        else
          cat("@blockSize: -NULL-\n")
        
        if(length(object@minAllele)!=0)
          cat("@minAllele: Brake points with minimum allele frequency of ", object@minAllele, "\n \n")
        else
          cat("@minAllele: -NULL-\n\n")
          
        cat("-end showing HaploCode- \n")

    }
)
