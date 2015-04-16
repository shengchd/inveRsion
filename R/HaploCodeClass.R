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
function(objectGenoDat,blockSize,minAllele,saveRes=TRUE,file=NULL,ROI,intSNP=FALSE, phasing="inversion&BP")
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
          message("  ", paste(ROI[1],ROI[length(ROI)],sep="-"),"\n")
       }

    if(sum(phasing%in%c("inversion&BP","forward", "forward&inverted"))==0)
              stop("Error: phasing parameter should be: inversion&BP,forward or forward&inverted") 

              
    haploCode<-callEncode(objectGenoDat,blockSize,minAllele,file,ROI,intSNP, phasing="inversion&BP")
    
    attr(haploCode,"phasing")<-phasing

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
function(object,BlockSize,MinAlTh,file=NULL,ROI,intSNP, phasing)
{ 
  #select SNPs with minimum allele frequency
  message("  preparing data")
  
  
  if(!missing(object))
  {
    if(!missing(ROI))
    {
       roi<-ROI*10^6
       objectGenoDat<-ROIGenoDat(object,roi)
    }else{
       objectGenoDat<-object
    }
      
    #for snp density
    if(intSNP)
      db<-quantile(diff(unlist(object@lociPos),lag=BlockSize),0.995)/10^6  
    else
      db<-NaN
    
    alleleSum<-objectGenoDat@alleleSum
    noMissCount<-objectGenoDat@noMissCount
  
    selMinorAl<- alleleSum>MinAlTh*noMissCount
    selMinorAl<-as.vector(selMinorAl)
  
    lociPos<-objectGenoDat@lociPos
    lociSel<-(1:length(lociPos))[selMinorAl] 
    lp<-unlist(lociPos)


    #slecte SNPs with a distance af at least blockSize from the borders
    lociNum<-length(objectGenoDat@lociPos)
    selBLoci<-lociSel>=BlockSize & lociSel<=(lociNum-BlockSize)
    lociSel<-lociSel[selBLoci]

    message("  Number of brakepoints: ", length(lociSel))
    #run hapCodeInt for selected SNPs
    if(phasing=="forward")
      ans<-sapply(lociSel,encodeGenoAcross,objectGenoDat,BlockSize,lp,db)


    if(phasing=="inversion&BP" | phasing=="forward&inverted")
      ans<-sapply(lociSel,encodeGeno,objectGenoDat,BlockSize,lp,db)

  }else{
  
    haploDat<-read.table(file=file,header=FALSE)
    subNum<-dim(haploDat)[1]-1
    
    lociPos<-as.vector(haploDat[1,])
    haploDat<-haploDat[-1,]


    #lable frequent allele=0, non frequent=1 
    cat("... labeling frequent allele=0, variant allele =1 \n")
    lablealleles<-function(x)
    {
        ans<-rep(1,NROW(haploDat))
        lab<-sort(-table(unlist(haploDat[,x])))
        nmLB<-names(lab)
        sel1<-haploDat[,x]==as.character(nmLB)[1]
        ans[sel1]<-0
        ans
    }
   
       
    hapLab<-sapply(1:NCOL(haploDat),lablealleles)           
    haploDat<-hapLab

    if(intSNP)
      db<-quantile(diff(unlist(lociPos),lag=BlockSize),0.995)/10^6  
    else
      db<-NaN
    
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
     
    lp<-unlist(lociPos)
     
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
    selBLoci<-lociSel>=BlockSize & lociSel<=(lociNum-BlockSize)
    lociSel<-lociSel[selBLoci]
    
    message("  Number of brakepoints: ", length(lociSel))
    #run hapCodeInt for selected SNPs
    ans<-sapply(lociSel,encodeHaplo,haploDat,BlockSize,lp,db)
       
      
  }

  ans<-matrix(ans,ncol=2*length(lociSel))
 
  lim<-ans[1,]
  ans<-ans[-1,]     
       
  nb<-as.character(lociPos[lociSel]/10^6)
  nc<-as.character(lociPos[lociSel+1]/10^6)
      
  colnames(ans)<-as.vector(rbind(as.character(nb),as.character(nc)))

  sel<-!is.na(ans[2,])
  ans<-ans[,sel] 
  attr(ans,"locilim")<-t(matrix(lim,nrow=2))    
  ans
}

#encode haplotypes into decimal integers (binaires with 2*blocksize digits)-phased data#
encodeHaplo <-
function(BrakePoint,haploDat,BlockSize,lp,db)
{
  cat(".")
    
  if(is.nan(db)) 
      blockLocB12<-(BrakePoint-BlockSize+1):(BrakePoint+BlockSize)
  else 
    blockLocB12<-getblockLocB12(BrakePoint,BlockSize,lp,db)
  
  if(length(blockLocB12)!=BlockSize*2)
  {
    #return NaN
    nn<-as.vector(rep(NaN,NROW(haploDat)))
    ans<-cbind(c(blockLocB12[1],nn),c(-sort(-blockLocB12)[1],nn))
    ans
    
  }else{
  
    blockB1<-haploDat[,blockLocB12[1:BlockSize]]
    blockB2<-haploDat[,blockLocB12[(BlockSize+1):(2*BlockSize)]]
  
    #define binary base for coding the haplotypes as decimal integers
    base2<-sapply((NCOL(blockB1)-1):0,function(i) 2^(i))
    base2Mat<-as.matrix(base2)[,rep(1,NROW(blockB1))]
  
  
    b1<-as.vector(c(rowSums(blockB1*t(base2Mat))))
    b2<-as.vector(c(rowSums(blockB2*t(base2Mat))))

    ans<-cbind(c(blockLocB12[1],b1),c(-sort(-blockLocB12)[1],b2))

    ans

  }  
}

#get local haplotyping of genotypes at each side of the break points, and then code them into decimal integers#
encodeGeno <-
function(BrakePoint,objectGenoDat,BlockSize,lp,db)
{
  cat(".")
  #set up input data
  GenoDat<-objectGenoDat@genoDat
  LociPos<-objectGenoDat@lociPos
  
  if(is.nan(db)) 
    blockLocB12<-(BrakePoint-BlockSize+1):(BrakePoint+BlockSize)
  else 
    blockLocB12<-getblockLocB12(BrakePoint,BlockSize,lp,db)
  
  if(length(blockLocB12)!=BlockSize*2)
  {
  
    #return NaN
    nn<-as.vector(rep(NaN,2*NROW(GenoDat)))
    ans<-cbind(c(blockLocB12[1],nn),c(-sort(-blockLocB12)[1],nn))
    ans
     
  }else{

    #prepare data for haplo.stats, binary data indimessageing precese of minor allele 
    
    blockLocB12Left<-blockLocB12[1:BlockSize]
    blockLocB12Right<-blockLocB12[(BlockSize+1):(2*BlockSize)]

    blockLocList<-list(blockLocB12Left,blockLocB12Right)
    
    res<-c()
    for(side in c(1,2))
    {
    
      blockB12Left<-GenoDat[,blockLocList[[side]]]

      blockB12Labels<-LociPos[blockLocList[[side]]]
      genoDat1<-(blockB12Left!=0)*1
      genoDat2<-(blockB12Left==2)*1

      blockB12bin<-cbind(genoDat1,genoDat2)
      sel<-as.vector(rbind(1:(BlockSize),(BlockSize+1):(2*BlockSize)))
   
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

 
      hapCode1b1<-as.vector(rowSums(hap1[,1:(NCOL(hap1))]*t(base2Mat)))
      hapCode1b1<-hapCode1b1[locMaxProb]

      hapCode2b1<-as.vector(rowSums(hap2[,1:(NCOL(hap2))]*t(base2Mat)))
      hapCode2b1<-hapCode2b1[locMaxProb]

       b1<-as.vector(rbind(hapCode1b1,hapCode2b1))

       res<-cbind(res,b1)
    }

      

    ans<-rbind(c(blockLocB12[1],-sort(-blockLocB12)[1]),res)
    ans
  }
  
}



#get local haplotyping of genotypes across the break points and then oncode them into decimal integers#
encodeGenoAcross <-
function(BrakePoint,objectGenoDat,BlockSize,lp,db)
{

  cat(".")
  #set up input data
  GenoDat<-objectGenoDat@genoDat
  LociPos<-objectGenoDat@lociPos
  
  if(is.nan(db)) 
    blockLocB12<-(BrakePoint-BlockSize+1):(BrakePoint+BlockSize)
  else 
    blockLocB12<-getblockLocB12(BrakePoint,BlockSize,lp,db)
  
  if(length(blockLocB12)!=BlockSize*2)
  {
  
    #return NaN
    nn<-as.vector(rep(NaN,2*NROW(GenoDat)))
    ans<-cbind(c(blockLocB12[1],nn),c(-sort(-blockLocB12)[1],nn))
    ans
     
  }else{

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
    base2<-sapply((NCOL(hap1)/2-1):0,function(i) 2^(i))
    base2Mat<-as.matrix(base2)[,rep(1,NROW(hap1))]

    #select most frequent haplotypes
    subIdFac<-as.factor(hap.f$subj.id)
    locMaxProb<-unlist(tapply(-hap.f$post,subIdFac,order))==1

 
    hapCode1b1<-as.vector(rowSums(hap1[,1:(NCOL(hap1)/2)]*t(base2Mat)))
    hapCode1b1<-hapCode1b1[locMaxProb]

 
    hapCode1b2<-as.vector(rowSums(hap1[,(NCOL(hap1)/2+1):NCOL(hap1)]*t(base2Mat)))
    hapCode1b2<-hapCode1b2[locMaxProb]

    hapCode2b1<-as.vector(rowSums(hap2[,1:(NCOL(hap2)/2)]*t(base2Mat)))
    hapCode2b1<-hapCode2b1[locMaxProb]

 
    hapCode2b2<-as.vector(rowSums(hap2[,(NCOL(hap2)/2+1):NCOL(hap2)]*t(base2Mat)))
    hapCode2b2<-hapCode2b2[locMaxProb]

 
    b1<-as.vector(rbind(hapCode1b1,hapCode2b1))
    b2<-as.vector(rbind(hapCode1b2,hapCode2b2))
    ans<-cbind(c(blockLocB12[1],b1),c(-sort(-blockLocB12)[1],b2))
    ans
  }
  
}



# In EncodeHaplo y EncodeGeno: select SNPs evenly spread aver a fixed block size (given in coordinates) to get a similar SNP density 
#across the crhomosome the block size is such that most of the SNPs (95%) have at least "BlockSize" number of 
#SNPs in it
getblockLocB12<-function(BrakePoint,BlockSize,lp,db)
{
  
  #select rightBlock SNPs within the coordinate block db
  selr<-lp>=lp[BrakePoint]  & lp<= lp[BrakePoint]+db*10^6  
  int<-1:length(lp)
  ln<-length(int[selr])
  if(ln==0)
  {
     selrB<-0
   }else{
   #get SNPs that are close to the limits of "BlockSize" intervals
    part<-seq(lp[selr][1],lp[selr][ln],length.out=BlockSize)

    s<-c(lp[selr][1],lp[selr][ln])
    for (i in (2:(length(part)-1)))
    {
        curr<-lp[selr][-s]
        ans<-curr[order(abs(part[i]-curr))[1]]
        s<-c(s,ans)
     }

     selrB<-int[lp%in%s]
  }

  #select the leftBlocks
  sell<-lp>=lp[BrakePoint]-db*10^6  & lp< lp[BrakePoint]  
  ln<-length(int[sell])
  if(ln==0)
  {
      sellB<-0
   }else{
     #get SNPs that are close to the limits of "BlockSize" intervals
     part<-seq(lp[sell][1],lp[sell][ln],length.out=BlockSize)

     s<-c(lp[sell][1],lp[sell][ln])
    for (i in (2:(length(part)-1)))
    {
       curr<-lp[sell][-s]
       ans<-curr[order(abs(part[i]-curr))[1]]
       s<-c(s,ans)
     }

     sellB<-int[lp%in%s]
   }

  ans<-c(sellB,selrB)
}


##Show HaploCode##
setMethod("show","HaploCode",
    function(object)
    {
    
        cat("-Showing object of class: HaploCode- \n")
        cat("\n")
        cat("@haploCode: Binary code for haplotypes of brake-points flanked by SNP blocks\n")

        if(length(object@haploCode)!=0)
        {
             cat("   *", class(object@haploCode[1,1]),  class(object@haploCode),"~", NROW(object@haploCode), "chromosomes by ", NCOL(object@haploCode)/2, "brake points\n")

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
