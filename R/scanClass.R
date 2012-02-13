##Define Class##
setClass(Class="scan",
	     representation=representation(leftBP="matrix",rightBP="matrix",leftBP2="matrix",rightBP2="matrix",LogLike="matrix",prob="matrix",ent="matrix",entTh="matrix",bic="matrix",window="numeric"))
           
      
##InveRsion Contructor##
scanInv <-
function(objectHaploCode,window,maxSteps=30,geno=FALSE,saveRes=TRUE,saveBlocks=TRUE)
{
   message("-Scan for inversions-")    

   if(missing(window))
       { 
            warning("Default search window of (0.5Mb) used \n")
            window<-0.5
        }

    if(!is(objectHaploCode,"HaploCode"))
        stop("Error: first argument must be class HaploCode \n") 
  
      b1b2<-objectHaploCode@haploCode
               
      Coorb1b2<-as.numeric(colnames(b1b2))
      numBrakePoints<-NCOL(b1b2)/2
     
      #set up the model  
      seqbrakepoints<-2*(1:numBrakePoints)-1
      leftCoorb1b2<-Coorb1b2[seqbrakepoints]
      rightCoorb1b2<-Coorb1b2[seqbrakepoints+1]
      BlockSize<-objectHaploCode@blockSize
      ls<-attr(objectHaploCode@haploCode,"locilim")

      #build b1b2 blocks, use only left brake point limit to ckeck window search  
      message("-computing inversion model-")
      invResults<-iterateInversionModel(b1b2,window,BlockSize,ls,maxSteps,geno,leftCoorb1b2,rightCoorb1b2)
    
      message("")  
      message("-computation done-")  
      scanRes<- new("scan", 
        leftBP=invResults[[1]],
        rightBP=invResults[[2]],
        leftBP2=invResults[[3]],
        rightBP2=invResults[[4]],
        LogLike=invResults[[5]],
        prob=invResults[[6]],
        ent=invResults[[7]],
        entTh=invResults[[8]],
        bic=invResults[[9]],
        window=window
        )
        
      if(saveRes)
        save(scanRes,file="scanRes.RData")
    
      scanRes 
}

#auxilary functions#
#1. docodes the haplotypes blocks around brake-points into two blocks of size "blockSize" each.
getb1b2 <-
function(BrakePoint,objectHaploCode)
{

  cat(".")
  #set up input data
  HaploDat<-objectHaploCode@haploCode
  BlockSize<-objectHaploCode@blockSize 

  a<-t(sapply(HaploDat[,BrakePoint],int2bit,BlockSize))

  b1<-a[,1:BlockSize] 
  b2<-a[,(BlockSize+1):(2*BlockSize)]

  #define base for coding the haplotypes as decimal integers 
  base2<-sapply((NCOL(b1)-1):0,function(i) 2^(i)) 
  base2Mat<-as.matrix(base2)[,rep(1,NROW(b1))]

  b1code<-as.vector(rowSums(b1*t(base2Mat)))
  b2code<-as.vector(rowSums(b2*t(base2Mat)))

  b12<-cbind(b1code,b2code)
  b12
}

#-within getb1b2 tranformation of integers into binaries -recovers habloptyes.
int2bit <-
function(Num,BlockSize)
{
  
  #set up binary basis
  base2<-sapply((2*BlockSize-1):0,function(i) 2^(i))
  basetemp<-base2
  bin<-rep(0,(2*BlockSize))
  locbin<-1:(2*BlockSize)
  locbintemp<-locbin
  sel<-vector()

  #residual
  res<-Num
  
  while(res>0)
  {
    a<-(basetemp<=res)
    sel<-c(sel,locbintemp[a][1])
    basetemp<-basetemp[a]
	locbintemp<-locbintemp[a]
    res<-res-basetemp[1]
  }

  bin[sel]<-1
  bin
}


#2. call inversion model for all candidate brake-points depending on window type and candidate brake points
iterateInversionModel <-function(b1b2,window,BlockSize,ls,maxSteps,geno,Coor,Coor2,candidatePoints,getInvclass)
{
    if(missing(getInvclass)) 
      getInvclass<-FALSE

    nr<-NROW(b1b2)    
  
    #look for scan of size window (fixed size)
    if(missing(candidatePoints))
    {
      posbrakePoints<-Coor
      int<-length(posbrakePoints)
     
   
      message("") 
      message("  Maximum brake-points to be tested ",int) 
      
      #initialize output variables
      like<-matrix(ncol=1,nrow=int)
      prob<-matrix(ncol=1,nrow=int)
      entInv<-matrix(ncol=1,nrow=int)
      entTh<-matrix(ncol=1,nrow=int)
      bic<-matrix(ncol=1,nrow=int)
      R<-list()

      rightind<-c()
      sel<-(posbrakePoints-posbrakePoints[1])>window
      nr<-NROW(b1b2);


      pnt<-1
      while(sum(sel)!=0 & pnt!=int)
      {
          cat(".")
        rightPnt<-(1:length(posbrakePoints))[sel][1]
        
        if(((posbrakePoints[rightPnt]-posbrakePoints[pnt])<2*window) & (ls[pnt,2]<ls[rightPnt,1]) )
        {            
          dat<-matrix(0,ncol=4,nrow=NROW(b1b2))
          dat[,1:2]<-b1b2[,(2*pnt-1):(2*pnt)]
          dat[,3:4]<-b1b2[,(2*(rightPnt)-1):(2*rightPnt)]
          Nr<-nr+1
          maxindx<-sum(sapply(0:(BlockSize-1), function(x) 2^x))
     
          #finish phasing depending on the strategy in haploCode
          if(geno)
          {
              
              if(attr(b1b2,"phasing")=="forward&inverted" | attr(b1b2,"phasing")=="forward") 
              {
                #datFor
                dat1<-dat[,1]
                dat2<-dat[,2]
                indx12<-flipGeno(dat1,dat2,maxSteps)

                dat3<-dat[,3]
                dat4<-dat[,4]
                indx34<-flipGeno(dat3,dat4,maxSteps)
              
                dat1234<-cbind(dat[,1],dat[indx12,2],dat[,3],dat[indx34,4])
                indx1234<-flipGeno(cbind(dat[,1],dat[indx12,2]),cbind(dat[,3],dat[indx34,4]),maxSteps)

                Datfor<-cbind(dat1234[,1],dat1234[,2],dat1234[indx1234,3],dat1234[indx1234,4])
                Datfor<-rbind(rep(maxindx+1,4),Datfor)

                if(attr(b1b2,"phasing")=="forward")
                { 
                  Datinv<-Datfor[,c(1,3,2,4)]
                }else{ 

                 #datInv
                  dat1<-dat[,1]
                  dat3<-dat[,3]
                  indx13<-flipGeno(dat1,dat3,maxSteps)

                  dat2<-dat[,2]
                  dat4<-dat[,4]
                  indx24<-flipGeno(dat2,dat4,maxSteps)
              
                  dat1324<-cbind(dat[,1],dat[indx13,3],dat[,2],dat[indx24,4])
                  indx1324<-flipGeno(cbind(dat[,1],dat[indx13,3]),cbind(dat[,2],dat[indx24,4]),maxSteps)

                  Datinv<-cbind(dat1324[,1],dat1324[,2],dat1324[indx1324,3],dat1324[indx1324,4])
                  Datinv<-rbind(rep(maxindx+1,4),Datinv)
                }
              }  

              if(attr(b1b2,"phasing")=="inversion&BP")
              {
                #finish phasing
                #phase inverted region
                dat2<-dat[,2]
                dat3<-dat[,3]
                indx23<-flipGeno(dat2,dat3,maxSteps)
                dat23<-cbind(dat[,2],dat[indx23,3])
      
              #inverted region with left BP
                dat1<-dat[,1]
                indx231<-flipGeno(dat23,dat1,maxSteps)               

              #inverted region with right BP
                dat4<-dat[,4]
                indx234<-flipGeno(dat23,dat4,maxSteps)
     
              #all phased data
                Datflipped<-cbind(dat[indx231,1],dat23,dat[indx234,4])
              
              #add fisrt row as reference
                Datfor<-rbind(rep(maxindx+1,4),Datflipped)
                Datinv<-rbind(rep(maxindx+1,4),Datflipped[,c(1,3,2,4)]) 
             }

          }else{
          
              #add fisrt row as reference
              Datfor<-rbind(rep(maxindx+1,4),dat)
              Datinv<-rbind(rep(maxindx+1,4),dat[,c(1,3,2,4)])
          }
                   
          ResInv<-.C("inversionModel",
                    as.double(unlist(Datfor)), 
                    as.double(unlist(Datinv)), 
                    as.integer(maxSteps), 
                    as.integer(Nr),
                    as.double(rep(0,5)), 
                    as.double(rep(0,Nr)), PACKAGE="inveRsion")

                 
          if(getInvclass)
          {          
             if(ResInv[[6]][1]<0.5)
                R[[pnt]]<-(ResInv[[6]])[-1]
              else
                R[[pnt]]<-(1-ResInv[[6]])[-1]
          }       
                                     
          ResInv<-ResInv[[5]]               
          like[pnt,1]<-ResInv[1]
          prob[pnt,1]<-ResInv[2]
          entInv[pnt,1]<-ResInv[3]
          entTh[pnt,1]<-ResInv[4]
          bic[pnt,1]<-ResInv[5]
                   

        }else{

          if(getInvclass)          
            R[[pnt]]<-rep(0,nr)
                
          like[pnt,1]<-NA
          prob[pnt,1]<-NA
          entInv[pnt,1]<-NA
          entTh[pnt,1]<-NA
          bic[pnt,1]<-NA
          

        }  
                
       #control variables (while loop)
       rightind<-c(rightind,rightPnt)
       pnt<-pnt+1
       sel<-posbrakePoints-posbrakePoints[pnt]>window

      }
      
      
    #write output
    pnt<-pnt-1  
    c1<-Coor[1:pnt]
    c2<-Coor[rightind] 
    leftBP<-matrix(c1,ncol=1,nrow=pnt)
    rightBP<-matrix(c2,ncol=1,nrow=pnt)

    c1<-Coor2[1:pnt]
    c2<-Coor2[rightind]
    leftBP2<-matrix(c1,ncol=1,nrow=pnt)
    rightBP2<-matrix(c2,ncol=1,nrow=pnt)
    
    #if candidate brake points are provided  
    }else{
      
      #variables that define search in left and right brakepoint sets    
      posbrakePoints<-Coor
      int<-length(Coor)
      
      lenL<-length(candidatePoints[[1]])
      lenR<-length(candidatePoints[[2]])
      
      intLFT<-sapply(1:lenL, function(x) (1:int)[posbrakePoints==candidatePoints[[1]][x]])
      intRGT<-sapply(1:lenR, function(x) (1:int)[posbrakePoints==candidatePoints[[2]][x]])
      
      #initialize output variables
      like<-matrix(ncol=1,nrow=lenL*lenR)
      prob<-matrix(ncol=1,nrow=lenL*lenR)
      entInv<-matrix(ncol=1,nrow=lenL*lenR)
      entTh<-matrix(ncol=1,nrow=lenL*lenR)
      bic<-matrix(ncol=1,nrow=lenL*lenR)
      R<-list()

      count<-0
      for(lft in intLFT)
      {
          for(rgt in intRGT)
          {
               
              dd<-posbrakePoints[rgt]-posbrakePoints[lft]            
              if(lft<rgt & dd>window &(ls[lft,2]<ls[rgt,1]))
               count<-count+1  
          }
      }


      message("") 
      message("  Number of brake-points to be tested ",count) 

      rightind<-c()
      leftind<-c()

      pnt<-1
      for(lft in intLFT)
      {
          for(rgt in intRGT)
          {
               
              dd<-posbrakePoints[rgt]-posbrakePoints[lft] 
              if( (lft<rgt) & (dd>window) & (ls[lft,2]<ls[rgt,1]) )
              {  
                 cat(".")
                 dat<-matrix(0,ncol=4,nrow=NROW(b1b2))
                 dat[,1:2]<-b1b2[,(2*lft-1):(2*lft)]
                 dat[,3:4]<-b1b2[,(2*rgt-1):(2*rgt)]
                 nr<-NROW(dat)
                 Nr<-nr+1
                 maxindx<-sum(sapply(0:(BlockSize-1), function(x) 2^x))
     
             #finish phasing depending on the strategy in haploCode
          if(geno)
          {
              
              if(attr(b1b2,"phasing")=="forward&inverted" | attr(b1b2,"phasing")=="forward") 
              {
                #datFor
                dat1<-dat[,1]
                dat2<-dat[,2]
                indx12<-flipGeno(dat1,dat2,maxSteps)

                dat3<-dat[,3]
                dat4<-dat[,4]
                indx34<-flipGeno(dat3,dat4,maxSteps)
              
                dat1234<-cbind(dat[,1],dat[indx12,2],dat[,3],dat[indx34,4])
                indx1234<-flipGeno(cbind(dat[,1],dat[indx12,2]),cbind(dat[,3],dat[indx34,4]),maxSteps)

                Datfor<-cbind(dat1234[,1],dat1234[,2],dat1234[indx1234,3],dat1234[indx1234,4])
                Datfor<-rbind(rep(maxindx+1,4),Datfor)

                if(attr(b1b2,"phasing")=="forward")
                { 
                  Datinv<-Datfor[,c(1,3,2,4)]
                }else{ 

                 #datInv
                  dat1<-dat[,1]
                  dat3<-dat[,3]
                  indx13<-flipGeno(dat1,dat3,maxSteps)

                  dat2<-dat[,2]
                  dat4<-dat[,4]
                  indx24<-flipGeno(dat2,dat4,maxSteps)
              
                  dat1324<-cbind(dat[,1],dat[indx13,3],dat[,2],dat[indx24,4])
                  indx1324<-flipGeno(cbind(dat[,1],dat[indx13,3]),cbind(dat[,2],dat[indx24,4]),maxSteps)

                  Datinv<-cbind(dat1324[,1],dat1324[,2],dat1324[indx1324,3],dat1324[indx1324,4])
                  Datinv<-rbind(rep(maxindx+1,4),Datinv)
                }
              }  

              if(attr(b1b2,"phasing")=="inversion&BP")
              {
                #finish phasing
                #phase inverted region
                dat2<-dat[,2]
                dat3<-dat[,3]
                indx23<-flipGeno(dat2,dat3,maxSteps)
                dat23<-cbind(dat[,2],dat[indx23,3])
      
              #inverted region with left BP
                dat1<-dat[,1]
                indx231<-flipGeno(dat23,dat1,maxSteps)               

              #inverted region with right BP
                dat4<-dat[,4]
                indx234<-flipGeno(dat23,dat4,maxSteps)
     
              #all phased data
                Datflipped<-cbind(dat[indx231,1],dat23,dat[indx234,4])
              
              #add fisrt row as reference
                Datfor<-rbind(rep(maxindx+1,4),Datflipped)
                Datinv<-rbind(rep(maxindx+1,4),Datflipped[,c(1,3,2,4)]) 
             }

          }else{
          
              #add fisrt row as reference
              Datfor<-rbind(rep(maxindx+1,4),dat)
              Datinv<-rbind(rep(maxindx+1,4),dat[,c(1,3,2,4)])
          }
               
                    
            ResInv<-.C("inversionModel",
                    as.double(unlist(Datfor)), 
                    as.double(unlist(Datinv)), 
                    as.integer(maxSteps), 
                    as.integer(Nr),
                    as.double(rep(0,5)), 
                    as.double(rep(0,Nr)), PACKAGE="inveRsion")

           
                 if(ResInv[[6]][1]<0.5)
                     R[[pnt]]<-(ResInv[[6]])[-1]
                 else
                     R[[pnt]]<-(1-ResInv[[6]])[-1]
                
                 ResInv<-ResInv[[5]]
                 like[pnt,1]<-ResInv[1]
                 prob[pnt,1]<-ResInv[2]
                 entInv[pnt,1]<-ResInv[3]
                 entTh[pnt,1]<-ResInv[4]
                 bic[pnt,1]<-ResInv[5]
                 
               
                 pnt<-pnt+1
                 
                 rightind<-c(rightind,rgt)
                 leftind<-c(leftind,lft)

              }
          }
    
       }
       
    #write output
    pnt<-pnt-1
    c1<-Coor[leftind]
    c2<-Coor[rightind] 
    leftBP<-matrix(c1,ncol=1,nrow=pnt)
    rightBP<-matrix(c2,ncol=1,nrow=pnt)

    c1<-Coor2[leftind]
    c2<-Coor2[rightind]
    leftBP2<-matrix(c1,ncol=1,nrow=pnt)
    rightBP2<-matrix(c2,ncol=1,nrow=pnt)
 
    }
      
    like<-matrix(like[1:pnt,1],ncol=1,nrow=pnt)
    prob<-matrix(prob[1:pnt,1],ncol=1,nrow=pnt)
    ent<-matrix(entInv[1:pnt,1],ncol=1,nrow=pnt)
    entTh<-matrix(entTh[1:pnt,1],ncol=1,nrow=pnt)
    bic<-matrix(bic[1:pnt,1],ncol=1,nrow=pnt)

             
    if(!missing(candidatePoints) | getInvclass)
    {
      ans<-list(leftBP,rightBP,leftBP2,rightBP2,like=like,prob=prob,ent=ent,ent=entTh,bic=bic,R=R)
    }else{
      ans<-list(leftBP,rightBP,leftBP2,rightBP2,like=like,prob=prob,ent=ent,ent=entTh,bic=bic)
    }   
}

#local haplotyping leaves unmatched the coupling between b1b2 and b3b4 for the 2 chromosomes of a sinlgle subject. 
#within iterateInversionModel, run inversion Model for identyfying most probable match of haplotypes corresponding 
#to left and right brake-points. Fuction call from findInveR setting parameter geno=TRUE
flipGeno<-function(datLeft,datRight,maxSteps)
{
  if(NCOL(datLeft)==2)
  { 
    s1<-sapply(1:NROW(datLeft),function(x) paste(datLeft[x,1],datLeft[x,2],sep=""))
  }else{
    s1<-datLeft
  }

  if(NCOL(datRight)==2)
  { 
    s2<-sapply(1:NROW(datRight),function(x) paste(datRight[x,1],datRight[x,2],sep=""))
  }else{
    s2<-datRight
  }
  
  
  ev<-2*(1:(length(s2)/2))
  odd<-ev-1

  #re-cast the inversion model as a chromosome mathcing task
  datNew<-matrix(ncol=4,nrow=NROW(s1)/2)
  datNew[,1]<-as.numeric(s1[odd])
  datNew[,2]<-as.numeric(s2[odd])
  datNew[,3]<-as.numeric(s2[ev])
  datNew[,4]<-as.numeric(s1[ev])

  nr<-NROW(datNew);
  
  Datfor<-datNew
  Datinv<-datNew[,c(1,3,2,4)]

                    
  inR<-.C("inversionModel",
              as.double(unlist(Datfor)), 
              as.double(unlist(Datinv)), 
              as.integer(maxSteps), 
              as.integer(nr),
              as.double(rep(0,5)), 
              as.double(rep(0,nr)), PACKAGE="inveRsion")

  
 
  R1<-inR[[6]]
   
  ev<-2*(1:nr)
  odd<-ev-1
    
  sel.ev.flip<-ev[R1<0.5]
  sel.odd.flip<-odd[R1<0.5]
    
  flip<-as.vector(rbind(sel.ev.flip,sel.odd.flip))
  no.flip<-as.vector(rbind(sel.odd.flip,sel.ev.flip))   
   
  #ansLeft<-1:nr 
  #ansRight<-1:nr

#  ansRight[no.flip]<-ansRight[flip]
  
  ans<-1:NROW(cbind(datLeft,datRight))
  ans[no.flip]<-ans[flip]
  
  ans 
}    



##plot scan##
setMethod(f="plot",signature=c(x="scan"),
    definition=function(x,y,which="bic",thBic=-Inf,Like,...)
	{
      wplot<-charmatch(which,   c("log", "prob","ent","bic"))
   
      a<-getInv(x,thBic=thBic,Like=Like)
       
      mm<-NROW(a)

     lab<-switch(wplot,"LogLike Ratio",
                  "Probability",
                  "entropy",
                  "BIC Difference")

      plot(c(min(a[1:mm,1]), max(a[1:mm,2])),c(min(a[1:mm,wplot+2]),max(a[1:mm,wplot+2])),ylab=lab, xlab="Segments Tested",pch="",...)
 
      for(ss in 1:mm)
      {
         lines(c(a[ss,1],a[ss,2]),c(a[ss,wplot+2],a[ss,wplot+2]),...)
       }  

      }
	
)


##Show scan#  
setMethod("show","scan",
    function(object)
    {
       cat("-Showing object of class: scan- \n\n")

       if (length(object@leftBP)==0 | length(object@rightBP)==0 | length(object@LogLike)==0 | length(object@prob)==0| length(object@ent)==0 | length(object@entTh)==0 | length(object@bic)==0 | length(object@window)==0)
       { 
             cat("Void: No scan computed \n")
       }else{
       
             cat("Top 10 brake-points with highest Likelihood ratio:\n\n")
             a<-round(as.vector(object@LogLike),2)
             a[is.na(a)]<-0 
             L<-formatC(as.vector(object@leftBP),format="f",digits=,5)
             R<-formatC(as.vector(object@rightBP),format="f",digits=,5)
             L2<-formatC(as.vector(object@leftBP2),format="f",digits=5)
             R2<-formatC(as.vector(object@rightBP2),format="f",digits=5)
             L<-paste(L,L2,sep="-")
             R<-paste(R,R2,sep="-")

             P<-round(as.vector(object@prob),3)
             E<-round(as.vector(object@ent),3)
             ET<-round(as.vector(object@entTh),3)
             B<-round(as.vector(object@bic),3)

             or<-order(-a)[1:10]


             res<-data.frame(LeftBP=L[or],RightBP=R[or],LogLike=a[or],Prob=P[or],BicDiff=B[or])
             row.names(res)<-as.character(1:10) 
              
             print(res)
             cat("\n")
   
             cat("others:\n")
             cat("window: length window (",object@window,") for searching inversion segments \n")
             
       }  
    }
)



##Generic Methods## 
#get inveRsion results into an array$          
setGeneric("getInv",function(object,thBic,rnd,Like){standardGeneric("getInv")})


setMethod("getInv","scan", 
         function(object,thBic,rnd,Like)
         {
            a<-as.vector(object@LogLike)
            a[is.na(a)]<-0
            L<-as.vector(object@leftBP)
            R<-as.vector(object@rightBP)
            P<-as.vector(object@prob)
            E<-as.vector(object@ent)
            ET<-as.vector(object@entTh)
            B<-as.vector(object@bic)
  
            if(missing(rnd))
              rnd<-TRUE
                                
            if(missing(thBic))
              thBic<--Inf
                              
            if(missing(Like))
              Like<-c(0,Inf)
              
              
            sel<-a>Like[1]  & a<Like[2]

            sel<-B>thBic & sel
           
            if(sum(sel)==0)
               stop("Error: no segments selected")
   
                       
            a<-a[sel]
            L<-L[sel]
            R<-R[sel]
            P<-P[sel]
            E<-E[sel]
            ET<-ET[sel]
            B<-B[sel]              

            or<-order(-a)
            if(rnd)
              res<-round(cbind(L[or],R[or],a[or],P[or],E[or],B[or],ET[or]),3)
            else
              res<-cbind(L[or],R[or],a[or],P[or],E[or],B[or],ET[or])
                 
            colnames(res)<-c("  LeftBP","  RightBP","  LogLike", "  Prob", "  Ent","  BIC", "  NumHap")
            row.names(res)<-as.character(1:length(or))
            res
        }
)


#get regions of interest defined by tested segments (of fixed window size) with bic>thBic and overlap among themselves
setGeneric("getROIs",function(object,thBic){standardGeneric("getROIs")})

 
setMethod("getROIs","scan", 
         function(object,thBic)
         {
                  
          if(missing(thBic))
            thBic<-0
           
          a<-getInv(object,thBic=thBic,rnd=FALSE)
  
          #select candidate segments
          or<-order(a[,1])
          s1<-a[or,1]
          s2<-a[or,2]
          d1<-s1[-1]-s1[-length(s1)]
          
          #segments for which the left brake points are no further than a window size appart
          sel<-d1>2*object@window/2

          int<- 1:length(sel)


          LeftBsup<-c(s1[-length(s1)][sel],s1[length(s1)])
          LeftBinf<-c(s1[1],s1[int[sel]+1])

          RightBsup<-c(s2[-length(s1)][sel],s2[length(s2)])
          RightBinf<-c(s2[1],s2[int[sel]+1])

        
          ans<-cbind(LeftBinf,LeftBsup,RightBinf,RightBsup)
          rownames(ans)<-1:NROW(ans)
          message("  ",NROW(ans), " ROI extracted  \n")
          ans

         }
)

      
##listInversion Contructor which is also a method for scan##
setGeneric("listInv",function(object,hapCode,geno,ROI,saveRes,thBic,all,saveROI){standardGeneric("listInv")})

 
setMethod("listInv","scan", 
         function(object,hapCode,geno,ROI,saveRes,thBic,all,saveROI)
         {
          
          if(missing(hapCode)) 
             stop("\n haploCode object from previous results needed\n")
          
          if(missing(geno))
          { 
             message("\n assuming phased data \n")  
             geno<-FALSE
          } 
          

          if(missing(thBic))
            thBic<-0


          if(missing(ROI))
             ROI<-getROIs(object,thBic=thBic)

          if(missing(saveRes))
            saveRes<-TRUE


          if(missing(all))
             all<-TRUE

          if(missing(saveROI))
             saveROI<-TRUE

          #minimum number of candidate brake points (left and right) for not doing an "all" re-run.
          minBP<-5
              
        
          #get previous results  
          leftCoorb1b2<-sort(union(object@leftBP,object@rightBP)) 
          a<-getInv(object,thBic=thBic,rnd=FALSE)
        
    
          results<-list()
            
          if(is.matrix(ROI))
          {
            LCinf<-ROI[,1]
            LCsup<-ROI[,2]
            nROI<-NROW(ROI)
          }else{
            LCinf<-ROI[1]
            LCsup<-ROI[2]
            nROI<-1
          }
          
          numbp<-1 
          invList<-new("inversionList")                 
          for(bp in 1:nROI) 
          {
          
              message("\n")
              message("  doing ROI: ", bp)
              
              #select scanned segments in the ROI
              selbp<-a[,1]>=LCinf[bp] & a[,1]<=LCsup[bp]
             
              P<-a[selbp,4]
              LB<-a[selbp,1]
              RB<-a[selbp,2]
                                                              
              #select models for consistent prediction of subject scan                                                                                                                                                                                                                                                                                                                                                                                     
              #implemented as a list so LPNT and RPNT can take a range of values and all paired combintations can be tested      
                
              if(!all) #fix window
              {
                 ss<-leftCoorb1b2>=min(LB) & leftCoorb1b2<=max(RB)
              }else{
                if(min(RB)<max(LB))
                {        
                   slb<-leftCoorb1b2>=min(LB) & leftCoorb1b2<=min(RB)
                   srb<-leftCoorb1b2>=max(LB) & leftCoorb1b2<=max(RB)
                }else{
                   slb<-leftCoorb1b2>=min(LB) & leftCoorb1b2<=max(LB) & leftCoorb1b2<=min(LB)+object@window/2
                   srb<-leftCoorb1b2>=min(RB) & leftCoorb1b2<=max(RB) & leftCoorb1b2>=max(RB)-object@window/2
              
                }
              
                canLBP<-leftCoorb1b2[slb]
                canRBP<-leftCoorb1b2[srb]

                candidatePoints<-list(canLBP,canRBP)
                 
                ss<-slb | srb
                
           
              }
              
              canp<-leftCoorb1b2[ss]

              selh<-rep(ss,each=2)
                 
              ROIb1b2<-hapCode@haploCode[,selh]              
              
              attr(ROIb1b2,"phasing")<-attr(hapCode@haploCode,"phasing")

              BlockSize<-hapCode@blockSize
              ls<-attr(hapCode@haploCode,"locilim")
              

              #compute all possible combinations in ROI/dist>window 
              if(all)  
                invResults<-iterateInversionModel(ROIb1b2,window=object@window,BlockSize=BlockSize,ls=ls,maxSteps=30,geno=geno,Coor=canp,Coor2=canp,candidatePoints=candidatePoints)
              else
                invResults<-iterateInversionModel(ROIb1b2,window=object@window,BlockSize=BlockSize,ls=ls,maxSteps=30,geno=geno,Coor=canp,Coor2=canp,getInvclass=TRUE)
                            
              RR<-invResults[[10]]
                
              #select only possitive bics
              sb<-invResults[[9]]>thBic
              
              if(length(sb)==0)
              {    
                warning("\n no pair of brakepoints selected for ROI:", bp,"\n")                
                
              }else{

                int<-(1:length(RR))[sb]
                  
                #sort inversion reference with respect to segment of maximum BIC


                r1<-sapply(1:length(RR[[1]]), function(y) mean(sapply(int,function(x) RR[[x]][y]<0.5)))              
                  
                LB<-c(min(LB),max(LB))
                names(LB)<-c("min","max")
              
                RB<-c(min(RB),max(RB))
                names(RB)<-c("min","max")
                             
                results[[numbp]]<- new("inversion", 
              	    classification=r1,
                    leftBP=invResults[[1]][sb],
                    rightBP=invResults[[2]][sb],
                    bic=invResults[[9]][sb],
                    intLeftBP=LB, 
                    intRightBP=RB,
                    invFreq=mean(r1>0.5),
                    RR=lapply(int, function(x) RR[[x]])
                    )
                    
                numbp<-numbp+1

              }
        }
       
       invList<-new("inversionList",results=results)
       
       if(saveRes)
         save(invList,file="invList.RData")                         
          
       invList              
      }
            
)


invROI<-
function(hapCode,geno,ROI,saveRes,thBic)
{
          
          if(missing(hapCode)) 
             stop("\n haploCode object from previous results needed\n")
          
          if(missing(geno))
          { 
             message("\n assuming phased data \n")  
             geno<-FALSE
          } 
          

          if(missing(thBic))
            thBic<--Inf


          if(missing(ROI))
             stop("\n needs ROI \n")
        

          if(missing(saveRes))
            saveRes<-TRUE

        
          #get previous results  
          leftCoorb1b2<-as.numeric(colnames(hapCode@haploCode))[2*((1:(NCOL(hapCode@haploCode)/2)))-1]
                
          results<-list()
            
          #invList<-new("inversionList")                 
              
          slb<-leftCoorb1b2>=ROI[1] & leftCoorb1b2<=ROI[2]
          srb<-leftCoorb1b2>=ROI[3] & leftCoorb1b2<=ROI[4]
              
          canLBP<-leftCoorb1b2[slb]
          canRBP<-leftCoorb1b2[srb]

          candidatePoints<-list(canLBP,canRBP)
                 
          ss<-slb | srb
                        
          canp<-leftCoorb1b2[ss]
          
          selh<-rep(ss,each=2)
                 
          ROIb1b2<-hapCode@haploCode[,selh]              
              
          BlockSize<-hapCode@blockSize
          ls<-attr(hapCode@haploCode,"locilim")

           #add first row as refernce
           maxindx<-sum(sapply(0:(BlockSize-1), function(x) 2^x))
           datNew<-matrix(ncol=NCOL(ROIb1b2),nrow=(NROW(ROIb1b2)+1))
           datNew[1,]<-rep(maxindx+1,NCOL(ROIb1b2))
           datNew[-1,]<-ROIb1b2
           ROIb1b2<-datNew

           attr(ROIb1b2,"phasing")<-attr(hapCode@haploCode,"phasing")

              
          #compute all possible combinations in ROI/dist>window 
          invResults<-iterateInversionModel(ROIb1b2,window=0,BlockSize=BlockSize,ls=ls,maxSteps=30,geno=geno,Coor=canp,Coor2=canp,candidatePoints=candidatePoints)
                            
          RR<-invResults[[10]]
            
                
          #select only possitive bics
          sb<-invResults[[9]]>thBic
                               
          if(length(sb)==0)
          {    
             warning("\n no pair of brakepoints selected for the ROI \n")                
                
          }else{

                int<-(1:length(RR))[sb]

                 for(ii in int)
                   if(RR[[ii]][1]<0.5)
                      RR[[ii]]<-1-RR[[ii]]

                r1<-sapply(1:length(RR[[1]]), function(y) mean(sapply(int,function(x) RR[[x]][y]<0.5)))              
                  
                LB<-c(ROI[1],ROI[2])
                names(LB)<-c("min","max")
              
                RB<-c(ROI[3],ROI[4])
                names(RB)<-c("min","max")
                             
                results[[1]]<- new("inversion", 
              	    classification=r1,
                    leftBP=invResults[[1]][sb],
                    rightBP=invResults[[2]][sb],
                    bic=invResults[[9]][sb],
                    intLeftBP=LB, 
                    intRightBP=RB,
                    invFreq=mean(r1>0.5),
                    RR=lapply(int, function(x) RR[[x]])
                    )
          }
        
       
       invList<-new("inversionList",results=results)
       
       if(saveRes)
         save(invList,file="invList.RData")                         
          
       invList              
}
            

