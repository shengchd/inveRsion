##Define inversion##
setClass(Class="inversion",
	     representation=representation(classification="vector",leftBP="vector",rightBP="vector",bic="vector",intLeftBP="vector", intRightBP="vector",invFreq="numeric",RR="list"))

##Define inversionList##
setClass(Class="inversionList",
	     representation=representation(results="list"))


##Show inversionList#  
setMethod("show","inversionList",
    function(object)
    {
       cat("-Showing object of class: inversionList- \n\n")

       mat<-c()
       if (length(object@results)==0)
       { 
             cat("Void: No inversions computed \n")
       }else{
         for(nroi in 1:length(object@results))
         {
             if(class(object@results[[nroi]])!="inversion")
             {
                cat("List entry is not of class inversion")
                mm<-rep(NA,6)
             }else{
                mm<-c(object@results[[nroi]]@intLeftBP,object@results[[nroi]]@intRightBP,max(object@results[[nroi]]@bic),object@results[[nroi]]@invFreq,length(object@results[[nroi]]@bic))
                names(mm)<-c()
             }
              
             mat<-rbind(mat,mm)  
         }    
         
         rownames(mat)<-1:length(object@results)
         colnames(mat)<-c("LBPmin","LBPmax","RBPmin","RBPmax","MaxBic","invFreq","ModelNum")
         print(mat)   
       }
   }
)   

##plot inversionList##
setMethod(f="plot",signature=c(x="inversionList"),
    definition=function(x,y,wROI,...)
	{
       if(missing(wROI))
       {
          a<-c()
          for(nroi in 1:length(x@results))
          {
            mm<-cbind(x@results[[nroi]]@leftBP,x@results[[nroi]]@rightBP,x@results[[nroi]]@bic)
            names(mm)<-c()
            a<-rbind(a,mm)  
          }    

         lab<-"BIC"
         nr<-NROW(a)
         plot(c(min(a[1:nr,1]), max(a[1:nr,2])),c(min(a[1:nr,3]),max(a[1:nr,3])),ylab=lab, xlab="Segments Tested",pch="",...)
 
         for(ss in 1:nr)
         {
            lines(c(a[ss,1],a[ss,2]),c(a[ss,3],a[ss,3]))
          }  
        }else{
        
         hist(x@results[[wROI]]@classification,xlab="average responsibilities of subjects across models",main=c("Histogram of subject responsibilities for ROI: ",
         paste(round(x@results[[wROI]]@intLeftBP[1],3),round(x@results[[wROI]]@intRightBP[2],3),sep="-")),...)
      
        }
    
    }
)


#generic functions
setGeneric("getClassif",function(object,thBic,wROI,bin,geno,id){standardGeneric("getClassif")})


setMethod("getClassif","inversionList",         
         function(object,thBic,wROI,bin,geno,id)
         {

            if(missing(wROI))
             wROI<-1
            
            if(missing(thBic))
             thBic<-0

            if(missing(bin))
             bin<-TRUE

            if(missing(geno))
             geno<-TRUE

            if(missing(id))
             id<-paste("sub",1:length(object@results[[wROI]]@RR[[1]]),sep="")
         
            RR<-object@results[[wROI]]@RR

            sb<-object@results[[wROI]]@bic>thBic
            ii<-(1:length(RR))[sb]

            #for each subject see the average classification given across models

            r<-sapply(1:length(RR[[1]]), function(y) mean(sapply(ii,function(x) RR[[x]][y]<0.5)))
           
            message("\n getting mean classification across ", length(ii), " models \n")

           
            if(geno)
            {
            
              even<-2*(1:(length(r)/2))
              odd<-even-1
              homInv<-r[even]*r[odd]
              het<-(1-r[even])*r[odd]+r[even]*(1-r[odd])
              hom<-(1-r[even])*(1-r[odd])
              id<-id
              ans<-data.frame(id=id, hom=hom, het=het, homInv=homInv) 
              ans
            
            }else{
            
              if(bin)
              { 
                ans<-as.numeric(r>0.5)
                names<-rep(id,each=2)
                ans
              }else{
                ans<-r
                names<-rep(id,each=2)
                ans
              }
            }
         }

)                


##Define Class##
setClass(Class="accuracy", representation=representation(out="matrix"))


##Class Contructor which is also a method for inversionList##
setGeneric("accBic",function(object,mem,classFile,nsub,npoints,geno,wROI){standardGeneric("accBic")})

 
setMethod("accBic","inversionList", 
         function(object,mem,classFile,nsub,npoints,geno,wROI)
         {

            if(missing(mem))
            { 
              if(missing(classFile)) 
                stop("\n please provide either vector or file with the inversion staus for the chromosomes of all subjects")

              mem<-read.table(file=classFile,header=FALSE)
              mem<-unlist(mem)
            }

            if(missing(nsub))
              stop("\n plase provide number of subjects in nsub")           
              
            if(missing(geno))
              geno<-FALSE  

            if(missing(wROI))
              wROI<-1
              
            
            nchr<-1:(2*nsub)
            ev<-2*(1:nsub)
            odd<-ev-1
            
            #genotype classification
            popGeno<-odd%in%nchr[mem]+ev%in%nchr[mem]

            #retrive classification for each probe sequence (fixed window model) from object
            RR<-object@results[[wROI]]@RR
            
            ac<-c()
            prob<-c()
 
            bicInt<- seq(0,max(object@results[[1]]@bic),length.out=npoints)
            bicInt<-bicInt[-length(bicInt)]
             message("\n computing accuracy for ", npoints," bic thresholds: ")

            for(thBic in bicInt)
            {
                #select models with bic>thBic
                cat(".")
                sb<-object@results[[wROI]]@bic>thBic
                ii<-(1:length(RR))[sb]
                
                #for each subject see the average classification given across models 
                r<-sapply(1:length(RR[[1]]), function(y) mean(sapply(ii,function(x) RR[[x]][y]<0.5)))
                
                if(geno)
                {
                  #binirize average classification into a final classification for each subject 
                  PHits<-odd%in%nchr[r>0.5]+ev%in%nchr[r>0.5] 
                
                  #compare with real classification
                  ac<-c(ac, mean(popGeno==PHits) )
                }else{
                  PHits<-nchr%in%nchr[mem]
                  popHaplo<-nchr%in%nchr[r>0.5]
                  ac<-c(ac,mean(popHaplo==PHits))
                
                }   
                  
                prob<-c(prob,mean(r>0.5))

            }
            
            
           ans<-cbind(bicInt,prob,ac) 
    
           new("accuracy", out=ans)   


    }
)


##plot accuracy##
setMethod(f="plot",signature=c(x="accuracy"),
    definition=function(x,y,w="a",...)
	{
        wplot<-charmatch(w,c("f", "a"))
 
        lab<-switch(wplot,"Inversion Frequency","accuracy")

        plot(x@out[,1],x@out[,1+wplot],t="l",xlab="BIC",ylab=lab,...)
        points(x@out[,1],x@out[,1+wplot])

    }
)

         
