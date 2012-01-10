/*
Function to test linking Lapack library with R.
1. create a file called Makevars in the directory and write in it the single line 
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

2.compile with R CMD SHLIB inversionModel.c

3. open R and do 


check 1:
dyn.load("inversionModel.so")
load(file="testblockFreq.Rdata")
out<-.C("blockFreq", as.double(freq), as.double(block), as.integer(lev), as.integer(nlev),as.integer(nsub),as.double(rep(0,nsub)))


check 2:
dyn.load("inversionModel.so")
load(file="testNewFreq.Rdata")
out<-.C("newFreq", as.double(Resp), as.double(block), as.double(lev), as.integer(nlev),as.integer(nsub), as.double(rep(0,nlev)))

check 3:
dyn.load("inversionModel.so")
load(file="testInversionModel.Rdata")
nlev<-1
out<-.C("blockAndLev", as.double(unlist(dat)), as.integer(nr),as.double(rep(0,nr)), as.integer(2), as.integer(4), as.integer(nlev))


out2<-.C("getFreq", as.double(out[[3]]), as.integer(nr), as.integer(out[[6]]), as.double(rep(0,out[[6]])),  as.double(rep(0,out[[6]])))

CALL:

dyn.load("inversionModel.so")
load(file="testInversionModel.Rdata")
out<-.C("inversionModel", as.double(unlist(dat)), as.integer(maxSteps), as.integer(nr),as.double(rep(0,2)), as.double(rep(0,nr)))


*/

#include "inversion.h"
#include <stdio.h>
#include <R.h>
#include <R_ext/Lapack.h> 
#include <R_ext/BLAS.h>


void blockFreq(double *freq, double *block, double *lev, int *nlev, int *nsub, double *out)
{
  int level, subject;
  
  for (subject=0; subject<*nsub; ++subject) 
  {
    for (level=0; level<*nlev; ++level) 
    {
      if(block[subject]==lev[level])
      {
        out[subject]=freq[level];
      }
     }
   }
  return;
}


/*********************/


void newFreq(double *AA, double *Resp, double *block, double *lev, int *nlev, int *nsub, double *out)
{

  int level, subject, row, col, inx;  
  double freqLev[*nlev];
  double bb[*nlev];
 
  double maxfreq;
  
  /* compute freqLev */
  for (level=0; level<*nlev; ++level) 
    freqLev[level]=0.0;

      
  for (subject=0; subject<*nsub; ++subject) 
  {
    for (level=0; level<*nlev; ++level) 
    {
      if(block[subject]==lev[level])
      {
       freqLev[level]=freqLev[level]+Resp[subject];       
      }
     }
   }

  /*register max freq*/
  maxfreq=-1.0;
  inx=0; 
  for (level=0; level<*nlev; ++level) 
  {
    if(freqLev[level]>maxfreq)
    {
      maxfreq=freqLev[level];
      inx=level;  
    }
   }
      
   /*initialize AA , we costruc the transpose of the matrix since the linear algebra routine takes the transpose*/
   
    for (row=0; row<*nlev; ++row) 
    {
        for (col=0; col<*nlev; ++col) 
        {
          if (row==col)
          {
            AA[row*(*nlev)+col]=1.0;
          }
          else
          {
            AA[row*(*nlev)+col]=0.0;
          }
          
          if (row==inx)
          {
            AA[row*(*nlev)+col]=-freqLev[col]/maxfreq;
          }
          
          if (col==inx)
          {
            AA[row*(*nlev)+col]=1;
          }
        }  
    
    }
        
    /*define bb*/
    for (row=0; row<*nlev; ++row) 
    {
      bb[row]=0.0;
      if(row==inx)
        bb[inx]=1;
    }
    
    
    /*solve linear algebra*/
    int nr=1;
    int *pnr;
    pnr=&nr;
    *pnr=nr;
    
    int inf=0;
    int *pinf;
    pinf=&inf;
    *pinf=inf;
    
    int ipiv[*nlev];
   
    F77_NAME(dgesv)(nlev, pnr, AA, nlev, ipiv, bb, nlev, pinf);
                
    for (row=0; row<*nlev; ++row) 
     out[row]=bb[row];   
}

/**************/
/*reduces data to haplotype frequencies for block*/ 

void blockAndLev(double *dat, int *nr, double *block, int *col1, int *col2, int *nlev)
{

   int ncol, row, level, found, c1, c2;
   int levelleft[*nr], levelright[*nr];
   
   *nlev=1;
   
   c1=*col1-1;
   c2=*col2-1;
   
      
   levelleft[*nlev-1]=dat[(*nr)*(c1)+0];
   levelright[*nlev-1]=dat[(*nr)*(c2)+0];
   
   for (row=0; row<*nr; ++row) 
   {

      found=0;
      for(level=0; level<*nlev; ++level)
      {
        if(dat[(*nr)*(c1)+row]==levelleft[level] & dat[(*nr)*(c2)+row]==levelright[level])
        {
          block[row]=level+1;
          found=1;
        }
      }
      
      if(found==0)
      {  
        *nlev=*nlev+1;
        levelleft[*nlev-1]=dat[(*nr)*(c1)+row];
        levelright[*nlev-1]=dat[(*nr)*(c2)+row];
        block[row]=*nlev;
      }
       
   }

}


/*************/
/*get frequencies*/ 

void getFreq(double *block, int *nr, int *nlev, double *outlev, double *outfreq)
{
  int level, row;
  
  for (level=0; level<*nlev; ++level) 
      outlev[level]=level+1;

   /*compute frequency of haplotypes for block b1b2*/
   for(level=0; level<*nlev; ++level)
      outfreq[level]=0; 
         
   for(row=0; row<*nr; ++row)
   {
      for(level=0; level<*nlev; ++level)
      {
        if(block[row]==outlev[level])
           outfreq[level]=outfreq[level]+1; 
       }
    }
     
   for(level=0; level<*nlev; ++level)
     outfreq[level]=outfreq[level]/(*nr); 
}


/****/
void inversionModel(double *datforward, double *datinverted, int *maxSteps, int *nr, double *outLike, double *outR1)
{

   int ncol, row, level, i;
   int *col1, *col2, *nlevb1b2, *nlevb3b4, *nlevb1b3, *nlevb2b4, *nlevb1b2b3b4;
/*   double BB[(*nr)*(*nr)];*/
  double *BB= (double *)malloc(((*nr)*(*nr))*sizeof(*BB));

   /* compute block  b1b2, level its haplotypes and get their frequencies*/
   int c1=1, c2=2, n1=1;   
   nlevb1b2=&n1; col1=&c1; col2=&c2;
   *nlevb1b2=n1; *col1=c1; *col2=c2;

   double blockb1b2[*nr];

   blockAndLev(datforward, nr, blockb1b2, col1, col2, nlevb1b2);
   
   double levb1b2[*nlevb1b2];
   double nor0[*nlevb1b2];
   getFreq(blockb1b2, nr, nlevb1b2, levb1b2, nor0);
   /* */

   /* compute block  b3b4, level its haplotypes and get their frequencies*/
   c1=3; c2=4; 
   int n2=1;
   nlevb3b4=&n2; col1=&c1; col2=&c2;
   *nlevb3b4=n2; *col1=c1; *col2=c2;

   double blockb3b4[*nr];

   blockAndLev(datforward, nr, blockb3b4, col1, col2, nlevb3b4);
   
   double levb3b4[*nlevb3b4];
   double nor1[*nlevb3b4];
   getFreq(blockb3b4, nr, nlevb3b4, levb3b4, nor1);
   /* */
 
   /* compute block  b1b3, level its haplotypes and get their frequencies, or b1b2 of the datainverted*/
   c1=1; c2=2;   
   int n3=1;
   nlevb1b3=&n3; col1=&c1; col2=&c2;
   *nlevb1b3=n3; *col1=c1; *col2=c2;

   double blockb1b3[*nr];

   blockAndLev(datinverted, nr, blockb1b3, col1, col2, nlevb1b3);
   
   double levb1b3[*nlevb1b3];
   double inv0[*nlevb1b3];
   getFreq(blockb1b3, nr, nlevb1b3, levb1b3, inv0);
   /* */
   
   /* compute block  b2b4, level its haplotypes and get their frequencies, or b3b4 of the datainverted*/
   c1=3; c2=4;
   int n4=1;
   nlevb2b4=&n4; col1=&c1; col2=&c2;
   *nlevb2b4=n4; *col1=c1; *col2=c2;

   double blockb2b4[*nr];

   blockAndLev(datinverted, nr, blockb2b4, col1, col2, nlevb2b4);
   
   double levb2b4[*nlevb2b4];
   double inv1[*nlevb2b4];
   getFreq(blockb2b4, nr, nlevb2b4, levb2b4, inv1);
   /* */

   /* compute the log-likelihood of the non-inverted model*/
   double r01[*nr], r11[*nr];
   double LoglikeNor, BicNor; 
  
   blockFreq(nor0, blockb1b2, levb1b2, nlevb1b2, nr, r01);
   blockFreq(nor1, blockb3b4, levb3b4, nlevb3b4, nr, r11);  
   
   LoglikeNor=0.0;   
   for (row=0; row<*nr; ++row) 
   LoglikeNor=log((r01[row])*(r11[row]))+LoglikeNor;
      
   BicNor=-2*LoglikeNor+(*nlevb1b2+*nlevb3b4-2)*log(*nr); 

/*BicNor=-2*LoglikeNor+(*nlevb1b2+*nlevb3b4-2)*2; 
   
   /* compute block  b1b2b3b4, level its haplotypes and get their frequencies*/
   c1=1; c2=2;
   int n5=1;
   nlevb1b2b3b4=&n5; col1=&c1; col2=&c2;
   *nlevb1b2b3b4=n5; *col1=c1; *col2=c2;

   double blockb1b2b3b4[*nr];
   double datInv[2*(*nr)];   
   
    for(row=0; row<*nr; ++row)
   {
       datInv[row]=blockb1b2[row];
       datInv[row+(*nr)]=blockb3b4[row];
   }       

                    
   blockAndLev(datInv, nr, blockb1b2b3b4, col1, col2, nlevb1b2b3b4);
   
   double levb1b2b3b4[*nlevb1b2b3b4];
   double nor12[*nlevb1b2b3b4];
   getFreq(blockb1b2b3b4, nr, nlevb1b2b3b4, levb1b2b3b4, nor12);
   /* */

   /*compute the entropy of blocks b1b2b3b4*/
   
   double entNor12=0.0;
   for (level=0; level<*nlevb1b2b3b4; ++level)
      entNor12=log(nor12[level])*(nor12[level])+entNor12;
   
   outLike[2]= -entNor12;
  
    
   /*responsabilities*/
   double r1[*nr], r2[*nr], R1[*nr], R2[*nr];
   double nor0Old[*nlevb1b2], nor1Old[*nlevb3b4], inv0Old[*nlevb1b3], inv1Old[*nlevb2b4], prob0Old;
   
   double prob0=0.95;
   
   
   /* EM loop*/
   double tol=1.0, mintol=0.0000001;
   int steps=1;
   
   
   while(steps<*maxSteps & tol>mintol)
   {
     /*store old variables*/ 
      prob0Old=prob0;
      
      for (level=0; level<*nlevb1b2; ++level)
        nor0Old[level]=nor0[level];

      for (level=0; level<*nlevb3b4; ++level)
        nor1Old[level]=nor1[level];
        
      for (level=0; level<*nlevb1b3; ++level)
        inv0Old[level]=inv0[level];
        
      for (level=0; level<*nlevb2b4; ++level)
        inv1Old[level]=inv1[level];
      /* */   


      /*compute responsabilities*/
      blockFreq(nor0, blockb1b2, levb1b2, nlevb1b2, nr, r01);
      blockFreq(nor1, blockb3b4, levb3b4, nlevb3b4, nr, r11);  

      
      for (row=0; row<*nr; ++row) 
        r1[row]=prob0*(r01[row])*(r11[row]);

      blockFreq(inv0, blockb1b3, levb1b3, nlevb1b3, nr, r01);
      blockFreq(inv1, blockb2b4, levb2b4, nlevb2b4, nr, r11);  

      for (row=0; row<*nr; ++row) 
        r2[row]=(1-prob0)*(r01[row])*(r11[row]);
      
      for (row=0; row<*nr; ++row) 
      {
        R1[row]=r1[row]/(r1[row]+r2[row]);
        R2[row]=r2[row]/(r1[row]+r2[row]);
      }
      /* */


     /*new probability of no-inversion*/
     prob0=0.0;
     for (row=0; row<*nr; ++row) 
     {
       prob0 = R1[row]+prob0;
     }
     
     prob0=prob0/(*nr);

     /*new haplotype frequencies*/
     newFreq(BB,R1,blockb1b2,levb1b2,nlevb1b2,nr,nor0);
     newFreq(BB,R1,blockb3b4,levb3b4,nlevb3b4,nr,nor1);
     newFreq(BB,R2,blockb1b3,levb1b3,nlevb1b3,nr,inv0);
     newFreq(BB,R2,blockb2b4,levb2b4,nlevb2b4,nr,inv1);

     /*compute tolerance*/
     tol=(prob0-prob0Old)*(prob0-prob0Old);
     
     for (level=0; level<*nlevb1b2; ++level)    
       tol=tol+(nor0Old[level]-nor0[level])*(nor0Old[level]-nor0[level]);
       
     for (level=0; level<*nlevb3b4; ++level)
       tol=tol+(nor1Old[level]-nor1[level])*(nor1Old[level]-nor1[level]);
        
     for (level=0; level<*nlevb1b3; ++level)
       tol=tol+(inv0Old[level]-inv0[level])*(inv0Old[level]-inv0[level]);
        
     for (level=0; level<*levb2b4; ++level)
       tol=tol+(inv1Old[level]-inv1[level])*(inv1Old[level]-inv1[level]);


     for (level=0; level<*nlevb2b4; ++level)    
        outR1[level]=inv1[level];     
   
     tol=sqrt(tol);
     steps=steps+1;  
     
   }

   /*get last values to compute likelihood of the complete inversion model*/
   /*compute responsabilities*/
    blockFreq(nor0, blockb1b2, levb1b2, nlevb1b2, nr, r01);
    blockFreq(nor1, blockb3b4, levb3b4, nlevb3b4, nr, r11);  

      
    for (row=0; row<*nr; ++row) 
      r1[row]=prob0*(r01[row])*(r11[row]);

    blockFreq(inv0, blockb1b3, levb1b3, nlevb1b3, nr, r01);
    blockFreq(inv1, blockb2b4, levb2b4, nlevb2b4, nr, r11);  

    for (row=0; row<*nr; ++row) 
      r2[row]=(1-prob0)*(r01[row])*(r11[row]);

   double LoglikeInv,BicInv; 
   for (row=0; row<*nr; ++row) 
      LoglikeInv=log(r1[row]+r2[row])+LoglikeInv;
      
   BicInv=-2*LoglikeInv+(*nlevb1b2+*nlevb3b4+*nlevb1b3+*nlevb2b4+1-4)*log(*nr); 

/*BicInv=-2*LoglikeInv+(*nlevb1b2+*nlevb3b4+*nlevb1b3+*nlevb2b4+1-4)*2;*/
   
   
   /*log-likelihood*/
   outLike[0]=2.0*(LoglikeInv-LoglikeNor);
       
   /*BIC*/
   outLike[4]=(BicNor-BicInv);
   
   /*final probability*/
   outLike[1]=prob0;
   
    /*output number of forward haplotypes*/  
   outLike[3]=*nlevb1b3+*nlevb2b4-1;
   
   /*responsbility of no inversion*/
   for (row=0; row<*nr; ++row)    
   outR1[row]=R1[row];
   
   free(BB);
   BB=NULL; 
}



/*
dyn.load("inversionModel.so")

load(file="testInversionModel.Rdata")
out<-.C("inversionModel", as.double(unlist(dat)), as.integer(maxSteps), as.integer(nr),as.double(rep(0,2)), as.double(rep(0,nr)))
*/


