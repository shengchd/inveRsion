#ifndef inversionModel_H
#define inversionModel_H


void blockFreq(double *freq, double *block, double *lev, int *nlev, int *nsub, double *out);
void newFreq(double *AA, double *Resp, double *block, double *lev, int *nlev, int *nsub, double *out);
void blockAndLev(double *dat, int *nr, double *block, int *col1, int *col2, int *nlev);
void getFreq(double *block, int *nr, int *nlev, double *outlev, double *outfreq);


void inversionModel(double *dat, int *maxSteps, int *nr, double *outLike, double *outR1);
void writeGenoDat(int *geno, int *numprobes, int *numsub, int *lev, int *numlev, int *sumallele, int *caco, int *levcaco, int *numlevcaco, double *outInv, double *outAleleSum, double *outNoMissCount);


#endif

