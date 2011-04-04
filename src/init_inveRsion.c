#include <R.h>
#include <Rinternals.h>

#include "inversion.h"
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>


static R_NativePrimitiveArgType inversionModel_in[5] = {REALSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType writeGenoDat_in[12] = {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};

R_CMethodDef CEntries[]  = {
        {"inversionModel", (DL_FUNC) &inversionModel, 5, inversionModel_in},        
        {"writeGenoDat", (DL_FUNC) &writeGenoDat, 12, writeGenoDat_in}, 
        {NULL, NULL, 0}
};


void R_init_inveRsion(DllInfo *ddl){

  R_registerRoutines(ddl, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(ddl, FALSE);

  R_RegisterCCallable("inveRsion", "inversionModel", (DL_FUNC)inversionModel);
  R_RegisterCCallable("inveRsion", "writeGenoDat", (DL_FUNC)writeGenoDat);
}