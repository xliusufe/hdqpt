#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _QPT_Test(void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_QPT_Test", (DL_FUNC) &_QPT_Test, 4},
    {NULL, NULL, 0}
};

void R_init_hdqpt(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}