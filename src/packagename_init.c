#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void fnlnf_C(void *, void *, void *, void *, void *, void *, void *);
extern void fnMH_C(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"fnlnf_C", (DL_FUNC) &fnlnf_C, 7},
  {"fnMH_C",  (DL_FUNC) &fnMH_C,  7},
  {NULL, NULL, 0}
};

void R_init_AdMit(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}