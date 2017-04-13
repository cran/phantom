#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .Call calls */
  extern SEXP phantom_cummaxC(SEXP);
extern SEXP phantom_descend2C(SEXP, SEXP);
extern SEXP phantom_getNullHetParamC(SEXP, SEXP, SEXP);
extern SEXP phantom_paretoFrontTest2C(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"phantom_cummaxC",           (DL_FUNC) &phantom_cummaxC,           1},
  {"phantom_descend2C",         (DL_FUNC) &phantom_descend2C,         2},
  {"phantom_getNullHetParamC",  (DL_FUNC) &phantom_getNullHetParamC,  3},
  {"phantom_paretoFrontTest2C", (DL_FUNC) &phantom_paretoFrontTest2C, 7},
  {NULL, NULL, 0}
};

void R_init_phantom(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}
