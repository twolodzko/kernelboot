// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern SEXP kernelboot_cpp_rmvk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kernelboot_cpp_dmvn(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kernelboot_cpp_rmvn(SEXP, SEXP, SEXP, SEXP);
extern SEXP kernelboot_cpp_ruvk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"kernelboot_cpp_rmvk",  (DL_FUNC) &kernelboot_cpp_rmvk,  6},
  {"kernelboot_cpp_dmvn",  (DL_FUNC) &kernelboot_cpp_dmvn,  5},
  {"kernelboot_cpp_rmvn",  (DL_FUNC) &kernelboot_cpp_rmvn,  4},
  {"kernelboot_cpp_ruvk",  (DL_FUNC) &kernelboot_cpp_ruvk,  6},
  {NULL, NULL, 0}
};

void R_init_kernelboot(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
