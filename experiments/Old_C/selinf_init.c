#include <R.h>
#include <Rinternals.h>

/* . entry points */
extern void sample_truncnorm_white(double *state, double *U, double *directions, double *alphas, double *output, int *pnconstraint, int *pndirection, int *pnstate, int *pburnin, int *pndraw);
static R_NativePrimitiveArgType sample_truncnorm_white_t[] = {
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

static const R_CMethodDef CEntries[] = {
  {"sample_truncnorm_white", (DL_FUNC) &sample_truncnorm_white, 10},
  {NULL, NULL, 0}
};

void R_init_cubature(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
