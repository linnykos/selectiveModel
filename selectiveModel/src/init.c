#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _selectiveModel_c_consecutive_true(SEXP);
extern SEXP _selectiveModel_c_construct_midpoints(SEXP);
extern SEXP _selectiveModel_c_gibbs_step(SEXP, SEXP, SEXP, SEXP);
extern SEXP _selectiveModel_c_intersect_intervals(SEXP);
extern SEXP _selectiveModel_c_l2norm(SEXP);
extern SEXP _selectiveModel_c_sample_truncnorm_white(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _selectiveModel_c_theta_in_all_matrix(SEXP, SEXP);
extern SEXP _selectiveModel_c_theta_in_matrix(SEXP, SEXP);
extern SEXP _selectiveModel_c_unique_sort_native(SEXP);
extern SEXP _selectiveModel_c_unlist_native(SEXP);
extern SEXP _selectiveModel_c_which_native(SEXP);
extern SEXP _rcpp_module_boot_module(void);

static const R_CallMethodDef CallEntries[] = {
  {"_selectiveModel_c_consecutive_true",       (DL_FUNC) &_selectiveModel_c_consecutive_true,       1},
  {"_selectiveModel_c_construct_midpoints",    (DL_FUNC) &_selectiveModel_c_construct_midpoints,    1},
  {"_selectiveModel_c_gibbs_step",             (DL_FUNC) &_selectiveModel_c_gibbs_step,             4},
  {"_selectiveModel_c_intersect_intervals",    (DL_FUNC) &_selectiveModel_c_intersect_intervals,    1},
  {"_selectiveModel_c_l2norm",                 (DL_FUNC) &_selectiveModel_c_l2norm,                 1},
  {"_selectiveModel_c_sample_truncnorm_white", (DL_FUNC) &_selectiveModel_c_sample_truncnorm_white, 6},
  {"_selectiveModel_c_theta_in_all_matrix",    (DL_FUNC) &_selectiveModel_c_theta_in_all_matrix,    2},
  {"_selectiveModel_c_theta_in_matrix",        (DL_FUNC) &_selectiveModel_c_theta_in_matrix,        2},
  {"_selectiveModel_c_unique_sort_native",     (DL_FUNC) &_selectiveModel_c_unique_sort_native,     1},
  {"_selectiveModel_c_unlist_native",          (DL_FUNC) &_selectiveModel_c_unlist_native,          1},
  {"_selectiveModel_c_which_native",           (DL_FUNC) &_selectiveModel_c_which_native,           1},
  {"_rcpp_module_boot_module",                 (DL_FUNC) &_rcpp_module_boot_module,                 0},
  {NULL, NULL, 0}
};

extern void R_init_selectiveModel(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
