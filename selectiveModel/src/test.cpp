#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector test() {
  Function gc("gc");

  SEXP character = Rf_mkChar("ouch");
  String string = String("ouch");

  // force big allocation + gc
  Shield<SEXP> other(Rf_allocVector(INTSXP, 1E6));
  gc();

  Rprintf("CHARSXP:\n");
  Rf_PrintValue(character);

  Rprintf("STRSXP:\n");
  Rf_PrintValue(wrap(string));

  return wrap(other);
}
