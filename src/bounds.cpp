#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix bounds_cpp(NumericVector x, NumericVector t, NumericVector n) {
  // Create Basic Pieces
  NumericVector omx = 1 - x;
  NumericVector Nb = x * n;
  NumericVector Nw = omx * n;
  int p = x.length();

  // Basic calculations
  NumericVector tx = t / x;
  NumericVector tomx = t / omx;
  NumericVector tomxx = tx - (omx / x);
  NumericVector txx = tomx - (x / (1 - x));

  // Create Output Pieces
  NumericVector LbetaB(p); // lower bound B
  NumericVector UbetaB(p); // upper bound B
  NumericVector LbetaW(p); // lower bound W
  NumericVector UbetaW(p); // upper bound W

  // Fill output vectors
  for (int i = 0; i < p;  i++) {
    if (x(i) == 0) {
      LbetaB(i) = NA_REAL;
      UbetaB(i) = NA_REAL;
      LbetaW(i) = t(i);
      UbetaW(i) = t(i);
    } else if (x(i) == 1) {
      LbetaB(i) = t(i);
      UbetaB(i) = t(i);
      LbetaW(i) = NA_REAL;
      UbetaW(i) = NA_REAL;
    } else {
      LbetaB(i) = std::max(0.0, (double) tomxx(i));
      UbetaB(i) = std::min(1.0, (double) tx(i));
      LbetaW(i) = std::max(0.0, (double) txx(i));
      UbetaW(i) = std::min(1.0, (double) tomx(i));
    }
  }

  // Assemble output pieces
  NumericMatrix out = cbind(LbetaB, UbetaB, LbetaW, UbetaW);
  colnames(out) = CharacterVector::create("LbetaB", "UbetaB", "LbetaW", "UbetaW");
  return out;
}
