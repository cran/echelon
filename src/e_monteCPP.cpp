#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector e_monteCPP(const IntegerMatrix& x, const IntegerMatrix& rin, const double& K,
	const int& Kmin, const NumericVector& par1, const NumericVector& par2, const int& type){
 Function eFunc("e.main.monte", Rcpp::Environment::namespace_env("echelon"));

 int n = x.ncol();
 NumericVector result(n);
 for (int i = 0; i < n; ++i) {
  NumericVector cas = NumericVector(x(_, i).begin(), x(_, i).end());
  result[i] = as<double>(eFunc(cas, rin, K, Kmin, par1, par2, type));
 }
 return result;
}
