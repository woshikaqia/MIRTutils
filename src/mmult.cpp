#include <Rcpp.h>
using namespace Rcpp;

//' A specialized matrix multiplication
//' @description A specialized matrix multiplication function used in person_fit() while calculating the variance of
//' loglikelihodd for cluster items
//' @param m, a numeric matrix
//' @param v, a numeric vector
//' @author Zhongtian Lin lzt713@gmail.com
//' @export
// [[Rcpp::export]]
NumericMatrix mmult(NumericMatrix m , NumericVector v) {
  if( m.ncol() != v.size() ) stop("Non-conformable arrays") ;
  NumericMatrix out(m) ;
  for (int j = 0; j < m.ncol(); j++) {
    for (int i = 0; i < m.nrow(); i++) {
      out(i,j) = m(i,j) * v[j];
    }
  }
  return out ;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

