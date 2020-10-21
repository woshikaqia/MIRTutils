#' Using Rcpp to create a specialized matrix multiplication function used in person_fit()
#' @description Create specialized matrix multiplication function used in person_fit() while calculating the varaince of
#' loglikelihodd for cluster items
#' @return a function named \code{mmult()}
#' @author Zhongtian Lin lzt713@gmail.com
#' @export

func <- 'NumericMatrix mmult( NumericMatrix m , NumericVector v) {
  if( m.ncol() != v.size() ) stop("Non-conformable arrays") ;
  NumericMatrix out(m) ;
  for (int j = 0; j < m.ncol(); j++) {
    for (int i = 0; i < m.nrow(); i++) {
      out(i,j) = m(i,j) * v[j];
    }
  }
  return out ;
}'
#  Make it available
Rcpp::cppFunction( func )
