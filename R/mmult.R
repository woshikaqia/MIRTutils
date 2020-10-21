# specilzied multiplication function used in person_fit()
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
