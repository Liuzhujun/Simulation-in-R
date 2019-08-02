#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector col_mean(NumericMatrix X) {
     int m = X.nrow();
     int n = X.ncol();
     NumericVector out(n);
     for(int i = 0; i < m; ++i) {
       out = out + X(i,_);
     }
     out = out / m;
     return out;
   }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
mean_R <- function(X){
  n = dim(X)[1]
  m = dim(X)[2]
  cm <- rep(0,m)
  for (i in 1:n) {
       cm = cm +X[i,]
  }
  return( cm / n)
}
*/
