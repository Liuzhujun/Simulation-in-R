sim_single <- function(n,beta,SLP=FALSE,mydist,...){
  p = length(beta)
  e=mydist(n,...)
  X = matrix(rnorm(n*p),n,p)
  Y = X%*%beta0+e
  if(SLP) {
for (i in 1:p) {
    solve(matrix(rnorm(n*n),n,n)+diag(n))

}
    }
  return(as.numeric(lm(Y~X+0)$coef))
}
