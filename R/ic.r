##############################################
# calculcate the information criteria
##############################################

get.ic.from.mse <- function(mse, n, p, m, cp.penalty = 1, ic){
  ic.factor <- NULL
  if(ic == "bic"){
    ic.factor <- log(n)
  } else if(ic == "aic") {
    ic.factor <- 2  
  } 
  effective.param <-(p+m*cp.penalty)
  return (log(mse*n/(n-effective.param))*n+effective.param*ic.factor)
}