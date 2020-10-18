##############################################
# helper functions 
##############################################
get.all.not.ths.n.cps <- function(not.mod) {
  all.cps <- append(0, not.mod$solution.path$cpt) # 0 means no change point, and 
  # all.cps[[length(all.cps)]] <- NULL # do not consider the last set of change points
  max.th <- 100000 # hard code for now
  all.ths <- c(max.th, not.mod$solution.path$th)
  return (list(all.cps = all.cps,
               all.ths = all.ths))
}

get.left.n.right.end.pts <- function(cps, n) {
  if(any(cps < 0 | cps > n)) {
    stop("change points should not be smaller than 0 or larger than n")
  }
  # makes sure change points are sorted
  cps <- sort(cps)
  # remove 0 
  right.end.pt <- cps[which(cps!= 0)]
  # add n 
  if( !(n%in%right.end.pt) ) {
    right.end.pt <- c(right.end.pt, n)
  }
  num.end.pt <- length(right.end.pt)
  left.end.pt <- c(1,(cps[-num.end.pt]+1))
  return (list(left = left.end.pt,
               right = right.end.pt))
}

mse.from.pred <- function(pred, test.d) {
  return (mean((pred-test.d)^2))
}