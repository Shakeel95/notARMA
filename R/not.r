# ========== NOT functions =========
get.not.mod <- function(d, contrast) {
  set.seed(1)
  not.mod <- not(d, contrast = contrast, rand.intervals = FALSE, intervals = get.all.not.intervals(length(d)))
  return (not.mod)
}

get.all.not.intervals <- function(n) {
  M <- (n-1)*n/2
  ind <- matrix(0, M, 2)
  ind[,1] <- rep(1:(n-1), (n-1):1)
  ind[,2] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)
  return (ind)
}