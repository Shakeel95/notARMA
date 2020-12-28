##############################################
# get ARMA-NOT-1 (or 2) with threshold selected by BIC
##############################################
get.not.arma.bic.mod <- function(not.mod, is.detrend, return.all.mods = FALSE) {
  all.not.ths.n.cps <- get.all.not.ths.n.cps(not.mod)
  all.cps <- all.not.ths.n.cps$all.cps
  all.ths <- all.not.ths.n.cps$all.ths
  
  # storage for return 
  result <- list(mod = NULL,
                 th = NULL,
                 cps = c(),
                 bic = Inf)
  if(return.all.mods) {
    result$all.mods <- list()
    result$all.ths <- c()
  }
  
  # only consider the first 30 (or more if original NOT selects more)
  count <- 0
  max.total.count <- min(length(all.cps), max(which(all.ths == not::features(not.mod)$th), 30))
  max.count <- 30
  for(i in 1:max.total.count) {
    cps <- sort(all.cps[[i]])  
    if(is.detrend) {
      current.mod <- get.not.arma.2.given.cps(not.mod = not.mod, cps = cps)
    } else{
      current.mod <- get.not.arma.1.given.cps(not.mod = not.mod, cps = cps)
    }
    current.bic <- get.bic(current.mod, not.mod)
    current.th <- all.ths[i] 
    
    # update if bic is smaller than previous one
    if(current.bic < result$bic) {
      count <- 0
      result$mod <- current.mod
      result$th <- current.th
      result$bic <- current.bic
    } else {
      count <- count + 1
    }
    
    # record the mods if needed
    if(return.all.mods) {
      result$all.mods[[i]] <- current.mod
      result$all.ths[i] <- current.th
    }
    
    # exit if try more than max count
    if(count >= max.count & is.finite(result$bic)) {
      break
    }
  }
  
  # return results
  if(length(result$mod) > 1) {
    result$cps <- result$arma.mod[[length(result$mod)-1]]$end
  }
  return (result)
}