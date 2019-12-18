#' RMSE
#'
#' @export
#' 


calc_RMSE <- function(obs,sim){
  
  N <- length(obs)
  RMSE <- sqrt(sum((obs-sim)^2)/(N))
  return(RMSE)
}

#' bias
#'
#' @export
#' 
#' 
calc_bias <- function(obs,sim){
  
  N <- length(obs)
  bias <- sum(sim-obs)/(N)
  return(bias)
}

#' SEPC
#'
#' @export
#' 
calc_SEPC <- function(obs,sim){
  
  N <- length(obs)
  bias <- calc_bias(obs,sim)
  SEPC <- sqrt(sum((sim-obs-bias)^2)/(N))
  return(SEPC)
}


#' All metrics
#'
#' @export
#' 
calc_metrics <- function(obs,sim){
  
  mat <- cbind(obs,sim)
  mat <- mat[!apply(is.na(mat),1,any),]
  obs <- mat[,1]
  sim <- mat[,2]
  
  return(c(rmse = calc_RMSE(obs,sim), bias = calc_bias(obs,sim), SEPC = calc_SEPC(obs,sim)))
}