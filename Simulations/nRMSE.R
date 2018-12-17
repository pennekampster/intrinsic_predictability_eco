nRMSE <- function(x, step_ahead){
  #browser()
  x <- x[, c("obs", "pred")] #now x is two columns, first is observed, second is predicted
  x <- x[complete.cases(x), ] #now cases with no predictions are removed
  
  RMSE <- sqrt(sum((x[,1]-x[,2])^2)/step_ahead)
  nRMSE <- RMSE/(max(x[,1]-min(x[,1])))
  
  return(nRMSE)
}

RMSE <- function(x, step_ahead){
  #browser()
  x <- x[, c("obs", "pred")] #now x is two columns, first is observed, second is predicted
  x <- x[complete.cases(x), ] #now cases with no predictions are removed
  
  RMSE <- sqrt(sum((x[,1]-x[,2])^2)/step_ahead)

  return(RMSE)
}