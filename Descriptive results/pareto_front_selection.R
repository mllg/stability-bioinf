######################################################
# Implementation of the epsilon constraint selection
######################################################


# iterative algorithm for selecting a good item out of a pareto front
# input:
# data: matrix or data.frame containing information about the items
# criteria_indices: numeric vector containing the indices of columns 
#                   of the criteria used for pareto optimality in order
#                   of importance
#                   note: all criteria must be minimized
# criteria_tol: numeric vector containing the tolerances for the corresponding indices
# delta: tolerance reducing factor
# verbose: logical, should information be printed?
pfs <- function(data, criteria_indices, criteria_tol = rep(0.05, length(criteria_indices)),
                delta = 0.9, verbose = FALSE)
{
  rows <- 1:nrow(data)
  
  stop_now <- function()
  {
    crit1 <- length(rows) == 1
    uniques <- apply(data[rows, criteria_indices], 2, function(x) length(unique(x)))
    crit2 <- max(uniques) == 1
    return(crit1 || crit2)
  }
  
  while(!stop_now())
  {
    for(i in 1:length(criteria_indices))
    {
      # best performance in criterion with index i
      best_performance_i <- min(data[rows, criteria_indices[i]])
      
      # keep only those items that lie within the tolerance
      keep <- data[rows, criteria_indices[i]] <= best_performance_i + criteria_tol[i]
      rows <- rows[keep]
      
      if(verbose) print(c(i, best_performance_i + criteria_tol[i]))
      
      if(stop_now()) break
    }
    criteria_tol <- criteria_tol * delta
  }
  return(rows)
}
