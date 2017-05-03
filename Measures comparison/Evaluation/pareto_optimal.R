#######################################################
# Determine Pareto optimal configurations for the 
# comparison of stability measures
######################################################


library(mco)
library(plyr)
library(checkmate)
library(data.table)


# returns indices of observations that are pareto optimal
# splits: named list of splits to be made
pareto_front_index <- function(data, targets, splits = NULL)
{
  # check types
  assertDataFrame(data, min.rows = 1, min.cols = 1)
  if(!is.null(splits)) assertList(splits, min.len = 1, any.missing = FALSE)
  assertCharacter(targets, min.len = 1, any.missing = FALSE)
  
  # check meaningful input
  assertSubset(targets, colnames(data))
  if(!is.null(splits)) assertNamed(splits)
  if(!is.null(splits)) assertSubset(names(splits), colnames(data))
  
  rownames(data) <- 1:nrow(data)
  
  # neccessary part of data matrix
  cols <- which(colnames(data) %in% targets)
  rows <- 1:nrow(data)
  if(!is.null(splits))
  {
    for(i in 1:length(splits))
    {
      index <- which(colnames(data) == names(splits)[i])
      rows_i <- which(data[, index] == splits[[i]])
      if(length(rows_i) == 0) warning(paste0("There are no observations with",
                                             names(splits)[i], " = ", splits[[i]], "!"))
      rows <- intersect(rows, rows_i)
    }
    if(length(rows) == 0) return(NULL)
  }
  
  data_part <- data[rows, cols, drop = FALSE]
  
  # remove rows with missing values
  rows_with_missing_values <- which(apply(data_part, 1, function(x) any(is.na(x))))
  if(length(rows_with_missing_values) > 0)
  {
    data_part <- data_part[-rows_with_missing_values, , drop = FALSE]
  }
  
  if(nrow(data_part) == 0)
  {
    return(NULL)
  }
  else
  {
    pareto_front <- paretoFilter(as.matrix(data_part))
    return(as.numeric(rownames(pareto_front)))
  }
}

pf <- function(dataset, target1, target2, target3)
{
  pfi <- pareto_front_index(data = results_flip, targets = c(target1, target2, target3),
                            splits = list(prob = dataset))
  return(list(pareto_indices = pfi, targets = c(target1, target2, target3), prob = dataset))
}


# Data preparation
load("results_measures_comparison.RData")
# indices of stability measures in result data
measures_indices <- which(colnames(results) %in%
                            c("cor_pearson", "davis_0", "davis_1", "davis_2", "davis_10", "dice", "jaccard",
                              "lustgarten", "novovicova", "ochiai", "somol", "zucknick"))

measures <- results[, measures_indices]
measure_names <- c("SC", "SD-0", "SD-1", "SD-2", "SD-10", "SD", "SJ", "SL", "SN", "SO", "SS", "SZ")

# calculate 1- stability and remove NAs for minimization
results_flip <- results
nas <- rowSums(is.na(measures)) > 0
results_flip[!nas, measures_indices] <- 1 - results_flip[!nas, measures_indices]
colnames(results_flip)[measures_indices] <- paste("1 -", colnames(results_flip)[measures_indices])





stabilities <- paste("1 -", c("cor_pearson", "davis_0", "davis_1", "davis_2", "davis_10", 
                              "dice", "jaccard", "lustgarten", 
                              "novovicova", "ochiai", "somol", "zucknick"))

# Pareto optimal results for AP_Colon_Kidney
results_pareto <- Map(pf, target1 = "mmce.test.mean", target2 = "mean_size", 
                      target3 = stabilities, dataset = "AP_Colon_Kidney")
names(results_pareto) <- 1:12
save(results_pareto, file = "results_pareto_CK.RData")

# Pareto optimal results for AP_Breast_Ovary
results_pareto <- Map(pf, target1 = "mmce.test.mean", target2 = "mean_size", 
                      target3 = stabilities, dataset = "AP_Breast_Ovary")
names(results_pareto) <- 1:12
save(results_pareto, file = "results_pareto_BO.RData")

# Pareto optimal results for Stomach
results_pareto <- Map(pf, target1 = "mmce.test.mean", target2 = "mean_size", 
                      target3 = stabilities, dataset = "Stomach")
names(results_pareto) <- 1:12
save(results_pareto, file = "results_pareto_St.RData")
