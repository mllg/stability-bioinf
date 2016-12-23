##############################################
# Evaluation of the random search experiments
# (includes creation of tables)
##############################################


library(data.table)
library(xtable)

load("results_stability_random.RData")
method <- strsplit(as.character(results$algo), "_")
algo <- factor(sapply(method, function(x) ifelse(length(x) == 2, x[1], paste(x[1], x[2], sep = "_"))),
               levels = c("glmboost", "logreg", "random_forest", "svm"))
filter <- factor(sapply(method, function(x) x[length(x)]),
                 levels = c("auc", "fmrmr", "variance"))

results <- cbind(prob = results[, 1], algo = algo, fw.method = filter, results[, 3:ncol(results)])


summary(results)


# calculate 1- stability for minimization
results_flip <- results
results_flip$stability_train <- 1 - results_flip$stability_train
stability_index <- which(colnames(results) == "stability_train")
colnames(results_flip)[stability_index] <- "1 - stability_train"



########## optimal models ##########################


# only one optimal model per data set
source("..\\Descriptive results\\pareto_front_selection.R")

results_flip_split <- split(results_flip, results_flip$prob)

percentages <- data.frame(performance = c(0.05, 0.025, 0.01, 0.05, 0.05),
                          stability = c(0.05, 0.05, 0.05, 0.05, 0.075),
                          size = c(0.05, 0.05, 0.05, 0.075, 0.1))

# AP_Colon_Kidney
best_indices_CK <- apply(percentages, 1, function(p) pfs(
  data = results_flip_split$AP_Colon_Kidney,
  criteria_indices = c(which(colnames(results) == "performance_train.mmce.test.mean"), 
                       which(colnames(results) == "stability_train"),
                       which(colnames(results) == "mean_size_train")),
  criteria_tol = p
))

best_CK <- lapply(best_indices_CK, function(i) results_flip_split$AP_Colon_Kidney[i,])



# best classification accuracy
best_index_acc_CK <- which(results_flip_split$AP_Colon_Kidney$performance_train.mmce.test.mean ==
                             min(results_flip_split$AP_Colon_Kidney$performance_train.mmce.test.mean)) 

best_acc_CK <- results_flip_split$AP_Colon_Kidney[best_index_acc_CK, ]


# AP_Breast_Ovary
best_indices_BO <- apply(percentages, 1, function(p) pfs(
  data = results_flip_split$AP_Breast_Ovary,
  criteria_indices = c(which(colnames(results) == "performance_train.mmce.test.mean"), 
                       which(colnames(results) == "stability_train"),
                       which(colnames(results) == "mean_size_train")),
  criteria_tol = p
))

best_BO <- lapply(best_indices_BO, function(i) results_flip_split$AP_Breast_Ovary[i,])


# best classification accuracy
best_index_acc_BO <- which(results_flip_split$AP_Breast_Ovary$performance_train.mmce.test.mean ==
                             min(results_flip_split$AP_Breast_Ovary$performance_train.mmce.test.mean)) 

best_acc_BO <- results_flip_split$AP_Breast_Ovary[best_index_acc_BO, ]



# Stomach
best_indices_S <- apply(percentages, 1, function(p) pfs(
  data = results_flip_split$Stomach,
  criteria_indices = c(which(colnames(results) == "performance_train.mmce.test.mean"), 
                       which(colnames(results) == "stability_train"),
                       which(colnames(results) == "mean_size_train")),
  criteria_tol = p
))

best_S <- lapply(best_indices_S, function(i) results_flip_split$Stomach[i,])


# best classification accuracy
best_index_acc_S <- which(results_flip_split$Stomach$performance_train.mmce.test.mean ==
                             min(results_flip_split$Stomach$performance_train.mmce.test.mean)) 

best_acc_S <- results_flip_split$Stomach[best_index_acc_S, ]



# table with models and table with performance criteria
tables <- function(percentages, best, best_acc, n_vars, digits = 3)
{
  # make data frame out of lists with best methods
  tab <- best_acc
  for(i in 1:length(best))
  {
    tab <- rbind(tab, best[[i]])
  }
  
  # include epsilon values
  perc <- matrix(0, ncol = ncol(percentages), nrow = nrow(tab))
  colnames(perc) <- colnames(percentages)
  perc <- as.data.frame(perc)
  
  # only accuracy
  n_best_acc <- nrow(best_acc)
  perc$stability[1:n_best_acc] <- Inf
  perc$size[1:n_best_acc] <- Inf
  
  # real epsion values
  row_counter <- n_best_acc
  for(i in 1:length(best))
  {
    n_i <- nrow(best[[i]])
    perc[row_counter + 1:n_i, ] <- percentages[i, ]
    row_counter <- row_counter + n_i
  }
  
  tab <- cbind(perc, tab)
  
  # fw.abs instead of fw.perc
  tab$fw.perc <- round(tab$fw.perc * n_vars, 0)
  colnames(tab)[colnames(tab) == "fw.perc"] <- "fw.abs"
  
  # stability instead of 1 - stability
  tab$'1 - stability_train' <- 1 - tab$'1 - stability_train' 
  colnames(tab)[colnames(tab) == "1 - stability_train"] <- "stability_train"

  tab_mod <- tab[,  c(5, 8, 7, 6, 15:19)]
  
  tab_perf <- tab[, c(1:3, 10, 9, 11, 13, 12, 14)]
  tab_perf <- round(tab_perf , digits)
  
  return(list(all = tab, methods = tab_mod, performances = tab_perf))
}


tab_BO <- tables(percentages, best_BO, best_acc_BO, n_vars = 10935)
tab_CK <- tables(percentages, best_CK, best_acc_CK, n_vars = 10935)
tab_S <- tables(percentages, best_S, best_acc_S, n_vars = 10000)


# give same IDs to identical wrapper methods
ID_finder <- function(data, offset = 0)
{
  ID <- rep(0L, nrow(data))
  dup <- duplicated(data)
  ID[!dup] <- as.integer(1:sum(!dup) + offset)
  
  # no duplicates
  if(sum(dup) == 0) return(ID)
  
  for(i in which(dup))
  {
    for(j in 1:(i-1))
    {
      if(all(data[i, ] == data[j, ], na.rm = TRUE))
      {
        ID[i] <- ID[j]
        break
      }
    }
  }
  return(ID)
}

ID_BO <- ID_finder(tab_BO$methods, 0)
ID_CK <- ID_finder(tab_CK$methods, max(ID_BO))
ID_S <- ID_finder(tab_S$methods, max(ID_CK))

performances <- rbind(tab_BO$performances, tab_CK$performances, tab_S$performances)
methods <- rbind(tab_BO$methods, tab_CK$methods, tab_S$methods)

# joint IDs for all data sets
performances <- cbind(ID = c(ID_BO, ID_CK, ID_S), performances)
methods$ID <- cbind(ID = c(ID_BO, ID_CK, ID_S), methods)

# performance table
print(xtable(performances, digits = 3), include.rownames = FALSE)

# methods table
methods2 <- data.frame(ID = c(ID_BO, ID_CK, ID_S), 
                       Classifier = c("GLM boosting", "Lasso Log. Reg.", "Random Forest", "SVM")[as.numeric(methods$algo)],
                       Parameters = ifelse(methods$algo == "glmboost", paste0("mstop = ", methods$mstop),
                                           ifelse(methods$algo == "logreg", paste0("lambda = ", round(methods$lambda, 5)),
                                                  ifelse(methods$algo == "random_forest", 
                                                         paste0("num.trees = ", methods$num.trees, ", min.node.size = ", methods$min.node.size),
                                                         paste0("sigma = ", methods$sigma, ", C = ", methods$C)))),
                       Filter = c("Variance", "AUC", "MRMR")[as.numeric(methods$fw.method)],
                       fw.abs = as.integer(methods$fw.abs)
                       )
methods2 <- methods2[!duplicated(methods2),]

print(xtable(methods2, digits = 5), include.rownames = FALSE)
