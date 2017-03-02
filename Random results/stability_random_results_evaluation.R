##############################################
# Evaluation of the random search experiments
##############################################


# show Pareto fronts
library(data.table)
library(mco)
library(ggplot2)
load("results_stability_random_all.RData")

# 1- stability for minimization
results_all[, 16:37] <- 1 - results_all[, 16:37]
colnames(results_all)[16:37] <- paste("1 -", colnames(results_all)[16:37])

measure_names <- paste("1 -", c("SC", "SD-0", "SD-1", "SD-2", "SD-10", "SD", "SJ", "SL", "SN", "SO", "SS"))

# Pareto front of criteria first and second with the constraint that
# the value of third is not greater than the best achieved value of third
# + third_epsilon
# returns Pareto front
pareto_front_2d <- function(dataset, first, second, third, third_epsilon)
{
  data <- split(results_all, results_all$prob)[[dataset]]
  sub <- data[, c(first, second, third)]
  
  nas <- apply(sub, 1, function(x) any(is.na(x)))
  sub <- sub[!nas, ]
  
  best_third <- min(sub[, 3])
  sub_best <- sub[sub[, 3] <= best_third + third_epsilon, 1:2]
  pareto_front <- paretoFilter(as.matrix(sub_best))
  rows_pf <- rownames(pareto_front)
  pf_all <- data[rownames(data) %in% rows_pf, ]
  return(pf_all)
}

# paretos_stab_size_BO <- lapply(1:11, function(i){
#   pareto_front_2d("AP_Breast_Ovary", "mean_size_train", colnames(results_all)[16:26][i],
#                   "performance_train.mmce.test.mean", 0.05)
# })
# 
# paretos_stab_size_CK <- lapply(1:11, function(i){
#   pareto_front_2d("AP_Colon_Kidney", "mean_size_train", colnames(results_all)[16:26][i],
#                   "performance_train.mmce.test.mean", 0.05)
# })
# 
# paretos_stab_size_St <- lapply(1:11, function(i){
#   pareto_front_2d("Stomach", "mean_size_train", colnames(results_all)[16:26][i],
#                   "performance_train.mmce.test.mean", 0.05)
# })
# 
# save(paretos_stab_size_BO, paretos_stab_size_CK, paretos_stab_size_St,
#      file = "pareto_fronts_random.RData")

load("pareto_fronts_random.RData")


# joint pareto plots per data set
joint_pareto <- function(paretolist, choice_indices, title, xnames, yname, third, 
                         xlab, ylab, legendlab, add_rows = NULL, log = FALSE)
{
  paretolist <- paretolist[choice_indices]
  # xnames <- xnames[choice_indices]
  
  ggdata <- lapply(1:length(paretolist), function(i)
  {
    pf <- paretolist[[i]]
    
    # if addrows != NULL non Pareto optimal rows are added
    if(!is.null(add_rows))
    {
      optimal <- c(rep("Yes", nrow(pf)), rep("No", nrow(add_rows)))
      pf <- rbind(pf, add_rows)
      pf <- cbind(pf, optimal = optimal)
    }
    else
    {
      pf <- cbind(pf, optimal = rep("Yes", nrow(pf)))
    }
    
    colnames(pf)[16:26] <- measure_names
    melt(pf, id.vars = c(yname, third, "optimal"), measure.vars = measure_names[choice_indices][i])
  })
  
  alldata <- Reduce(rbind, ggdata)
  
  if(!(is.null(add_rows)))
  {
    ggall <- ggplot(alldata, aes(y = get(yname), x = value, color = get(third), group = "optimal")) +
      geom_point(aes(shape = optimal), size = 4, show.legend = TRUE) + 
      scale_shape_manual(values = c(17, 19), name = "", labels = c("Accuracy\noptimal", "Pareto\noptimal")) + 
      guides(shape = guide_legend(order = 1))
  }
  else
  {
    ggall <- ggplot(alldata, aes(y = get(yname), x = value, color = get(third))) + 
      geom_point(size = 4)
  }
  ggall <- ggall + labs(x = xlab, y = ylab, title = title) + 
    theme(legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          axis.title = element_text(size = 17),
          axis.text = element_text(size = 14),
          strip.text.x = element_text(size = 15),
          title = element_text(size = 14),
          legend.position = "bottom",
          legend.key.width = unit(1.25,"cm")) +
    scale_colour_gradient(name = legendlab, low = "darkred", high = "lightgoldenrod1",
                          guide = guide_colourbar(title.vjust = 0.8)) +
    guides(fill = guide_legend(order = 0)) + 
    facet_wrap( ~ variable, ncol = 2)
  if(log) ggall <- ggall + scale_y_log10(breaks = 10^(1:4))
  
  ggall
}

# determines methods with highest accuracy on training data
best_accuracy <- function(dataset)
{
  data <- results_all[results_all$prob == dataset, ]
  best_acc <- min(data$performance_train.mmce.test.mean)
  data[data$performance_train.mmce.test.mean == best_acc, ]
}

pdf("..\\Plots\\pareto_stab_size_joint.pdf", height = 6.5, width = 8)
print(joint_pareto(paretos_stab_size_BO, c(7, 3, 1, 5), "AP_Breast_Ovary", 
                   colnames(results_all)[16:26], "mean_size_train", "performance_train.mmce.test.mean",
                   "1 - Stability value", "Mean size of chosen features sets", "Error", 
                   best_accuracy("AP_Breast_Ovary"), log = TRUE))
print(joint_pareto(paretos_stab_size_CK, c(7, 3, 1, 5), "AP_Colon_Kidney", 
                   colnames(results_all)[16:26], "mean_size_train", "performance_train.mmce.test.mean",
                   "1 - Stability value", "Mean size of chosen features sets", "Error",
                   best_accuracy("AP_Colon_Kidney"), log = TRUE))
print(joint_pareto(paretos_stab_size_St, c(7, 3, 1, 5), "Stomach", 
                   colnames(results_all)[16:26], "mean_size_train", "performance_train.mmce.test.mean",
                   "1 - Stability value", "Mean size of chosen features sets", "Error",
                   best_accuracy("Stomach"), log = TRUE))
dev.off()


# plot stability and mean size of the best (min_error + third_eps) methods
plot_best <- function(dataset, xnames, yname, third, third_eps, xlab, ylab, legendlab, log = FALSE)
{
  data <- results_all[results_all$prob == dataset, ]
  data_best <- min(data[, third])
  data <- data[data[, third] <= data_best + third_eps, ]
  colnames(data)[16:26] <- measure_names
  
  # radom order
  set.seed(1234)
  data <- data[sample(nrow(data)), ]
  
  plot_data <- NULL
  for(xn in xnames)
  {
    melt <- melt(data, id.vars = c(yname, "algo"), measure.vars = xn)
    plot_data <- rbind(plot_data, melt)
  }

  ggall <- ggplot(plot_data, aes(y = get(yname), x = value, color = algo)) + 
    geom_point(size = 1) + 
    labs(x = xlab, y = ylab, title = dataset) + 
    theme(legend.title = element_text(size = 17),
          legend.text = element_text(size = 17),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          title = element_text(size = 15)) +
    scale_color_manual(values = 
                         c("firebrick1", "firebrick3", "firebrick4",
                           "steelblue1", "steelblue3", "steelblue4", 
                           "darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4",
                           "mediumorchid1", "mediumorchid3", "mediumorchid4"), 
                       name = legendlab,
                       labels = 
                         c("GLM Boosting & AUC", "GLM Boosting & MRMR", "GLM Boosting & Variance",
                           "Lasso Log.Reg. & AUC", "Lasso Log.Reg & MRMR", "Lasso Log.Reg & Variance",
                           "Random Forest & AUC", "Random Forest & MRMR", "Random Forest & Variance",
                           "SVM & AUC", "SVM & MRMR", "SVM & Variance")) + 
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
    facet_wrap( ~ variable, ncol = 2)
  if(log) ggall <- ggall + scale_y_log10(breaks = 10^(1:4))
  
  ggall
}

pdf("..\\Plots\\best_BO.pdf", height = 6, width = 13.5)
plot_best("AP_Breast_Ovary", measure_names[c(7, 3, 1, 5)], "mean_size_train", "performance_train.mmce.test.mean", 0.05,
          "1 - Stability value", "Mean size of chosen features sets", "Wrapper method")
dev.off()


# organize Pareto fronts and sets and eventually print them in a txt file
pareto_relevant <- function(paretolist, choice_indices = 1:length(paretolist))
{
  paretolist <- paretolist[choice_indices]
  for(i in 1:length(choice_indices))
  {
    j <- choice_indices[i]
    paretolist[[i]] <- paretolist[[i]][, c(1:4, 6:7, 9:15, 15 + j, 26 + j)]
  }
  paretolist
}

pf_BO <- pareto_relevant(paretos_stab_size_BO, c(7, 3, 1, 5))
pf_CK  <- pareto_relevant(paretos_stab_size_CK, c(7, 3, 1, 5))
pf_St <- pareto_relevant(paretos_stab_size_St, c(7, 3, 1, 5))

sink("pareto_fronts.txt")
pf_BO
pf_CK
pf_St
sink()

# pretty Pareto fronts
pareto_tables <- function(paretorel, n)
{
  tab <- paretorel[[1]]
  if(length(paretorel) > 1)
  {
    for(i in 2:length(paretorel))
    {
      tab <- cbind(tab, NA, NA)
      new_cols <- ncol(tab) - 1:0
      colnames(tab)[new_cols] <- colnames(paretorel[[i]])[14:15]
      
      old_indices <- which(rownames(paretorel[[i]]) %in% rownames(tab))
      if(length(old_indices) > 0)
      {
        for(j in old_indices)
        {
          pos <- which(rownames(tab) == rownames(paretorel[[i]])[j])
          tab[pos, new_cols] <- paretorel[[i]][j, 14:15]
        }
      }
      
      new_indices <- which(!(rownames(paretorel[[i]]) %in% rownames(tab)))
      if(length(new_indices) > 0)
      {
        new_part <- paretorel[[i]][new_indices, ]
        na_cols <- setdiff(1:ncol(tab), c(1:13, new_cols))
        if(length(na_cols) > 0)
        {
          new_part <- cbind(new_part[, 1:13], 
                            matrix(NA, nrow = nrow(new_part), ncol = length(na_cols)),
                            new_part[, 14:15])
        }
        tab <- rbind(tab, new_part)
      }
    }
  }
  
  
  tab$fw.perc <- tab$fw.perc * n
  tab[, -(1:2)] <- round(tab[, -(1:2)], 3)
  tab_order <- c(2, 4, 3, 9:13, 5, 7, 6, 8, 14:15, 18:19, 16:17, 20:21)
  tab <- tab[, tab_order]
  colnames(tab) <- c("Method", "n.feats", "lambda", "mstop", "C", "sigma", "num.trees", "min.node.size", 
                     "Error (train)", "Error (test)", "Size (train)", "Size (test)", "1 - SJ (train)",
                     "1 - SJ (test)", "1 - SC (train)", "1 - SC (test)", "1 - SD-1 (train)", "1 - SD-1 (test)",
                     "1 - SD-10 (train)", "1 - SD-10 (test)")
  tab
}

sink("pareto_tables.txt")
pareto_tables(pf_BO, 10935)
pareto_tables(pf_CK, 10935)
pareto_tables(pf_St, 10000)
sink()

best_acc_nicer <- function(tab, n)
{
  tab$fw.perc <- tab$fw.perc * n
  tab[, -(1:2)] <- round(tab[, -(1:2)], 3)
  tab_order <- c(2, 4, 3, 11:15, 6, 9, 7, 10, 22, 33, 16, 27, 18, 29, 20, 31)
  tab <- tab[, tab_order]
  colnames(tab) <- c("Method", "n.feats", "lambda", "mstop", "C", "sigma", "num.trees", "min.node.size",
                     "Error (train)", "Error (test)", "Size (train)", "Size (test)", "1 - SJ (train)",
                     "1 - SJ (test)", "1 - SC (train)", "1 - SC (test)", "1 - SD-1 (train)", "1 - SD-1 (test)",
                     "1 - SD-10 (train)", "1 - SD-10 (test)")
  tab
}

sink("best_accuracies.txt")
best_acc_nicer(best_accuracy("AP_Breast_Ovary"), 10935)
best_acc_nicer(best_accuracy("AP_Colon_Kidney"), 10935)
best_acc_nicer(best_accuracy("Stomach"), 10000)
sink()

# features of Pareto optimal methods for AP_Colon_Kidney
load("features_train_stability_random.RData")
rows_CK <- rownames(pareto_tables(pf_CK, 10935))
feat_rows_CK <- features_train[rows_CK]
feats_CK_unique <- lapply(feat_rows_CK, function(x) unique(unlist(x)))
any(lengths(feats_CK_unique) > 1) # always one feature
feats_CK <- table(unlist(feats_CK_unique))
feats_CK
# V10841  V1448  V8237 
#      1      9     57 

# orignal names of genes for AP_Colon_Kidney
data_CK <- foreign::read.arff("..\\Data\\AP_Colon_Kidney.arff")
data_CK <- subset(data_CK, select = -c(Tissue, ID_REF))
feats_CK_pos <- as.numeric(substring(names(feats_CK), first  = 2))
feats_orig <- colnames(data_CK)[feats_CK_pos]
feats_orig
# "46323_at"    "201839_s_at" "224596_at"  

load("..\\Data\\AP_Colon_Kidney.RData")
ggboxplot <- function(gene, gene_name)
{
  ggplot(AP_Colon_Kidney, aes(x = target, y = get(gene))) + 
    geom_boxplot() +
    labs(y = gene_name, x = "Class") +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 15))
}

pdf(file = "..\\Plots\\boxplot_genes.pdf", width = 5, height = 5)
print(ggboxplot(names(feats_CK)[1], feats_orig[1]))
print(ggboxplot(names(feats_CK)[2], feats_orig[2]))
print(ggboxplot(names(feats_CK)[3], feats_orig[3]))
dev.off()


AUC_scores <- function(x) {
  pred <- ROCR::prediction(predictions = x, labels = AP_Colon_Kidney$target)
  score <- ROCR::performance(pred, "auc")@y.values[[1L]]
  abs(0.5 - score)
}

aucs <- apply(AP_Colon_Kidney[, -1], 2, AUC_scores)
auc_ranks <- rank(aucs)

10936 - auc_ranks[c("V10841", "V1448", "V8237")] # make 1 the best value
