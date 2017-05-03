##############################################
# Evaluation of the experiments for the 
# random search for desirable configurations
##############################################


# show Pareto fronts
library(data.table)
library(mco)
library(ggplot2)
library(xtable)

if(!file.exists("..\\Plots")) dir.create("..\\Plots")


load("results_model_finding.RData")
results <- results[, c("prob", "algo", "fw.perc", "lambda", "mstop", "C", "sigma", "num.trees", "min.node.size",
                       "performance_train.mmce.test.mean", "mean_size_train",
                       "performance_test.mmce.test.mean", "mean_size_test",
                       paste0(c("cor_pearson", "davis_0", "davis_1", "davis_2", "davis_10", "dice",
                                "jaccard", "lustgarten", "novovicova", "ochiai", "somol", "zucknick"), "_train"),
                       paste0(c("cor_pearson", "davis_0", "davis_1", "davis_2", "davis_10", "dice",
                                "jaccard", "lustgarten", "novovicova", "ochiai", "somol", "zucknick"), "_test"))]


# 1- stability for minimization
results[, 14:37] <- 1 - results[, 14:37]
colnames(results)[14:37] <- paste("1 -", colnames(results)[14:37])

measure_names <- paste("1 -", c("SC", "SD-0", "SD-1", "SD-2", "SD-10", "SD", "SJ", "SL", "SN", "SO", "SS", "SZ"))

# Pareto front of criteria first and second with the constraint that
# the value of third is not greater than the best achieved value of third
# + third_epsilon
# returns Pareto front
# pareto_front_2d <- function(dataset, first, second, third, third_epsilon)
# {
#   data <- split(results, results$prob)[[dataset]]
#   sub <- data[, c(first, second, third)]
# 
#   nas <- apply(sub, 1, function(x) any(is.na(x)))
#   sub <- sub[!nas, ]
# 
#   best_third <- min(sub[, 3])
#   sub_best <- sub[sub[, 3] <= best_third + third_epsilon, 1:2]
#   pareto_front <- paretoFilter(as.matrix(sub_best))
#   rows_pf <- rownames(pareto_front)
#   pf_all <- data[rownames(data) %in% rows_pf, ]
#   return(pf_all)
# }
# 
# paretos_stab_size_BO <- lapply(1:12, function(i){
#   pareto_front_2d("AP_Breast_Ovary", "mean_size_train", colnames(results)[14:25][i],
#                   "performance_train.mmce.test.mean", 0.05)
# })
# 
# paretos_stab_size_CK <- lapply(1:12, function(i){
#   pareto_front_2d("AP_Colon_Kidney", "mean_size_train", colnames(results)[14:25][i],
#                   "performance_train.mmce.test.mean", 0.05)
# })
# 
# paretos_stab_size_St <- lapply(1:12, function(i){
#   pareto_front_2d("Stomach", "mean_size_train", colnames(results)[14:25][i],
#                   "performance_train.mmce.test.mean", 0.05)
# })
# 
# names(paretos_stab_size_BO) <- names(paretos_stab_size_CK) <- names(paretos_stab_size_St) <- colnames(results)[14:25]
# 
# save(paretos_stab_size_BO, paretos_stab_size_CK, paretos_stab_size_St,
#      file = "pareto_fronts_model_finding.RData")

load("pareto_fronts_model_finding.RData")


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
    
    colnames(pf)[14:25] <- measure_names
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
  data <- results[results$prob == dataset, ]
  best_acc <- min(data$performance_train.mmce.test.mean)
  data[data$performance_train.mmce.test.mean == best_acc, ]
}

pdf("..\\Plots\\pareto_stab_size_joint.pdf", height = 6.5, width = 8)
print(joint_pareto(paretos_stab_size_BO, c(7, 3, 1, 5), "AP_Breast_Ovary", 
                   colnames(results)[14:25], "mean_size_train", "performance_train.mmce.test.mean",
                   "1 - Stability value", "Mean number of chosen features", "Error", 
                   best_accuracy("AP_Breast_Ovary"), log = TRUE))
print(joint_pareto(paretos_stab_size_CK, c(7, 3, 1, 5), "AP_Colon_Kidney", 
                   colnames(results)[14:25], "mean_size_train", "performance_train.mmce.test.mean",
                   "1 - Stability value", "Mean number of chosen features", "Error",
                   best_accuracy("AP_Colon_Kidney"), log = TRUE))
print(joint_pareto(paretos_stab_size_St, c(7, 3, 1, 5), "Stomach", 
                   colnames(results)[14:25], "mean_size_train", "performance_train.mmce.test.mean",
                   "1 - Stability value", "Mean number of chosen features", "Error",
                   best_accuracy("Stomach"), log = TRUE))
dev.off()


# plot stability and mean size of the best (min_error + third_eps) methods
plot_best <- function(dataset, xnames, yname, third, third_eps, xlab, ylab, legendlab, log = FALSE)
{
  data <- results[results$prob == dataset, ]
  data_best <- min(data[, third])
  data <- data[data[, third] <= data_best + third_eps, ]
  colnames(data)[14:25] <- measure_names
  
  cat("Number of displayed methods: ")
  cat(nrow(data))
  
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
          "1 - Stability value", "Mean number of chosen features", "Augmented method")
dev.off()




########### Pareto fronts and best accuracies for all stability measures ################
organizer <- function(tab, n, upper_bound_mean_size_train, measures1, measures2, measures3,
                      names1, names2, names3)
{
  tab <- tab[order(tab$mean_size_train, decreasing = TRUE), ]
  if(!is.null(upper_bound_mean_size_train))
  {
    tab <- tab[tab$mean_size_train <= upper_bound_mean_size_train, ]
  }
  
  tab$fw.perc <- round(tab$fw.perc * n)
  tab <- cbind(ID = 1:nrow(tab), tab)
  parameters <- c("lambda", "mstop", "C", "sigma", "num.trees", "min.node.size")
  
  splitted <- strsplit(as.character(tab$algo), "_", fixed = TRUE)
  classifier <- unlist(lapply(splitted, function(x) paste(x[-length(x)], collapse = "_")))
  classifier[classifier == "glmboost"] <- "GLM Boosting"
  classifier[classifier == "logreg"] <- "Lasso. Log. Reg."
  classifier[classifier == "random_forest"] <- "Random Forest"
  classifier[classifier == "svm"] <- "SVM"
  
  filter <- unlist(lapply(splitted, function(x) paste(x[length(x)])))
  filter[filter == "auc"] <- "AUC"
  filter[filter == "fmrmr"] <- "MRMR"
  filter[filter == "variance"] <- "Variance"
  
  options(scipen = 100)  # omit exponential notation
  pars <- apply(tab[, parameters], 1, function(x){
    existing <- which(!is.na(x))
    tmp <- sapply(existing, function(y){
      value <- round(x[y], 3)
      r <- 3
      while(value == 0)
      {
        r <- r + 1
        value <- round(x[y], r)
      }
      
      par_names <- parameters[y]
      
      # nicer latex formatting
      par_names[par_names == "mstop"] <- "m_{\\text{stop}}"
      par_names[par_names == "lambda"] <- "\\lambda"
      par_names[par_names == "num.trees"] <- "\\texttt{num.trees}"
      par_names[par_names == "min.node.size"] <- "\\texttt{min.node.size}"
      par_names[par_names == "sigma"] <- "\\sigma"
      
      paste0("$", par_names, " = ", value, "$")
    })
    if(length(tmp) > 1)
    {
      tmp <- paste(tmp, collapse = ", ")
    }
    tmp
  })
  
  
  tab0 <- cbind(data.frame(ID = tab$ID), Classifier = classifier, Parameters = pars, Filter = filter,
                "\\texttt{n.feats}" = round(tab$fw.perc, 0), "Error (train)" = round(tab$performance_train.mmce.test.mean, 3),
                "Error (test)" = round(tab$performance_test.mmce.test.mean, 3), "Size (train)" = tab$mean_size_train,
                "Size (test)" = tab$mean_size_test)
  
  cols1 <- paste0(rep(paste0("1 - ", measures1), each = 2), rep(c("_train", "_test"), length(measures1)))
  tab1 <- tab[, c("ID", cols1)]
  tab1[, -1] <- round(tab1[, -1], 3)
  colnames(tab1)[-1] <- paste0(rep(paste0("1 - ", names1), each = 2), rep(c(" (train)", " (test)"), length(names1)))
  
  cols2 <- paste0(rep(paste0("1 - ", measures2), each = 2), rep(c("_train", "_test"), length(measures2)))
  tab2 <- tab[, c("ID", cols2)]
  tab2[, -1] <- round(tab2[, -1], 3)
  colnames(tab2)[-1] <- paste0(rep(paste0("1 - ", names2), each = 2), rep(c(" (train)", " (test)"), length(names2)))
  
  cols3 <- paste0(rep(paste0("1 - ", measures3), each = 2), rep(c("_train", "_test"), length(measures3)))
  tab3 <- tab[, c("ID", cols3)]
  tab3[, -1] <- round(tab3[, -1], 3)
  colnames(tab3)[-1] <- paste0(rep(paste0("1 - ", names3), each = 2), rep(c(" (train)", " (test)"), length(names3)))
  
  return(list(tab0, tab1, tab2, tab3))
}

all_tables <- function(paretolist, n, upper_bound_mean_size_train = NULL,
                       measures1 = c("jaccard", "cor_pearson", "davis_1", "davis_10"),
                       measures2 = c("davis_0", "davis_2", "dice", "lustgarten"),
                       measures3 = c("novovicova", "ochiai", "somol", "zucknick"),
                       names1 = c("SJ", "SC", "SD-1", "SD-10"),
                       names2 = c("SD-0", "SD-2", "SD", "SL"),
                       names3 = c("SN", "SO", "SS", "SZ"))
{
  for(i in 1:length(paretolist))
  {
    other_measures_train <- setdiff(paste0("1 - ", c(measures1, measures2, measures3), "_train"), names(paretolist)[i])
    other_measures_test <- paste0(unlist(strsplit(other_measures_train, "train", fixed = TRUE)), "test")
    drop_cols <- which(colnames(paretolist[[i]]) %in% c(other_measures_train, other_measures_test))
    paretolist[[i]] <- paretolist[[i]][, -drop_cols, drop = FALSE]
  }
  
  tab <- unique(Reduce(function(a, b) merge(a, b, all = TRUE, sort = FALSE), paretolist))
  organizer(tab, n, upper_bound_mean_size_train, measures1, measures2, measures3, names1, names2, names3)
}

all_BO <- all_tables(paretos_stab_size_BO, 10935)
all_CK <- all_tables(paretos_stab_size_CK, 10935)
all_St <- all_tables(paretos_stab_size_St, 10000)


all_accuracies <- function(dataname, n, upper_bound_mean_size_train = NULL,
                           measures1 = c("jaccard", "cor_pearson", "davis_1", "davis_10"),
                           measures2 = c("davis_0", "davis_2", "dice", "lustgarten"),
                           measures3 = c("novovicova", "ochiai", "somol", "zucknick"),
                           names1 = c("SJ", "SC", "SD-1", "SD-10"),
                           names2 = c("SD-0", "SD-2", "SD", "SL"),
                           names3 = c("SN", "SO", "SS", "SZ"))
{
  tab <- best_accuracy(dataname)
  organizer(tab, n, upper_bound_mean_size_train, measures1, measures2, measures3, names1, names2, names3)
}

acc_BO <- all_accuracies("AP_Breast_Ovary", 10935)
acc_CK <- all_accuracies("AP_Colon_Kidney", 10935)
acc_St <- all_accuracies("Stomach", 10000)

###########

# print as xtables
xtab <- function(data, data_long, begin_id = 1)
{
  all <- get(paste0("all_", data))
  acc <- get(paste0("acc_", data))
  
  both <- lapply(1:4, function(i) rbind(all[[i]], acc[[i]]))
  both <- lapply(1:4, function(i){
    both[[i]]$ID <- 1:nrow(both[[i]]) + begin_id - 1
    both[[i]]
  })
  
  opts <- list(include.rownames = FALSE, 
               NA.string = "---",
               hline.after = c(-1, 0, nrow(all[[1]]), nrow(both[[1]])),
               sanitize.text.function = function(x) x)
  
  labels <- paste(c("configs", paste0("stabilities", 1:3)), data, sep = "_")
  xt <- xtable(both[[1]], digits = c(0, 0, 0, 0, 0, 0, 3, 3, 1, 1),
               caption = paste0("All Pareto optimal configurations for data set ", data_long,
                                " (above the horizontal line) and accuracy optimal configurations (below the horizontal line). ",
                                "The ID column references Tables~\\ref{", labels[2], "}, \\ref{", labels[3], "}, and \\ref{", 
                                labels[4], "}."),
               label = labels[1],
               align = c(">{\\raggedright\\arraybackslash}p{0cm}",">{\\raggedleft\\arraybackslash}p{0.5cm}", 
                         ">{\\raggedright\\arraybackslash}p{2.5cm}", ">{\\raggedright\\arraybackslash}p{3cm}", 
                         ">{\\raggedright\\arraybackslash}p{1.5cm}", ">{\\raggedleft\\arraybackslash}p{1.25cm}", 
                         ">{\\raggedleft\\arraybackslash}p{1.25cm}", ">{\\raggedleft\\arraybackslash}p{1.25cm}", 
                         ">{\\raggedleft\\arraybackslash}p{1.25cm}", ">{\\raggedleft\\arraybackslash}p{1.25cm}"))

  opts2 <- c(list(x = xt), opts)
  do.call(print, opts2)
  
  invisible(lapply(2:4, function(i){
    colwidths <- c(0, 0.5, rep(1.3, ncol(both[[i]]) - 1))
    aligns <- paste0(">{\\raggedleft\\arraybackslash}p{", colwidths, "cm}")
    xt <- xtable(both[[i]], digits = c(0, 0, rep(3, ncol(both[[i]]) - 1)),
                 caption = paste0("Stability values for the Pareto optimal configurations in Table~\\ref{", labels[1],
                                  "}. The ID column references Table~\\ref{", labels[1], 
                                  "}. \\enquote{---} means that the configuration is not Pareto optimal for the corresponding stability measure."),
                 label = labels[i],
                 align = aligns)
    opts2 <- c(list(x = xt), opts)
    do.call(print, opts2)
  }))
  
  invisible()
}

sink("..\\Plots\\complete_pareto_tables.txt")
xtab("BO", "AP\\_Breast\\_Ovary", 1)
xtab("CK", "AP\\_Colon\\_Kidney", nrow(all_BO[[1]]) + nrow(acc_BO[[1]]) + 1)
xtab("St", "Stomach", nrow(all_BO[[1]]) + nrow(acc_BO[[1]]) + nrow(all_CK[[1]]) + nrow(acc_CK[[1]]) + 1)
sink()

#################################################################################################
# boxplot performance values CK
pareto_perf_CK <- cbind(all_CK[[1]][, 6:9], all_CK[[2]][, -1])
boxplot(pareto_perf_CK)

pp_CK_melt <- cbind(id = rownames(pareto_perf_CK), pareto_perf_CK)
pp_CK_melt <- melt(pp_CK_melt, id.vars = "id")


ggbox1 <- ggplot(pp_CK_melt, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(y = "Target criteria values", x = "Target criteria") +
  scale_x_discrete(expand = c(0, 0), 
                   labels = c("Error\n(train)", "Error\n(test)", "Size\n(train)", "Size\n(test)", 
                              "1 - SJ\n(train)", "1 - SJ\n(test)",
                              "1 - SC\n(train)", "1 - SC\n(test)", "1 - SD-1\n(train)",
                              "1 - SD-1\n(test)", "1 - SD-10\n(train)", "1 - SD-10\n(test)")) +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 10))


pdf("..\\Plots\\pareto_optimal_CK_perf.pdf", height = 3, width = 10)
print(ggbox1)
dev.off()



###################################################################################################

# features of Pareto optimal methods for AP_Colon_Kidney
load("features_train_model_finding.RData")
rows_CK <- rownames(paretos_stab_size_CK[[1]])
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

10935 - auc_ranks[c("V10841", "V1448", "V8237")] # make 1 the best value
