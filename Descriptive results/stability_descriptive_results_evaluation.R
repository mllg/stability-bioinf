##############################################
# Evaluation of the descriptive experiments
# (includes creation of plots)
##############################################



set.seed(123)

load("results_stability_descriptive.RData")
ck <- which(results$prob == "AP_Colon_Kidney")
results <- results[ck, ]

library(gplots)
library(ggplot2)
library(data.table)
library(GGally)

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
colnames(results_flip)[measures_indices] <- paste("1 -",colnames(results_flip)[measures_indices])

################# overview ######################
summary(measures)

# Boxplots of measure values
measures_melt <- cbind(id = rownames(measures), as.data.frame(measures))
measures_melt <- melt(measures_melt, id.vars = "id")
ggbox <- ggplot(measures_melt, aes(x = variable, y = value)) + 
  geom_boxplot() +
  labs(y = "Stability values", x = "Stability measures") +
  scale_x_discrete(expand = c(0, 0), labels = measure_names) +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 10))

pdf("..\\Plots\\measures_boxplots.pdf", width = 6.5, height = 3, useKerning = FALSE)
print(ggbox)
dev.off()  

#################### heatmap #######################
# Heatmaps
source("ggplot_functions.R")

cor <- cor(measures, use = "pairwise.complete.obs")
gghm <- heatmap_ggplot(cor, measure_names)

pdf("..\\Plots\\heatmap_all.pdf", height = 8, width = 8, useKerning = FALSE)
print(gghm)
dev.off()

# do this analysis seperately for all classifiers
measures_split <- split(measures, results$algo)
for(i in 1: length(measures_split))
{
  cor_split <- cor(measures_split[[i]], use = "pairwise.complete.obs")
  gghm <- heatmap_ggplot(cor_split, measure_names, names(measures_split)[i])
  print(gghm)
}

# do this analysis seperately for all classifiers and filter methods
measures_split <- split(measures, results$algo)
names(measures_split) <- c("GLM Boosting", "Lasso Log. Reg.", "Random Forest", "SVM")
filter_split <- split(results$fw.method, results$algo)

for(i in 1:length(measures_split))
{
  measures_split2 <- split(measures_split[[i]], filter_split[[i]])
  names(measures_split2) <- c("AUC", "MRMR", "Variance")
  
  for(j in 1:length(measures_split2))
  {
    cor_split <- cor(measures_split2[[j]], use = "pairwise.complete.obs")
    title <- paste0(names(measures_split)[i], " & ", names(measures_split2)[j])
    gghm <- heatmap_ggplot(cor_split, measure_names, title)

    file <- paste0("..\\Plots\\heatmap_", i,"_" , j, ".pdf")
    pdf(file, height = 8, width = 8, useKerning = FALSE)
    print(gghm)
    dev.off()
  }
}

########### analysis of different classifiers ###########

gg_scatter_matrix <- function(indices, title = element_blank())
{
  order <- c(9, 10, 6, 12, 7, 2, 3, 4, 5, 8, 1, 11)

  ms <- measures
  colnames(ms) <- measure_names
  data <- cbind(id = indices, ms[indices, order], size = results$mean_size[indices])
  
  # randomize order for plotting
  data <- data[sample(1:nrow(data)), ]
  
  measures_order <- measure_names[order]
  data_melt <- NULL
  for(i in 1:length(measures_order))
  {
    m <- melt(data, id.vars = c("id", "size", measures_order[i]), measure.vars = measures_order)
    m <- cbind(m[, 1:2], variable2 = colnames(m)[3], value2 = m[,3], m[, 4:5])
    data_melt <- rbind(data_melt, m)
  }

  scatter_matrix <- ggplot(data_melt, aes(x = value2, y = value, color = size)) +
    geom_point(size = 0.5) + 
    labs(x = "Stability value", y = "Stability value", title = title) + 
    facet_grid(variable ~ variable2, labeller = "label_parsed") + 
    scale_colour_gradient(high = "darkred", low = "lightgoldenrod1", limits = c(0, 10935), name = "Size") + 
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    coord_equal(ratio = 1)
  
  return(scatter_matrix)
}

pdf(file = "..\\Plots\\scatter_glmboost.pdf", width = 20, height = 20)
print(gg_scatter_matrix(which(results$algo == "glmboost"), "GLM Boosting"))
dev.off()

pdf(file = "..\\Plots\\scatter_logreg.pdf", width = 20, height = 20)
print(gg_scatter_matrix(which(results$algo == "logreg"), "Lasso Logistic Regression"))
dev.off()

pdf(file = "..\\Plots\\scatter_rf.pdf", width = 20, height = 20)
print(gg_scatter_matrix(which(results$algo == "random_forest"), "Random Forest"))
dev.off()

postscript(file = "..\\Plots\\scatter_svm.ps", width = 20, height = 20)
print(gg_scatter_matrix(which(results$algo == "svm"), "SVM"))
dev.off()

################### stability and size ###################
# correlations between measures and mean_size
apply(measures, 2, function(x) cor(x, results$mean_size, use = "pairwise.complete.obs"))

size_stability <- as.data.frame(measures)
colnames(size_stability) <- measure_names

# better order for the columns
c_order <-  c(9, 12, 3, 8, 10, 7, 4, 1, 6, 2, 5, 11)
size_stability <- size_stability[, c_order]

size_stability <- cbind(size = results$mean_size, size_stability, algo = results$algo)

# randomize order of rows
size_stability <- size_stability[sample(1:nrow(size_stability)), ]


size_stability <- melt(size_stability, id.vars = c("size", "algo"))

ggmatrix <- ggplot(size_stability, aes(x = size, y = value, color = algo)) +
  geom_point(size = 1) + 
  labs(x = "Mean size of chosen features sets", y = "Stability value") + 
  facet_wrap( ~ variable, ncol = 4, labeller = "label_parsed") + 
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        strip.text.x = element_text(size = 12)) + 
  scale_color_discrete(name = "Classifier", 
                       labels = c("GLM Boosting", "Lasso Log. Reg.", 
                                  "Random Forest", "SVM")) + 
  guides(colour = guide_legend(override.aes = list(size = 2)))

pdf("..\\Plots\\size_stability.pdf", height = 6, width = 7.5, useKerning = FALSE)
print(ggmatrix)
dev.off()

########## optimal models ##########################

# which models are Pareto optimal when considering different measures of stability
load("results_pareto_cK.RData")

# indices of methods that are pareto optimal for any stability measure
paretos <- sort(unique(unlist(lapply(results_pareto, function(x) x$pareto_indices))))

# 0-1-coding: which of the above methods (row) is pareto optimal 
# for which stability measure (col)
pareto_models <- matrix(0, nrow = length(paretos), ncol = length(results_pareto))
rownames(pareto_models) <- paretos
colnames(pareto_models) <- unlist(lapply(results_pareto, function(x) x$targets[3]))
for(i in 1:length(results_pareto))
{
  pareto_models[paretos %in% results_pareto[[i]]$pareto_indices, i] <- 1
}

ggbin <- binary_heatmap_ggplot(pareto_models, measure_names)

pdf("..\\Plots\\heatmap_optimal.pdf", height = 3, width = 7.25, useKerning = FALSE)
print(ggbin)
dev.off()

table(results[paretos, c(2,4)])





# only one optimal model per measure -> same for all measures
source("pareto_front_selection.R")

best <- lapply(1:length(measures_indices), function(i) 
  pfs(data = results_flip[results_pareto[[i]]$pareto_indices, ], 
      criteria_indices = c(which(colnames(results) == "mmce.test.mean"), measures_indices[i],
                          which(colnames(results) == "mean_size"))))
best_indices <- lapply(1:length(best), function(i) results_pareto[[i]]$pareto_indices[best[[i]]])

results_flip[best_indices[[1]],]

# find out features used by best model (takes a while)
# load("features_stability_descriptive.RData")
# features[best_indices[[1]]]


# best classification accuracy
summary(results$mmce.test.mean)
results_flip[results_flip$mmce.test.mean == min(results_flip$mmce.test.mean), ]


# consider only accuracy and size
dat <- subset(results, select = c("mmce.test.mean", "mean_size"))
plot(dat)
dat_mco <- mco::paretoFilter(as.matrix(dat))
points(dat_mco, col = "red")

plot(dat_mco, log = "y")
