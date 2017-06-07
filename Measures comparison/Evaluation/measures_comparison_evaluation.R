####################################################
# Evaluation of the measures comparison experiments
####################################################

library(gplots)
library(ggplot2)
library(data.table)
library(GGally)

if(!file.exists("..\\Plots")) dir.create("..\\Plots")

set.seed(123)

load("results_measures_comparison.RData")
dat <- which(results$prob == "AP_Colon_Kidney")
short <- "CK"
results <- results[dat, ]

levels(results$learner) <- c("Lasso Log. Reg.", "GLM Boosting", "SVM", "Random Forest")
results$algo <- factor(results$learner, levels = c("GLM Boosting", "Lasso Log. Reg.", "Random Forest", "SVM"))

levels(results$filter) <- c("AUC", "MRMR", "Variance")

results <- cbind(results[, 1:2], fw.method = results$filter, results[, 3:(ncol(results) - 2)])



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
boxplot_order <- c(7, 6, 10, 12, 8, 2:5, 9, 11, 1)
measures_melt <- cbind(id = rownames(measures), as.data.frame(measures)[, boxplot_order])
measures_melt <- melt(measures_melt, id.vars = "id")
ggbox <- ggplot(measures_melt, aes(x = variable, y = value)) + 
  geom_boxplot() +
  labs(y = "Stability values", x = "Stability measures") +
  scale_x_discrete(expand = c(0, 0), labels = measure_names[boxplot_order]) +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 10))

pdf("..\\Plots\\measures_boxplots.pdf", width = 6.5, height = 3, useKerning = FALSE)
print(ggbox)
dev.off()  


########### scatter plot matrix ###########
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
    scale_colour_gradient(low = "darkred", high = "lightgoldenrod1", limits = c(0, 10935), 
                          name = "Mean number of chosen features",
                          guide = guide_colourbar(title.vjust = 1)) + 
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    coord_equal(ratio = 1) +
    theme(legend.position = "bottom", 
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 11),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(size = 16),
      legend.key.width = unit(2,"cm"))
  
  return(scatter_matrix)
}

pdf(file = "..\\Plots\\scatter_all.pdf", width = 20, height = 20)
print(gg_scatter_matrix(1:12000))
dev.off()

################### stability and size ###################
# correlations between measures and mean_size
apply(measures, 2, function(x) cor(x, results$mean_size, use = "pairwise.complete.obs"))

size_stability <- as.data.frame(measures)
colnames(size_stability) <- measure_names

# better order for the columns
c_order <-  c(9, 12, 3, 8, 10, 7, 4, 1, 6, 2, 5, 11)
size_stability <- size_stability[, c_order]

size_stability <- cbind(size = results$mean_size, size_stability, algo = paste(results$algo, results$fw.method, sep = " & "))

# randomize order of rows
size_stability <- size_stability[sample(1:nrow(size_stability)), ]


size_stability <- melt(size_stability, id.vars = c("size", "algo"))

ggmatrix <- ggplot(size_stability, aes(x = size, y = value)) +
  geom_point(size = 1, color = "black", alpha = 0.25) + 
  labs(x = "Mean number of chosen features", y = "Stability value") + 
  facet_wrap( ~ variable, ncol = 4, labeller = "label_parsed") + 
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        strip.text.x = element_text(size = 12)) 

pdf("..\\Plots\\size_stability.pdf", height = 6, width = 8, useKerning = FALSE)
print(ggmatrix)
dev.off()


########## optimal models ##########################

# which models are Pareto optimal when considering different measures of stability
load(paste0("results_pareto_", short, ".RData"))

# indices of methods that are pareto optimal for any stability measure
paretos <- sort(unique(unlist(lapply(results_pareto, function(x) x$pareto_indices))))
length(paretos)

# 0-1-coding: which of the above methods (row) is pareto optimal 
# for which stability measure (col)
pareto_models <- matrix(0, nrow = length(paretos), ncol = length(results_pareto))
rownames(pareto_models) <- paretos
colnames(pareto_models) <- unlist(lapply(results_pareto, function(x) x$targets[3]))
for(i in 1:length(results_pareto))
{
  pareto_models[paretos %in% results_pareto[[i]]$pareto_indices, i] <- 1
}

ggbin <- binary_heatmap_ggplot(pareto_models, measure_names, rev(c(9, 10, 6, 12, 7, 2:5, 8, 1, 11)))

pdf("..\\Plots\\heatmap_optimal.pdf", width = length(paretos) / 20, height = 3, useKerning = FALSE)
print(ggbin)
dev.off()

