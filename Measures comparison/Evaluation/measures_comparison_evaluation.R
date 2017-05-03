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

# do this analysis seperately for all classifiers and filter methods
measures_split <- split(measures, results$algo)
filter_split <- split(results$fw.method, results$algo)

for(i in 1:length(measures_split))
{
  measures_split2 <- split(measures_split[[i]], filter_split[[i]])
  
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

########## correlation depending on size ###############
sizes <- (1 + 25):(max(results$mean_size) - 25)
cors <- sapply(sizes, function(s){
  res <- results[results$mean_size >= s - 25 & results$mean_size <= s + 25, , drop  = FALSE]
  if(nrow(res) < 2) return(NA)
  cor(res$lustgarten, res$cor_pearson)
})
plot(cors, type = "l")


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

pdf(file = "..\\Plots\\scatter_glmboost.pdf", width = 20, height = 20)
print(gg_scatter_matrix(which(results$algo == "GLM Boosting"), "GLM Boosting"))
dev.off()

pdf(file = "..\\Plots\\scatter_logreg.pdf", width = 20, height = 20)
print(gg_scatter_matrix(which(results$algo == "Lasso Log. Reg."), "Lasso Logistic Regression"))
dev.off()

pdf(file = "..\\Plots\\scatter_rf.pdf", width = 20, height = 20)
print(gg_scatter_matrix(which(results$algo == "Random Forest"), "Random Forest"))
dev.off()

pdf(file = "..\\Plots\\scatter_svm.pdf", width = 20, height = 20)
print(gg_scatter_matrix(which(results$algo == "SVM"), "SVM"))
dev.off()

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
  theme(# legend.position = "right", 
        # legend.title = element_text(size = 13),
        # legend.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        strip.text.x = element_text(size = 12)) 
  # scale_color_manual(values = 
  #                      c("firebrick1", "firebrick3", "firebrick4",
  #                        "steelblue1", "steelblue3", "steelblue4", 
  #                        "darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4",
  #                        "mediumorchid1", "mediumorchid3", "mediumorchid4"),
  #                    name = "Wrapper method") + 
  # guides(colour = guide_legend(override.aes = list(size = 2)))

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

ggbin <- binary_heatmap_ggplot(pareto_models, measure_names)

pdf("..\\Plots\\heatmap_optimal.pdf", width = length(paretos) / 20, height = 3, useKerning = FALSE)
print(ggbin)
dev.off()



###################################################################

# analysis why SL behaves so oddly in size_stability
inds <- which(results$prob == "AP_Breast_Ovary" & results$mean_size >= 10500 & results$lustgarten >= 0.525)
f <- features[[inds[1]]]

res <- matrix(0, ncol = 10, nrow = 10)
for(i in 1:9)
{
  for(j in (i+1):10)
  {
    res[i, j] <- length(intersect(f[[i]], f[[j]])) / results$mean_size[inds[1]]
  }
}
# for large values of mean size, the intersection is at least 200 smaller than the equally sizes sets
# -> analyse: cardinality(V_i) = cardinality(V_j) = x and cardinality(intersection(V_i, V_j)) = x - 200
# plot for different values of x explains odd behaviour!
f1 <- function(x) 1 - 200/x - x/10000
f2 <- function(x) (x - 200 - x^2 / 10000) / (10000 - x)
curve(f1, from = 200, to = 5000, xlim = c(200, 10000), main = "SL")
curve(f2, from = 5000, to = 10000, add = TRUE)

# same analysis for SC
f3 <- function(x) cor(c(rep(1, x), rep(0, 10000 - x)), c(rep(0, 100), rep(1, x), rep(0, 10000 - 100 - x)))
f3_values <- sapply(200:9800, f3)
plot(f3_values, x = 200:9800, type = "l", main = "SC")

# same analysis for SJ
f4 <- function(x) (x - 200) / x
curve(f4, from = 200, to = 10000, main = "SJ")

