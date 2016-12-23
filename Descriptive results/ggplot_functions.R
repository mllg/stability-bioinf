####################################################
# Functions for creating the correlation heatmaps
# and the plot of Pareto optimal models
# for the descriptive experiments using ggplot2
####################################################


heatmap_ggplot <- function(cor, names, title = "")
{
  library(gplots)
  library(ggplot2)
  library(data.table)
  
  hm <- heatmap.2(cor, dendrogram = "column", distfun = function(x) as.dist(1 - x), scale = "none",
                  density.info = "none", trace = "none",
                  hclustfun = function(x) hclust(x, method = "single"))
  
  
  cluster_order <- hm$rowInd
  cor_data <- cbind(measure = rownames(cor), as.data.frame(cor)) 
  cor_data <- melt(cor_data, id.vars = "measure")
  cor_data$measure <- factor(cor_data$measure, levels = rownames(cor)[cluster_order])
  cor_data$variable <- factor(cor_data$variable, levels = rownames(cor)[cluster_order])
  
  gg_heatmap <- ggplot(cor_data, aes(measure, variable)) +
    geom_tile(aes(fill = value), colour = "white") + 
    # scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", limits = c(-1, 1), name = "Correlation") +
    scale_fill_gradient(low = "white", high = "black", limits = c(-0.25, 1), name = "Correlation\n") +
    theme_grey() + 
    labs(x = "", y = "", title = title) + 
    scale_x_discrete(expand = c(0, 0), labels = names[cluster_order]) +
    scale_y_discrete(expand = c(0, 0), labels = names[cluster_order]) + 
    theme(axis.ticks = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          title =  element_text(size = 12))
  
  return(gg_heatmap)
}

binary_heatmap_ggplot <- function(cor, names, title = "")
{
  library(gplots)
  library(ggplot2)
  library(data.table)
  
  hm <- heatmap.2(cor, dendrogram = "column", scale = "none",
                  density.info = "none", trace = "none",
                  hclustfun = function(x) hclust(x, method = "single"))
  cluster_order <- hm$colInd
  
  cor_factor <- apply(cor, 1, function(x) factor(x, levels = 0:1))
  cor_data <- cbind(measure = colnames(cor), as.data.frame(cor_factor))
  cor_data <- melt(cor_data, id.vars = "measure")
  cor_data$measure <- factor(cor_data$measure, levels = colnames(cor)[cluster_order])
  cor_data$variable <- factor(cor_data$variable, levels = rownames(cor)[hm$rowInd])
  
  gg_heatmap <- ggplot(cor_data, aes(measure, variable)) +
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_grey(name = "Pareto\noptimal", labels = c("No", "Yes"), start = 0.8, end = 0.2) +
    theme_grey() + 
    labs(x = "Stability measures", y = "Methods", title = title) + 
    scale_x_discrete(expand = c(0, 0), labels = names[cluster_order]) +
    scale_y_discrete(expand = c(0, 0), labels = element_blank()) +
    theme(axis.ticks = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
          
  
  return(gg_heatmap)
}