####################################
# PCA plots of the three data sets
####################################


library(ggplot2)

load("AP_Colon_Kidney.RData")
load("AP_Breast_Ovary.RData")
load("Stomach.RData")

# microarray
data <- rbind(AP_Breast_Ovary, AP_Colon_Kidney)

m <- scale(data[, -1])
pca <- prcomp(m, scale = FALSE)

# explained variation
sum(pca$sdev[1:2]^2)/sum(pca$sdev^2)
exp <- pca$sdev[1:2]^2 / sum(pca$sdev^2)

col_pca <- c("red", "blue", "goldenrod1", "black")[as.numeric(data$target)]
pcs <- m %*% pca$rotation

pcs_part <- cbind(as.data.frame(pcs[, 1:2]), Class = data$target)
pcgg <- ggplot(pcs_part, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 1) + 
  coord_equal(ratio = 1) +
  labs(x = paste0("PC1 (", round(exp[1] * 100, 2), "% of data variation)"),
       y = paste0("PC2 (", round(exp[2] * 100, 2), "% of data variation)"))

pdf("..\\Plots\\pca_ap.pdf", height = 5, width = 6.5, useKerning = FALSE)
print(pcgg)
dev.off()


# RNA seq
data <- Stomach

m <- scale(data[, -1])
pca <- prcomp(m, scale = FALSE)

# explained variation
sum(pca$sdev[1:2]^2)/sum(pca$sdev^2)
exp <- pca$sdev[1:2]^2 / sum(pca$sdev^2)

col_pca <- c("red", "blue")[as.numeric(data$target)]
pcs <- scale(m, scale = FALSE) %*% pca$rotation

pcs_part <- cbind(as.data.frame(pcs[, 1:2]), Class = data$target)
pcgg <- ggplot(pcs_part, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 1) + 
  coord_equal(ratio = 1) +
  labs(x = paste0("PC1 (", round(exp[1] * 100, 2), "% of data variation)"),
       y = paste0("PC2 (", round(exp[2] * 100, 2), "% of data variation)"))

pdf("..\\Plots\\pca_stomach.pdf", height = 5, width = 7, useKerning = FALSE)
print(pcgg)
dev.off()


