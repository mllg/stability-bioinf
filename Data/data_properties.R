####################################
# PCA plots of the three data sets
####################################


library(ggplot2)

load("AP_Colon_Kidney.RData")
load("AP_Breast_Ovary.RData")
load("Stomach.RData")

if(!file.exists("Plots")) dir.create("Plots")


plot_pca <- function(data, pca, m)
{
  exp <- pca$sdev[1:2]^2 / sum(pca$sdev^2)
  pcs <- m %*% pca$rotation
  
  pcs_part <- cbind(as.data.frame(pcs[, 1:2]), Class = data$target)
  
  pcgg <- ggplot(pcs_part, aes(x = PC1, y = PC2, color = Class, shape = Class)) +
    geom_point(size = 2) + 
    coord_equal(ratio = 1) +
    labs(x = paste0("PC1 (", round(exp[1] * 100, 2), "% of data variation)"),
         y = paste0("PC2 (", round(exp[2] * 100, 2), "% of data variation)")) +
    theme(legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.position = "bottom")
  pcgg
}

# AP_Breast_Ovary
m_BO <- scale(AP_Breast_Ovary[, -1])
pca_BO <- prcomp(m_BO, scale = FALSE)

# explained variation
sum(pca_BO$sdev[1:2]^2)/sum(pca_BO$sdev^2)

pdf("Plots\\pca_BO.pdf", height = 5, width = 5, useKerning = FALSE)
print(plot_pca(AP_Breast_Ovary, pca_BO, m_BO))
dev.off()


# AP_Colon_Kidney
m_CK <- scale(AP_Colon_Kidney[, -1])
pca_CK <- prcomp(m_CK, scale = FALSE)

# explained variation
sum(pca_CK$sdev[1:2]^2)/sum(pca_CK$sdev^2)

pdf("Plots\\pca_CK.pdf", height = 5, width = 5, useKerning = FALSE)
print(plot_pca(AP_Colon_Kidney, pca_CK, m_CK))
dev.off()



# Stomach
m_St <- scale(Stomach[, -1])
pca_St <- prcomp(m_St, scale = FALSE)

# explained variation
sum(pca_St$sdev[1:2]^2)/sum(pca_St$sdev^2)

pdf("Plots\\pca_St.pdf", height = 5, width = 5, useKerning = FALSE)
print(plot_pca(Stomach, pca_St, m_St))
dev.off()


