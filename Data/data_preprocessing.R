#######################################
# Preprocessing of the three data sets
#######################################
library(OpenML)

preprocess_data <- function(data, target, ignore = character(0))
{
  # remove columns like id columns
  if(length(ignore) > 0)
  {
    remove <- which(colnames(data) %in% ignore)
    data <- subset(data, select = -remove)
  }
  
  # target
  target_index <- which(names(data) == target)
  
  # remove constant columns
  useful <- which(apply(data[, -target_index], 2, function(x) min(x) < max(x)))
  data <- cbind(as.factor(data[, target_index]), data[, -target_index][, useful])
  
  colnames(data) <- c("target", paste0("V", 1:length(useful)))
  return(data)
}

AP_Breast_Ovary <- getOMLDataSet(1165)$data
AP_Breast_Ovary <- preprocess_data(AP_Breast_Ovary, "Tissue", "ID_REF")
save(AP_Breast_Ovary, file = "AP_Breast_Ovary.RData")

AP_Colon_Kidney <- getOMLDataSet(1137)$data
AP_Colon_Kidney <- preprocess_data(AP_Colon_Kidney, "Tissue", "ID_REF")
save(AP_Colon_Kidney, file = "AP_Colon_Kidney.RData")


# NOTE: You need to download these data sets from the referenced paper and place them in 
# this directory
temp <- read.table("RPKM_Expression_Matrix.291Samples_GAF3genes.BCGSC.20131127.tsv", header = TRUE)
stomach_info <- xlsx::read.xlsx("20140110_STAD_Clinical_Data_Blacklisted_Cases_Removed.xlsx", 1, startRow = 2)

stomach <- as.data.frame(t(as.matrix(temp[, -1])))
stomach <- cbind(Patient = rownames(stomach), stomach, stringsAsFactors = FALSE)
#colnames(stomach)[-1] <- temp[, 1]

# patient information as in clinical file
stomach$Patient <- gsub(".", "-", stomach$Patient, fixed = TRUE)
stomach$Patient <- substr(stomach$Patient, 1, 12)

# merge gene expression and patient information
stomach_all <- merge(x = stomach, y = stomach_info, by.x = "Patient", by.y = "bcr_patient_barcode")

# transformation log(x+1)
stomach_all[, 2:29700] <- apply(stomach_all[, 2:29700], 2, function(x) log(x + 1))

# consider only types C16.1 (110) and C16.3 (105)
stomach_two <- stomach_all[stomach_all$icd_10 %in% c("C16.1", "C16.3"),]

target <- droplevels(stomach_two$icd_10)
Stomach <- cbind(target = target, stomach_two[, 2:29700])

# filter: 10 000 variables with highest variance
vars <- apply(Stomach[, -1], 2, var)
vars_order <- order(vars, decreasing = TRUE)[1:10000]

# keep order
greatest_vars <- sort(vars_order)
Stomach <- Stomach[, c(1, 1 + greatest_vars)]

save(Stomach, file = "Stomach.RData")
