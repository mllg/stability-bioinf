############################################################
# Collection of the results of the
# experiments for the random search for a good model
############################################################


results <- reduceResultsDataFrame(reg, ids = findDone(reg),
                                  fun = function(job, res, ...)
                                  {
                                    c(prob = job$prob.id, algo = job$algo.id, res$parameters,
                                      res$training, res$testing)
                                  })
save(results, file = "results_stability_random.RData")




features_train <- reduceResultsList(reg, ids = findDone(reg),
                                    fun = function(job, res, ...) res$features_train)
save(features_train, file = "features_train_stability_random.RData")

features_test <- reduceResultsList(reg, ids = findDone(reg),
                                   fun = function(job, res, ...) res$features_test)
save(features_test, file = "features_test_stability_random.RData")


##############################################################################################


# calculate the values of all stability measures additionally
# (except for Zucknick due to computation time)
library(BBmisc)

load("results_stability_random.RData")
load("features_train_stability_random.RData")
load("features_test_stability_random.RData")

source("..\\Experiments\\measures_of_stability.R")
zucknick <- function(...) NA

load("..\\Data\\AP_Breast_Ovary.RData")
load("..\\Data\\AP_Colon_Kidney.RData")
load("..\\Data\\Stomach.RData")

p_BO <- ncol(AP_Breast_Ovary) - 1
p_CK <- ncol(AP_Colon_Kidney) - 1
p_St <- ncol(Stomach) - 1

names_BO <- rownames(results[results$prob == "AP_Breast_Ovary", ])
names_CK <- rownames(results[results$prob == "AP_Colon_Kidney", ])
names_St <- rownames(results[results$prob == "Stomach", ])

ids_BO_train <- which(names(features_train) %in% names_BO)
ids_BO_test <- which(names(features_test) %in% names_BO)
ids_CK_train <- which(names(features_train) %in% names_CK)
ids_CK_test <- which(names(features_test) %in% names_CK)
ids_St_train <- which(names(features_train) %in% names_St)
ids_St_test <- which(names(features_test) %in% names_St)


stability_train_BO <- lapply(features_train[ids_BO_train], function(feats){
  all_stability_measures(features = feats, p = p_BO, correlations = NA)
})

stability_train_CK <- lapply(features_train[ids_CK_train], function(feats){
  all_stability_measures(features = feats, p = p_CK, correlations = NA)
})

stability_train_St <- lapply(features_train[ids_St_train], function(feats){
  all_stability_measures(features = feats, p = p_St, correlations = NA)
})

stability_test_BO <- lapply(features_test[ids_BO_test], function(feats){
  all_stability_measures(features = feats, p = p_BO, correlations = NA)
})

stability_test_CK <-lapply(features_test[ids_CK_test], function(feats){
  all_stability_measures(features = feats, p = p_CK, correlations = NA)
})

stability_test_St <- lapply(features_test[ids_St_test], function(feats){
  all_stability_measures(features = feats, p = p_St, correlations = NA)
})



names(stability_train_BO) <- names(stability_test_BO) <- names_BO
names(stability_train_CK) <- names(stability_test_CK) <- names_CK
names(stability_train_St) <- names(stability_test_St) <- names_St


stability_train_BO <- convertListOfRowsToDataFrame(stability_train_BO)
stability_test_BO <- convertListOfRowsToDataFrame(stability_test_BO)
stability_train_CK <- convertListOfRowsToDataFrame(stability_train_CK)
stability_test_CK <- convertListOfRowsToDataFrame(stability_test_CK)
stability_train_St <- convertListOfRowsToDataFrame(stability_train_St)
stability_test_St <- convertListOfRowsToDataFrame(stability_test_St)

colnames_train <- paste0(colnames(stability_train_BO), "_train")
colnames_test <- paste0(colnames(stability_test_BO), "_test")
colnames(stability_train_BO) <- colnames(stability_train_CK) <- colnames(stability_train_St) <- colnames_train
colnames(stability_test_BO) <- colnames(stability_test_CK) <- colnames(stability_test_St) <- colnames_test

stability_BO <- cbind(stability_train_BO, stability_test_BO)
stability_CK <- cbind(stability_train_CK, stability_test_CK)
stability_St <- cbind(stability_train_St, stability_test_St)

save(stability_BO, file = "stability_BO.RData")
save(stability_CK, file = "stability_CK.RData")
save(stability_St, file = "stability_St.RData")

stability_add <- rbind(stability_BO, stability_CK, stability_St)

match_order <- match(rownames(stability_add), rownames(results))
results_all <- cbind(results[match_order, ], stability_add[match_order, ])

results_all <- results_all[, -c(27, 39)]
save(results_all, file = "results_stability_random_all.RData")



