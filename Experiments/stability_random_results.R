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
