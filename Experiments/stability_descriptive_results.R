############################################################
# Collection of the results of the
# experiments for the descriptive comparision
# of the stability measures
############################################################


results <- reduceResultsDataFrame(reg, ids = findDone(reg),
                                  fun = function(job, res, ...)
                                    {
                                    c(prob = job$prob.id, algo = job$algo.id,
                                      job$algo.pars, res$results)
                                    })
save(results, file = "results_stability_descriptive.RData")

features <- reduceResultsList(reg, ids = findDone(reg),
                              fun = function(job, res, ...) res$features)
save(features, file = "features_stability_descriptive.RData")
