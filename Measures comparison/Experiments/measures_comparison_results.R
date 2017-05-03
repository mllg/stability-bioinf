##################################################################
# Collection of the results of the experiments for the
# comparison of stability measures
##################################################################

options(batchtools.progress = FALSE)

results <- reduceResultsList(ids = findDone(),
                                  fun = function(job, res, ...)
                                    {
                                    c(prob = job$prob.name, algo = job$algo.name,
                                      res$results, job$pars$algo.pars, res$parameters)
                                    })
results <- BBmisc::convertListOfRowsToDataFrame(results)

save(results, file = "results_measures_comparison.RData")

features <- reduceResultsList(ids = findDone(),
                              fun = function(job, res, ...) res$features)
save(features, file = "features_measures_comparison.RData")

