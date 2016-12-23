############################################################
# Definition of experiments for the descriptive comparision
# of the stability measures
############################################################


library(BatchExperiments)
reg <- makeExperimentRegistry(id = "stability_descriptive")

load("AP_Colon_Kidney.RData")
load("AP_Breast_Ovary.RData")
load("Stomach.RData")


addProblem(reg, id = "AP_Colon_Kidney", seed = 1,
           static = list(data = AP_Colon_Kidney),
           dynamic = function(static, cv_folds = 10, ...)
           {
             library(mlr)
             task <- makeClassifTask(data = static$data, target = "target")
             rdesc <- makeResampleDesc(method = "CV", stratify = TRUE, iters = cv_folds)
             rin <- makeResampleInstance(rdesc, task = task)
             return(list(task = task, rdesc = rdesc, rin = rin))
           }
)

addProblem(reg, id = "AP_Breast_Ovary", seed = 2,
           static = list(data = AP_Breast_Ovary),
           dynamic = function(static, cv_folds = 10, ...)
           {
             library(mlr)
             task <- makeClassifTask(data = static$data, target = "target")
             rdesc <- makeResampleDesc(method = "CV", stratify = TRUE, iters = cv_folds)
             rin <- makeResampleInstance(rdesc, task = task)
             return(list(task = task, rdesc = rdesc, rin = rin))
           }
)


algo_logreg <- function(static, dynamic, lambda, fw.method, fw.abs, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")

  lrn <- makeLearner("classif.LiblineaRL1LogReg", par.vals = list(cost = 1/lambda))
  lrn <- makeFilterWrapper(lrn, fw.method = fw.method, fw.abs = fw.abs)
  
  r <- resample(learner = lrn, task = dynamic$task, resampling = dynamic$rin, show.info = FALSE,
                extract = function(x) drop(x$learner.model$next.model$learner.model$W), models = FALSE)
  performance <- r$aggr
  features <- lapply(r$extract, function(x) names(x)[which(abs(x[-length(x)]) > 1e-16)])
  
  stability <- all_stability_measures(features, ncol(static$data) - 1,
                                      cor(subset(static$data, select = -get("target"))))
  
  time <- mean(r$pred$time, na.rm = TRUE)
  names(time) <- "time"
  
  mean_size <- mean(unlist(lapply(features, length)))
  names(mean_size) <- "mean_size"
  
  results <- list(results = c(stability, performance, time, mean_size), features = features)
  return(results)
}

algo_glmboost <- function(static, dynamic, mstop, fw.method, fw.abs, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  lrn <- makeLearner("classif.glmboost", par.vals = list(mstop = mstop))
  lrn <- makeFilterWrapper(lrn, fw.method = fw.method, fw.abs = fw.abs)
  
  r <- resample(learner = lrn, task = dynamic$task, resampling = dynamic$rin, show.info = FALSE,
                extract = function(x) unlist(x$learner.model$next.model$learner.model$coef()),
                models = FALSE)
  performance <- r$aggr
  features <- lapply(r$extract, function(x) setdiff(names(x), "(Intercept)"))
  
  stability <- all_stability_measures(features, ncol(static$data) - 1,
                                      cor(subset(static$data, select = -get("target"))))
  
  time <- mean(r$pred$time, na.rm = TRUE)
  names(time) <- "time"
  
  mean_size <- mean(unlist(lapply(features, length)))
  names(mean_size) <- "mean_size"
  
  results <- list(results = c(stability, performance, time, mean_size), features = features)
  return(results)
}

algo_svm <- function(static, dynamic, C, sigma, fw.method, fw.abs, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  lrn <- makeLearner("classif.ksvm", par.vals = list(C = C, sigma = sigma))
  lrn <- makeFilterWrapper(lrn, fw.method = fw.method, fw.abs = fw.abs)

  r <- resample(learner = lrn, task = dynamic$task, resampling = dynamic$rin, show.info = FALSE,
                models = TRUE)
  performance <- r$aggr
  features <- lapply(r$models, getFilteredFeatures)
  
  stability <- all_stability_measures(features, ncol(static$data) - 1,
                                      cor(subset(static$data, select = -get("target"))))
  
  time <- mean(r$pred$time, na.rm = TRUE)
  names(time) <- "time"
  
  mean_size <- mean(unlist(lapply(features, length)))
  names(mean_size) <- "mean_size"
  
  results <- list(results = c(stability, performance, time, mean_size), features = features)
  return(results)
}

algo_random_forest <- function(static, dynamic, num.trees, min.node.size, fw.method, fw.abs, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  lrn <- makeLearner("classif.ranger", par.vals = list(num.trees = num.trees, min.node.size = min.node.size,
                                                       importance = "impurity"))
  lrn <- makeFilterWrapper(lrn, fw.method = fw.method, fw.abs = fw.abs)

  r <- resample(learner = lrn, task = dynamic$task, resampling = dynamic$rin, show.info = FALSE,
                extract = function(x) drop(x$learner.model$next.model$learner.model$variable.importance),
                models = FALSE)
  performance <- r$aggr
  features <- lapply(r$extract, function(x) names(x)[which(x > 1e-16)])
  
  stability <- all_stability_measures(features, ncol(static$data) - 1,
                                      cor(subset(static$data, select = -get("target"))))
  
  time <- mean(r$pred$time, na.rm = TRUE)
  names(time) <- "time"
  
  mean_size <- mean(unlist(lapply(features, length)))
  names(mean_size) <- "mean_size"
  
  results <- list(results = c(stability, performance, time, mean_size), features = features)
  return(results)
}





par_fw.abs <- c(1, 2, 5, 10, 15, 20, 50, 100, 200, 500, 1000, 2000, 5000, ncol(AP_Colon_Kidney) - 1)

design_logreg <- makeDesign(id = "logreg",
                            design = expand.grid(lambda = 2^((-10):10),
                                                 fw.method = c("variance", "auc", "fmrmr"),
                                                 fw.abs = par_fw.abs,
                                                 stringsAsFactors = FALSE))

design_glmboost <- makeDesign(id = "glmboost",
                              design = expand.grid(mstop = c(50, 100, 250, 500),
                                                   fw.method = c("variance", "auc", "fmrmr"),
                                                   fw.abs = par_fw.abs,
                                                   stringsAsFactors = FALSE))

design_svm <- makeDesign(id = "svm",
                         design = expand.grid(C = 2^(-10:10),
                                              sigma = 2^(-10:10),
                                              fw.method = c("variance", "auc", "fmrmr"),
                                              fw.abs = par_fw.abs,
                                              stringsAsFactors = FALSE))

design_random_forest <- makeDesign(id = "random_forest",
                                   design = expand.grid(fw.method = c("variance", "auc", "fmrmr"),
                                                        fw.abs = par_fw.abs,
                                                        num.trees = c(250, 500, 1000),
                                                        min.node.size = 1:10,
                                                        stringsAsFactors = FALSE))


addAlgorithm(reg, id = "logreg", fun = algo_logreg)
addAlgorithm(reg, id = "glmboost", fun = algo_glmboost)
addAlgorithm(reg, id = "svm", fun = algo_svm)
addAlgorithm(reg, id = "random_forest", fun = algo_random_forest)



addExperiments(reg, algo.designs = list(design_logreg, design_glmboost, design_svm, design_random_forest))

chunked <- chunk(getJobIds(reg), n.chunks = 500, shuffle = FALSE)
submitJobs(reg, chunked, resources = list(memory = 16*1024), job.delay = TRUE)



# add Stomach

addProblem(reg, id = "Stomach", seed = 3,
           static = list(data = Stomach),
           dynamic = function(static, cv_folds = 10, ...)
           {
             library(mlr)
             task <- makeClassifTask(data = static$data, target = "target")
             rdesc <- makeResampleDesc(method = "CV", stratify = TRUE, iters = cv_folds)
             rin <- makeResampleInstance(rdesc, task = task)
             return(list(task = task, rdesc = rdesc, rin = rin))
           }
)

par_fw.abs2 <- c(1, 2, 5, 10, 15, 20, 50, 100, 200, 500, 1000, 2000, 5000, ncol(Stomach) - 1)

design_logreg2 <- makeDesign(id = "logreg",
                            design = expand.grid(lambda = 2^((-10):10),
                                                 fw.method = c("variance", "auc", "fmrmr"),
                                                 fw.abs = par_fw.abs2,
                                                 stringsAsFactors = FALSE))

design_glmboost2 <- makeDesign(id = "glmboost",
                              design = expand.grid(mstop = c(50, 100, 250, 500),
                                                   fw.method = c("variance", "auc", "fmrmr"),
                                                   fw.abs = par_fw.abs2,
                                                   stringsAsFactors = FALSE))

design_svm2 <- makeDesign(id = "svm",
                         design = expand.grid(C = 2^(-10:10),
                                              sigma = 2^(-10:10),
                                              fw.method = c("variance", "auc", "fmrmr"),
                                              fw.abs = par_fw.abs2,
                                              stringsAsFactors = FALSE))

design_random_forest2 <- makeDesign(id = "random_forest",
                                   design = expand.grid(fw.method = c("variance", "auc", "fmrmr"),
                                                        fw.abs = par_fw.abs2,
                                                        num.trees = c(250, 500, 1000),
                                                        min.node.size = 1:10,
                                                        stringsAsFactors = FALSE))

addExperiments(reg, prob.designs = "Stomach",
               algo.designs = list(design_logreg2, design_glmboost2, design_svm2, design_random_forest2))
