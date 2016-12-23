############################################################
# Definition of experiments for the random search for a
# good model
############################################################


library(BatchExperiments)
reg <- makeExperimentRegistry(id = "stability_random")


# load data sets
load("AP_Colon_Kidney.RData")
load("AP_Breast_Ovary.RData")
load("Stomach.RData")


dynamic <- function(static, ...)
{
  library(mlr)
  
  # data for training
  n <- nrow(static$data)
  
  # stratification
  class1_indices <- which(static$data$target == static$data$target[1])
  class2_indices <- setdiff(1:n, class1_indices)
  
  n1 <- length(class1_indices)
  n2 <- length(class2_indices)
  
  sample1 <- sample(class1_indices, ceiling(n1/2))
  sample2 <- sample(class2_indices, ceiling(n2/2))
  sample <- c(sample1, sample2)
  
  train <- static$data[sample, ]
  task_train <- makeClassifTask(data = train, target = "target")
  
  # data for testing
  test <- static$data[setdiff(1:n, sample), ]
  task_test <- makeClassifTask(data = test, target = "target")

  # resampling
  rdesc <- makeResampleDesc(method = "CV", stratify = TRUE, iters = 10)
  rin_train <- makeResampleInstance(rdesc, task = task_train)
  rin_test <- makeResampleInstance(rdesc, task = task_test)
  
  return(list(task_train = task_train, task_test = task_test,
              rin_train = rin_train, rin_test = rin_test))
}


addProblem(reg, id = "AP_Colon_Kidney", seed = 1,
           static = list(data = AP_Colon_Kidney),
           dynamic = dynamic)

addProblem(reg, id = "AP_Breast_Ovary", seed = 2,
           static = list(data = AP_Breast_Ovary),
           dynamic = dynamic)

addProblem(reg, id = "Stomach", seed = 3,
           static = list(data = Stomach),
           dynamic = dynamic)



algo_logreg <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("lambda", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1),
    makeDiscreteParam("fw.method", values = c("variance", "auc", "fmrmr"))
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.LiblineaRL1LogReg", par.vals = list(cost = 1 / pars$lambda))
  lrn <- makeFilterWrapper(lrn, fw.method = as.character(pars$fw.method), fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                     extract = function(x) drop(x$learner.model$next.model$learner.model$W), models = FALSE)
  
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) names(x)[which(abs(x[-length(x)]) > 1e-16)])
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))

  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) drop(x$learner.model$next.model$learner.model$W), models = FALSE)
  
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) names(x)[which(abs(x[-length(x)]) > 1e-16)])
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}

algo_glmboost <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("mstop", lower = 0, upper = 15, trafo = function(x) floor(2^x)),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1),
    makeDiscreteParam("fw.method", values = c("variance", "auc", "fmrmr"))
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.glmboost", par.vals = list(mstop = pars$mstop))
  lrn <- makeFilterWrapper(lrn, fw.method = as.character(pars$fw.method), fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      extract = function(x) unlist(x$learner.model$next.model$learner.model$coef()),
                      models = FALSE)
  
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) setdiff(names(x), "(Intercept)"))
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))

  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) unlist(x$learner.model$next.model$learner.model$coef()),
                     models = FALSE)
  
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) setdiff(names(x), "(Intercept)"))
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}

algo_svm <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("C", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("sigma", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1),
    makeDiscreteParam("fw.method", values = c("variance", "auc", "fmrmr"))
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.ksvm", par.vals = list(C = pars$C, sigma = pars$sigma))
  lrn <- makeFilterWrapper(lrn, fw.method = as.character(pars$fw.method), fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      models = TRUE)
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$models, getFilteredFeatures)
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))

  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     models = TRUE)
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$models, getFilteredFeatures)
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))

  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}

algo_random_forest <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("num.trees", lower = 0, upper = 15, trafo = function(x) floor(2^x)),
    makeNumericParam("min.node.size", lower = 0, upper = 5, trafo = function(x) floor(2^x)),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1),
    makeDiscreteParam("fw.method", values = c("variance", "auc", "fmrmr"))
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.ranger", par.vals = list(num.trees = pars$num.trees, min.node.size = pars$min.node.size,
                                                       importance = "impurity"))
  lrn <- makeFilterWrapper(lrn, fw.method = as.character(pars$fw.method), fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      extract = function(x) drop(x$learner.model$next.model$learner.model$variable.importance),
                      models = FALSE)
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) names(x)[which(x > 1e-16)])
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))

  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) drop(x$learner.model$next.model$learner.model$variable.importance),
                     models = FALSE)
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) names(x)[which(abs(x) > 1e-16)])
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))

  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}


addAlgorithm(reg, id = "logreg", fun = algo_logreg)
addAlgorithm(reg, id = "glmboost", fun = algo_glmboost)
addAlgorithm(reg, id = "svm", fun = algo_svm)
addAlgorithm(reg, id = "random_forest", fun = algo_random_forest)

##########################################################

algo_logreg_auc <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("lambda", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.LiblineaRL1LogReg", par.vals = list(cost = 1 / pars$lambda))
  lrn <- makeFilterWrapper(lrn, fw.method = "auc", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      extract = function(x) drop(x$learner.model$next.model$learner.model$W), models = FALSE)
  
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) names(x)[which(abs(x[-length(x)]) > 1e-16)])
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) drop(x$learner.model$next.model$learner.model$W), models = FALSE)
  
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) names(x)[which(abs(x[-length(x)]) > 1e-16)])
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}

algo_logreg_fmrmr <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("lambda", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.LiblineaRL1LogReg", par.vals = list(cost = 1 / pars$lambda))
  lrn <- makeFilterWrapper(lrn, fw.method = "fmrmr", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      extract = function(x) drop(x$learner.model$next.model$learner.model$W), models = FALSE)
  
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) names(x)[which(abs(x[-length(x)]) > 1e-16)])
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) drop(x$learner.model$next.model$learner.model$W), models = FALSE)
  
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) names(x)[which(abs(x[-length(x)]) > 1e-16)])
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}

algo_logreg_variance <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("lambda", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.LiblineaRL1LogReg", par.vals = list(cost = 1 / pars$lambda))
  lrn <- makeFilterWrapper(lrn, fw.method = "variance", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      extract = function(x) drop(x$learner.model$next.model$learner.model$W), models = FALSE)
  
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) names(x)[which(abs(x[-length(x)]) > 1e-16)])
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) drop(x$learner.model$next.model$learner.model$W), models = FALSE)
  
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) names(x)[which(abs(x[-length(x)]) > 1e-16)])
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}


algo_glmboost_auc <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("mstop", lower = 0, upper = 15, trafo = function(x) floor(2^x)),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.glmboost", par.vals = list(mstop = pars$mstop))
  lrn <- makeFilterWrapper(lrn, fw.method = "auc", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      extract = function(x) unlist(x$learner.model$next.model$learner.model$coef()),
                      models = FALSE)
  
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) setdiff(names(x), "(Intercept)"))
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) unlist(x$learner.model$next.model$learner.model$coef()),
                     models = FALSE)
  
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) setdiff(names(x), "(Intercept)"))
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}

algo_glmboost_fmrmr <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("mstop", lower = 0, upper = 15, trafo = function(x) floor(2^x)),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.glmboost", par.vals = list(mstop = pars$mstop))
  lrn <- makeFilterWrapper(lrn, fw.method = "fmrmr", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      extract = function(x) unlist(x$learner.model$next.model$learner.model$coef()),
                      models = FALSE)
  
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) setdiff(names(x), "(Intercept)"))
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) unlist(x$learner.model$next.model$learner.model$coef()),
                     models = FALSE)
  
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) setdiff(names(x), "(Intercept)"))
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}

algo_glmboost_variance <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("mstop", lower = 0, upper = 15, trafo = function(x) floor(2^x)),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.glmboost", par.vals = list(mstop = pars$mstop))
  lrn <- makeFilterWrapper(lrn, fw.method = "variance", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      extract = function(x) unlist(x$learner.model$next.model$learner.model$coef()),
                      models = FALSE)
  
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) setdiff(names(x), "(Intercept)"))
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) unlist(x$learner.model$next.model$learner.model$coef()),
                     models = FALSE)
  
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) setdiff(names(x), "(Intercept)"))
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}


algo_svm_auc <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("C", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("sigma", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.ksvm", par.vals = list(C = pars$C, sigma = pars$sigma))
  lrn <- makeFilterWrapper(lrn, fw.method = "auc", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      models = TRUE)
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$models, getFilteredFeatures)
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     models = TRUE)
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$models, getFilteredFeatures)
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}

algo_svm_fmrmr <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("C", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("sigma", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.ksvm", par.vals = list(C = pars$C, sigma = pars$sigma))
  lrn <- makeFilterWrapper(lrn, fw.method = "fmrmr", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      models = TRUE)
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$models, getFilteredFeatures)
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     models = TRUE)
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$models, getFilteredFeatures)
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}

algo_svm_variance <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("C", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("sigma", lower = -15, upper = 15, trafo = function(x) 2^x),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.ksvm", par.vals = list(C = pars$C, sigma = pars$sigma))
  lrn <- makeFilterWrapper(lrn, fw.method = "variance", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      models = TRUE)
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$models, getFilteredFeatures)
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     models = TRUE)
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$models, getFilteredFeatures)
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}


algo_random_forest_auc <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("num.trees", lower = 0, upper = 15, trafo = function(x) floor(2^x)),
    makeNumericParam("min.node.size", lower = 0, upper = 5, trafo = function(x) floor(2^x)),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.ranger", par.vals = list(num.trees = pars$num.trees, min.node.size = pars$min.node.size,
                                                       importance = "impurity"))
  lrn <- makeFilterWrapper(lrn, fw.method = "auc", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      extract = function(x) drop(x$learner.model$next.model$learner.model$variable.importance),
                      models = FALSE)
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) names(x)[which(x > 1e-16)])
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) drop(x$learner.model$next.model$learner.model$variable.importance),
                     models = FALSE)
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) names(x)[which(abs(x) > 1e-16)])
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}

algo_random_forest_fmrmr <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("num.trees", lower = 0, upper = 15, trafo = function(x) floor(2^x)),
    makeNumericParam("min.node.size", lower = 0, upper = 5, trafo = function(x) floor(2^x)),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.ranger", par.vals = list(num.trees = pars$num.trees, min.node.size = pars$min.node.size,
                                                       importance = "impurity"))
  lrn <- makeFilterWrapper(lrn, fw.method = "fmrmr", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      extract = function(x) drop(x$learner.model$next.model$learner.model$variable.importance),
                      models = FALSE)
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) names(x)[which(x > 1e-16)])
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) drop(x$learner.model$next.model$learner.model$variable.importance),
                     models = FALSE)
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) names(x)[which(abs(x) > 1e-16)])
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}

algo_random_forest_variance <- function(static, dynamic, ...)
{
  library(mlr)
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  ps <- makeParamSet(
    makeNumericParam("num.trees", lower = 0, upper = 15, trafo = function(x) floor(2^x)),
    makeNumericParam("min.node.size", lower = 0, upper = 5, trafo = function(x) floor(2^x)),
    makeNumericParam("fw.perc", lower = 1/10000, upper = 1)
  )
  
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  lrn <- makeLearner("classif.ranger", par.vals = list(num.trees = pars$num.trees, min.node.size = pars$min.node.size,
                                                       importance = "impurity"))
  lrn <- makeFilterWrapper(lrn, fw.method = "variance", fw.perc = pars$fw.perc)
  
  # resampling on train data set
  r_train <- resample(learner = lrn, task = dynamic$task_train, resampling = dynamic$rin_train, show.info = FALSE,
                      extract = function(x) drop(x$learner.model$next.model$learner.model$variable.importance),
                      models = FALSE)
  performance_train <- r_train$aggr
  features_train <- lapply(r_train$extract, function(x) names(x)[which(x > 1e-16)])
  stability_train <- jaccard(features_train)
  mean_size_train <- mean(unlist(lapply(features_train, length)))
  
  # resampling on test data set
  r_test <- resample(learner = lrn, task = dynamic$task_test, resampling = dynamic$rin_test, show.info = FALSE,
                     extract = function(x) drop(x$learner.model$next.model$learner.model$variable.importance),
                     models = FALSE)
  performance_test <- r_test$aggr
  features_test <- lapply(r_test$extract, function(x) names(x)[which(abs(x) > 1e-16)])
  stability_test <- jaccard(features_test)
  mean_size_test <- mean(unlist(lapply(features_test, length)))
  
  
  results <- list(training = c(stability_train = stability_train, performance_train = performance_train, mean_size_train = mean_size_train), 
                  testing = c(stability_test = stability_test, performance_test = performance_test, mean_size_test = mean_size_test),
                  parameters = pars, 
                  features_train = features_train, features_test = features_test)
  return(results)
}


addAlgorithm(reg, id = "logreg_auc", fun = algo_logreg_auc)
addAlgorithm(reg, id = "logreg_fmrmr", fun = algo_logreg_fmrmr)
addAlgorithm(reg, id = "logreg_variance", fun = algo_logreg_variance)

addAlgorithm(reg, id = "glmboost_auc", fun = algo_glmboost_auc)
addAlgorithm(reg, id = "glmboost_fmrmr", fun = algo_glmboost_fmrmr)
addAlgorithm(reg, id = "glmboost_variance", fun = algo_glmboost_variance)

addAlgorithm(reg, id = "svm_auc", fun = algo_svm_auc)
addAlgorithm(reg, id = "svm_fmrmr", fun = algo_svm_fmrmr)
addAlgorithm(reg, id = "svm_variance", fun = algo_svm_variance)

addAlgorithm(reg, id = "random_forest_auc", fun = algo_random_forest_auc)
addAlgorithm(reg, id = "random_forest_fmrmr", fun = algo_random_forest_fmrmr)
addAlgorithm(reg, id = "random_forest_variance", fun = algo_random_forest_variance)

#######################################################

addExperiments(reg, repls = 1000, skip.defined = TRUE)
submitJobs(reg, getJobIds(reg), resources = list(memory = 16*1024), job.delay = TRUE)


