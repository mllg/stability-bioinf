#####################################################################
# Definition of experiments for the comparison of stability measures
#####################################################################

library(batchtools)
library(mlr)

reg <- makeExperimentRegistry(file.dir = "measures_comparison")


make_mlr_task_and_rin <- function(data, target_name)
{
  require(mlr)
  
  task <- makeClassifTask(data = data, target = target_name)
  
  # resampling
  rdesc <- makeResampleDesc(method = "CV", stratify = TRUE, iters = 10)
  rin <- makeResampleInstance(rdesc, task = task)
  
  list(task = task, rin = rin)
}


load("AP_Breast_Ovary.RData")
load("AP_Colon_Kidney.RData")
load("Stomach.RData")


addProblem(name = "AP_Breast_Ovary", seed = 1,
           data = list(mlr = make_mlr_task_and_rin(AP_Breast_Ovary, "target")
           )
)

addProblem(name = "AP_Colon_Kidney", seed = 2,
           data = list(mlr = make_mlr_task_and_rin(AP_Colon_Kidney, "target")
           )
)

addProblem(name = "Stomach", seed = 3,
           data = list(mlr = make_mlr_task_and_rin(Stomach, "target")
           )
)

############################################################

algo_classif_wrapper <- function(job, data, instance, learner, filter, ...)
{
  library(mlr)
  library(BBmisc)
  
  source("measures_of_stability.R")
  source("filter_methods.R")
  
  source("determine_parset.R")
  source("extract_functions.R")
  
  # necessary type conversions
  learner <- as.character(learner)
  filter <- as.character(filter)
  
  ps <- determine_parset(learner)
  pars <- as.list(generateDesign(1, ps, trafo = TRUE))
  
  factors <- which(unlist(lapply(pars, is.factor)))
  if(length(factors) > 0)
  {
    pars[factors] <- lapply(pars[factors], as.character)
  }
  
  lrn <- makeLearner(paste("classif", learner, sep = "."))
  lrn <- makeFilterWrapper(lrn, fw.method = filter)
  lrn <- setHyperPars(lrn, par.vals = pars)
  
  extract <- choose_extract_function(paste("classif", learner, sep = "."))
  
  r <- resample(learner = lrn, task = data$mlr$task, resampling = data$mlr$rin, show.info = FALSE,
                      extract = extract, models = FALSE)
  
  performance <- r$aggr
  features <- r$extract
  
  stability <- all_stability_measures(features, getTaskNFeats(data$mlr$task),
                                      cor(getTaskData(data$mlr$task, target.extra = TRUE)$data))
  mean_size <- mean(unlist(lapply(features, length)))
  
  results <- list(results = c(stability, performance, mean_size = mean_size), 
                  features = features,
                  parameters = pars)
  results
}

addAlgorithm(reg, name = "classif_wrapper", fun = algo_classif_wrapper)


#######################################################

design_cw <- list(classif_wrapper = expand.grid(learner = c("glmboost", "ksvm", "LiblineaRL1LogReg", "ranger"),
                                                filter = c("auc", "fmrmr", "variance"),
                                                stringsAsFactors = FALSE))

addExperiments(algo.designs = design_cw, repls = 1000)

# ids[, chunk := chunk(job.id, chunk.size = 10)]

