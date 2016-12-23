###############################################
# Filter methods for wrapper methods:
# additional filters that are not part of mlr
###############################################


makeFilter(
  name = "auc",
  desc = "Simple AUC filter for numeric features",
  pkg  = "ROCR",
  supported.tasks = "classif",
  supported.features = "numerics",
  fun = function(task, nselect, ...) {
    data = getTaskData(task, target.extra = TRUE)
    score = vnapply(data$data, function(x, y) {
      pred = ROCR::prediction(predictions = x, labels = y)
      ROCR::performance(pred, "auc")@y.values[[1L]]
    }, y = data$target)
    abs(0.5 - score)
  }
)


makeFilter(
  name = "mrmr.classif",
  desc = "MRMR filter for binary classification",
  pkg  = "sideChannelAttack",
  supported.tasks = "classif",
  supported.features = "numerics",
  fun = function(task, nselect, ...) {
    data = getTaskData(task, target.extra = TRUE)
    filter.order = filter.mRMR(X = as.matrix(data$data),
                               Y = as.vector(data$target),
                               nbreVarX_ = nselect)$filter
    score.values = (nselect:1) / nselect
    scores = rep(NA, ncol(data$data))
    scores[filter.order] = score.values
    scores
  }
)


makeFilter(
  name = "fmrmr",
  desc = "MRMR filter for binary classification",
  pkg  = c("fmrmr", "ROCR"),
  supported.tasks = "classif",
  supported.features = "numerics",
  fun = function(task, nselect, ...) {
    data = getTaskData(task, target.extra = TRUE)
    
    # AUC scores
    score = vnapply(data$data, function(x, y) {
      pred = ROCR::prediction(predictions = x, labels = y)
      ROCR::performance(pred, "auc")@y.values[[1L]]
    }, y = data$target)
    auc.scores = abs(0.5 - score)
    
    # MRMR
    x = as.matrix(data$data)
    fmrmr::mrmr_with_relevance(auc.scores, x, nselect = nselect)
  }
)

