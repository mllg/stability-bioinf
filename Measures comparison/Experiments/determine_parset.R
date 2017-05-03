########################################################
# Determine the hyper parameters to be optimized 
# for the given learner
########################################################

determine_parset <- function(learner)
{
  library(checkmate)
  assert_choice(learner, choices = c("glmboost", "ksvm", "LiblineaRL1LogReg", "ranger"))
  
  if(learner == "glmboost")
  {
    ps <- makeParamSet(
      makeIntegerParam("mstop", lower = 1, upper = 2^15),
      makeNumericParam("fw.perc", lower = 0, upper = 1)
    )
  }
  else if(learner == "ksvm")
  {
    ps <- makeParamSet(
      makeNumericParam("C", lower = -15, upper = 15, trafo = function(x) 2^x),
      makeNumericParam("sigma", lower = -15, upper = 15, trafo = function(x) 2^x),
      makeNumericParam("fw.perc", lower = 0, upper = 1)
    ) 
  }
  else if(learner == "LiblineaRL1LogReg")
  {
    ps <- makeParamSet(
      makeNumericParam("cost", lower = -15, upper = 15, trafo = function(x) 1/(2^x)),
      makeNumericParam("fw.perc", lower = 0, upper = 1)
    )
  }
  else
  {
    ps <- makeParamSet(
      makeIntegerParam("num.trees", lower = 1, upper = 2^15),
      makeIntegerParam("min.node.size", lower = 1, upper = 2^5),
      makeDiscreteParam("importance", values = "impurity"),
      makeNumericParam("fw.perc", lower = 0, upper = 1)
    )
  }
  
  ps
}