################################################
# Functions for extracting the chosen features
# by the different learners
################################################


# glmboost
extract_glmboost <- function(x) 
{
  if(length(getFilteredFeatures(x)) == 0) return(character(0))
  coef <- unlist(x$learner.model$next.model$learner.model$coef() )
  setdiff(names(coef), "(Intercept)")
}

# ksvm
extract_ksvm <- function(x) getFilteredFeatures(x)

# LiblineaR
extract_LiblineaR <- function(x)
{
  w <- drop(x$learner.model$next.model$learner.model$W)
  if(is.null(w)) return(character(0))
  w <- w[-length(w)] # remove Intercept
  names(w)[which(abs(w) > 1e-16)]
}

# ranger
extract_ranger <- function(x)
{
  v <- drop(x$learner.model$next.model$learner.model$variable.importance)
  if(is.null(v)) return(character(0))
  names(v)[which(v > 1e-16)]     
}



# function to choose the right extract function
choose_extract_function <- function(learner)
{
  require(checkmate)
  assertChoice(learner, choices = c("classif.glmboost",
                                    "classif.ksvm",
                                    "classif.LiblineaRL1L2SVC",
                                    "classif.LiblineaRL2L1SVC",
                                    "classif.LiblineaRL2SVC",
                                    "classif.LiblineaRL1LogReg",
                                    "classif.LiblineaRL2LogReg",
                                    "classif.ranger"))
  
  # omit "classif."
  learner <- substring(learner, first = 9)
  
  if(grepl(pattern = "LiblineaR", x = learner, fixed = TRUE)) extract <- extract_LiblineaR
  else extract <- get(paste0("extract_", learner))
  
  extract
}