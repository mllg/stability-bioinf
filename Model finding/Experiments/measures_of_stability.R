######################################
# Measures for stability assessment
######################################



require(checkmate)

# input convention:
# features: (list of character vectors) each vector contains the names of the chosen features
#           all vectors must have strictly positive length
# p: (integerish(1)) number of variables in dataset, must be strictly positive
# penalty: (numeric(1)) penalty term, must be positive
# correlations: (numeric matrix(pxp)) of correlations

davis <- function(features, p, penalty = 0, ...)
{
  # make sure that input has desired form
  assertList(x = features, types = "character", min.len = 1)
  assertIntegerish(p, lower = 1, any.missing = FALSE, len = 1)
  assertNumeric(penalty, lower = 0, any.missing = FALSE, len = 1)
  
  # not all feature sets can be empty
  features_chosen <- unlist(lapply(features, length))
  if(!(sum(features_chosen == 0) < length(features))) return(NA) ###
  
  # number of feature sets
  n <- length(features)
  
  # all selected features
  all_features <- unlist(features)
  
  # frequencies of features
  freq <- table(all_features)
  
  # median number of selected features
  median_f <- median(unlist(lapply(features, length)))
  
  # calculate index
  part1 <- sum(freq)/(n * length(unique(all_features)))
  part2 <- penalty * median_f/p
  davis <- max(c(0, part1-part2))
  
  #if(is.na(davis)) davis <- 0
  return(davis)
}

dice <- function(features, ...)
{
  # make sure that input has desired form
  assertList(x = features, types = "character", min.len = 1)
  
  # not more than one feature set can be empty
  features_chosen <- unlist(lapply(features, length))
  if(!(sum(features_chosen == 0) < 2)) return(NA) ###
  
  # number of feature sets
  n <- length(features)
  
  # sizes of feature sets
  feat_sizes <- unlist(lapply(features, length))
  
  # calculate index
  index_values <- unlist(sapply(1:(n-1), function(i)
  {
    sapply((i+1):n, function(j)
    {
      intersection <- length(intersect(features[[i]], features[[j]]))
      return(2*intersection / (feat_sizes[i] + feat_sizes[j]))
    })
  }))
  
  dice <- sum(unlist(index_values)) * 2/(n*(n-1))
  #if(is.na(dice)) dice <- 0
  return(dice)
}

jaccard <- function(features, ...)
{
  # make sure that input has desired form
  assertList(x = features, types = "character", min.len = 1)
  
  # not more than one feature set can be empty
  features_chosen <- unlist(lapply(features, length))
  if(!(sum(features_chosen == 0) < 2)) return(NA) ###
  
  # number of feature sets
  n <- length(features)
  
  # calculate sizes of intersections
  intersections <- unlist(sapply(1:(n-1), function(i)
  {
    sapply((i+1):n, function(j) length(intersect(features[[i]], features[[j]])))
  }))
  
  # calculate sizes of unions
  unions <- unlist(sapply(1:(n-1), function(i)
  {
    sapply((i+1):n, function(j) length(union(features[[i]], features[[j]])))
  }))
  
  # calculate Jaccard index
  jaccard <- 2/(n*(n-1)) * sum(intersections/unions)
  
  #if(is.na(jaccard)) jaccard <- 0
  return(jaccard)
}

lustgarten <- function(features, p, ...)
{
  # make sure that input has desired form
  assertList(x = features, types = "character", min.len = 1)
  assertIntegerish(p, lower = 1, any.missing = FALSE, len = 1)
  
  # no feature set can be empty and not more than one can contain all p features
  features_chosen <- unlist(lapply(features, length))
  if(!(sum(features_chosen == 0) == 0)) return(NA) ###
  if(!(sum(features_chosen == p) < 2)) return(NA) ###
  
  # number of feature sets
  n <- length(features)
  
  # sizes of feature sets
  feat_sizes <- unlist(lapply(features, length))
  
  # calculate index
  index_values <- unlist(sapply(1:(n-1), function(i)
  {
    sapply((i+1):n, function(j)
    {
      intersection <- length(intersect(features[[i]], features[[j]]))
      part1 <- intersection - feat_sizes[i]*feat_sizes[j]/p
      part2 <- min(c(feat_sizes[i], feat_sizes[j])) - max(c(0, feat_sizes[i] + feat_sizes[j] - p))
      return(part1/part2)
    })
  }))
  
  lustgarten <- sum(unlist(index_values)) * 2/(n*(n-1))
  
  #if(is.na(lustgarten)) lustgarten <- 0
  return(lustgarten)
}

novovicova <- function(features, ...)
{
  # make sure that input has desired form
  assertList(x = features, types = "character", min.len = 1)
  
  # not all feature set can be empty
  features_chosen <- unlist(lapply(features, length))
  if(!(sum(features_chosen == 0) < length(features))) return(NA) ###
  
  # number of feature sets
  n <- length(features)
  
  # all selected features
  all_features <- unlist(features)
  
  # frequencies of features
  freq <- table(all_features)
  q <- sum(freq)
  
  # calculate index
  novovicova <- sum(freq * log2(freq)) / (q*log2(n))
  
  #if(is.na(novovicova)) novovicova <- 0
  return(novovicova)
}

ochiai <- function(features, ...)
{
  # make sure that input has desired form
  assertList(x = features, types = "character", min.len = 1)
  
  # no feature set can be empty
  features_chosen <- unlist(lapply(features, length))
  if(!(sum(features_chosen == 0) < 1)) return(NA) ###
  
  # number of feature sets
  n <- length(features)
  
  # sizes of feature sets
  feat_sizes <- unlist(lapply(features, length))
  
  # calculate index
  index_values <- unlist(sapply(1:(n-1), function(i)
  {
    sapply((i+1):n, function(j)
    {
      intersection <- length(intersect(features[[i]], features[[j]]))
      return(intersection / sqrt(feat_sizes[i] * feat_sizes[j]))
    })
  }))
  
  ochiai <- sum(unlist(index_values)) * 2/(n*(n-1))
  
  #if(is.na(ochiai)) ochiai <- 0
  return(ochiai)
}

somol <- function(features, p, ...)
{
  # make sure that input has desired form
  assertList(x = features, types = "character", min.len = 1)
  assertIntegerish(p, lower = 1, any.missing = FALSE, len = 1)
  
  # not all feature set can be empty or contain all p features
  features_chosen <- unlist(lapply(features, length))
  if(!(sum(features_chosen == 0) < length(features))) return(NA) ###
  if(!(sum(features_chosen == p) < length(features))) return(NA) ###
  
  # number of feature sets
  n <- length(features)
  
  # all selected features
  all_features <- unlist(features)
  
  # frequencies of features
  freq <- table(all_features)
  q <- sum(freq)
  
  # calculate index
  c_min <- (q^2 - p*(q - q%%p) - (q%%p)^2) / (p*q*(n-1))
  c_max <- ((q%%n)^2 + q*(n-1)- (q%%n)*n) / (q*(n-1))
  somol <- (sum(freq * (freq-1)) / (q*(n-1)) - c_min) / (c_max - c_min)
  
  #if(is.na(somol)) somol <- 0
  return(somol)
}

zucknick <- function(features, correlations, ...)
{
  # make sure that input has desired form
  assertList(x = features, types = "character", min.len = 1)
  assertMatrix(x = correlations, mode = "numeric", any.missing = FALSE, min.rows = 1, min.cols = 1)
  assert(nrow(correlations) == ncol(correlations))
  assertNumeric(x = correlations, lower = -1, upper = 1)
  
  # not more than one feature set can be empty
  features_chosen <- unlist(lapply(features, length))
  if(!(sum(features_chosen == 0) < 2)) return(NA) ###
  
  # number of feature sets
  n <- length(features)
  
  # sizes of feature sets
  feat_sizes <- unlist(lapply(features, length))
  
  # calculate sizes of intersections
  intersections <- unlist(sapply(1:(n-1), function(i)
  {
    sapply((i+1):n, function(j) length(intersect(features[[i]], features[[j]])))
  }))
  
  # calculate sizes of unions
  unions <- unlist(sapply(1:(n-1), function(i)
  {
    sapply((i+1):n, function(j) length(union(features[[i]], features[[j]])))
  }))
  
  # calculate added correlations
  features_indices <- lapply(features, function(x)
    {
    sapply(x, function(f) which(colnames(correlations) == f))
  })
  
  cs <- unlist(sapply(1:(n-1), function(i)
  {
    sapply((i+1):n, function(j)
    {
      j_without_i <- setdiff(features_indices[[j]], features_indices[[i]])
      i_without_j <- setdiff(features_indices[[i]], features_indices[[j]])
      
      if(length(j_without_i) == 0 && length(i_without_j) == 0)
      {
        c <- 0
      }
      else if(length(j_without_i) == 0)
      {
        grid_j <- as.matrix(expand.grid(features_indices[[j]], i_without_j))
        cors_j <- drop(correlations[grid_j])
        r <- median(cors_j)
        cors_j <- ifelse(cors_j > r, cors_j, 0) / feat_sizes[i]
        c <- sum(cors_j)
      }
      else if(length(i_without_j) == 0)
      {
        grid_i <- as.matrix(expand.grid(features_indices[[i]], j_without_i))
        cors_i <- drop(correlations[grid_i])
        r <- median(cors_i)
        cors_i <- ifelse(cors_i > r, cors_i, 0) / feat_sizes[j]
        c <- sum(cors_i)
      }
      else
      {
        grid_i <- as.matrix(expand.grid(features_indices[[i]], j_without_i))
        grid_j <- as.matrix(expand.grid(features_indices[[j]], i_without_j))
        cors_i <- drop(correlations[grid_i])
        cors_j <- drop(correlations[grid_j])
        r <- median(c(cors_i, cors_j))
        cors_i <- ifelse(cors_i > r, cors_i, 0) / feat_sizes[j]
        cors_j <- ifelse(cors_j > r, cors_j, 0) / feat_sizes[i]
        c <- sum(cors_i) + sum(cors_j)
      }
      return(c)
    })
  }))
  
  # calculate index
  zucknick <- 2/(n*(n-1)) * sum((intersections + cs)/unions)
  
  #if(is.na(zucknick)) zucknick <- 0
  return(zucknick)
}

cor_pearson <- function(features, p, ...)
{
  # make sure that input has desired form
  assertList(x = features, types = "character", min.len = 1)
  assertIntegerish(p, lower = 1, any.missing = FALSE, len = 1)
  
  # number of feature sets
  n <- length(features)
  
  # write feature sets in matrix coded by 0 and 1
  all_selected_features <- unique(unlist(features))
  feat_mat <- matrix(0, ncol = n, nrow = p)
  for(i in 1:n)
  {
    feature_indices <- which(all_selected_features %in% features[[i]])
    feat_mat[feature_indices, i] <- 1
  }
  
  # calculate correlations
  cors <- unlist(sapply(1:(n-1), function(i)
  {
    sapply((i+1):n, function(j) cor(feat_mat[, i], feat_mat[, j]))
  }))
  
  return(mean(cors))
}

all_stability_measures <- function(features, p, correlations)
{
  stability <- numeric(12)
  stability[1] <- cor_pearson(features, p)
  stability[2] <- davis(features, p, 0)
  stability[3] <- davis(features, p, 1)
  stability[4] <- davis(features, p, 2)
  stability[5] <- davis(features, p, 10)
  stability[6] <- dice(features)
  stability[7] <- jaccard(features)
  stability[8] <- lustgarten(features, p)
  stability[9] <- novovicova(features)
  stability[10] <- ochiai(features)
  stability[11] <- somol(features, p)
  stability[12] <- zucknick(features, correlations)
  names(stability) <- c("cor_pearson", "davis_0", "davis_1", "davis_2", "davis_10", "dice",
                        "jaccard", "lustgarten", "novovicova", "ochiai", "somol", "zucknick")
  return(stability)
}