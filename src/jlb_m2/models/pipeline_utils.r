# Making task
making_task <- function(data, dataname, target, positive) {
    data <- as.data.frame(data)
    data[, target] <- as.factor(data[, target])
    # Make task
    task <- TaskClassif$new(id = dataname,
                        backend = data,
                        target = target,
                        positive = positive)
    return(task)
}

# Preprocess
preprocess <- function(task,
                       removeConstant,
                       normalize,
                       filterFeatures) {

  # Remove constant features
  if (removeConstant == TRUE) {
    rcf <- po("removeconstants")
    task <- rcf$train(list(task = task))$output
    print("Constant features have been removed!")
  }
  # Normalizing features
  if (normalize == TRUE) {
    nf <- po("scale")
    task <- nf$train(input = list(task))$output
    print("Features have been normalized!")
  }
  # Filter features
  if (filterFeatures == TRUE) {
    filter <- po("filter",
                 filter = flt("kruskal_test"),
                 filter.frac = 0.01)
    task <- filter$train(list(task = task))$output
    print("Features have been filtered!")
  }

  return(task)
}