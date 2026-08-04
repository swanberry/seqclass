#' Use an existing model
#'
#' @param data A data frame or count matrix (Features as rows, samples as columns)
#' @param model A model generated from generate_model
#' @return A factor of predictions
#' @export
use_model <- function(data, model) {
  matrix <- as.matrix(data)
  rownames(matrix) <- rownames(data)
  data <- matrix
  dummy <- S4Vectors::DataFrame(colnames(data))
  data.S4 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data,
    colData = dummy,
    design = ~1
  )
  return(MLSeq::predict(model, data.S4))
}

#' Generate a voomNSC model using MLSeq
#'
#' @param data A data frame or count matrix (Features as rows, samples as columns)
#' @param class A binary vector of strings containing the condition factor
#' @param positive The string meaning positive for the class; default "T". Note this cannot be "TRUE".
#' @param ctrl A control function to modify the model; by default repeatedcv, 2 repeats x 5, tune Length 7
#' @param train_test_split Determine whether or not to split (70:30) into train and test data
#' @param negative Required if train_test_split is true-- the string meaning negative for the class; default "F"
#' @return An MLSeq model that can be passed into use_model. If train-test-split is true, then this outputs a named list with components:
#' @return \item{model}{The trained model}
#' @return \item{confusion_matrix}{The confusion matrix from the tts}
#' @export
generate_model <- function(
  data,
  class,
  positive = "T",
  ctrl = MLSeq::voomControl(
    method = "repeatedcv",
    number = 5,
    repeats = 2,
    tuneLength = 7
  ),
  train_test_split = FALSE,
  negative = "F"
) {
  # use the abstraction provided to cleanly implement tts
  if (train_test_split) {
    nTest <- floor(0.3 * ncol(data))
    ind <- sample(ncol(data), nTest)
    data_tr <- data[, ind]
    data_ts <- data[, -ind]
    class_tr <- class[ind]
    class_ts <- class[-ind]
    model <- generate_model(
      data = data_tr,
      class = class_tr,
      positive = positive,
      ctrl = ctrl
    )
    predictions <- use_model(data = data_ts, model = model)
    class_ts <- factor(class_ts, levels = c(negative, positive))
    confusion_matrix <- caret::confusionMatrix(
      data = predictions,
      reference = class_ts,
      positive = positive
    )
    return(list("model" = model, "confusion_matrix" = confusion_matrix))
  }

  class <- S4Vectors::DataFrame(
    condition = factor(class)
  )

  data <- as.matrix(data + 1)

  # convert to DESeq2 datasets
  data.S4 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data,
    colData = class,
    design = ~condition
  )

  model <- MLSeq::classify(
    data = data.S4,
    method = "voomNSC",
    ref = positive,
    control = ctrl,
    normalize = "deseq"
  )

  return(model)
}
