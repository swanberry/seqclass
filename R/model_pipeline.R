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
#' @param positive The string meaning positive for the class. Note this cannot be "TRUE".
#' @param use_svmRadial Determines whether or not an svmRadial model is used instead of voomNSC
#' @param voom_ctrl A control function to modify the model for voomNSC models
#' @param svm_ctrl A control function to modify the model for SVM models (use only if use_svmRadial is TRUE)
#' @param train_test_split Determine whether or not to split (70:30) into train and test data
#' @param normalize For voom-based clasifiers. Character string indicating type of normalization-- 'deseq', 'tmm', or 'none'
#' @param negative Required if train_test_split is true-- the string meaning negative for the class; default "F"
#' @param preProcessing For caret-based classifiers. Name of preprocessing method-- eg "deseq-vst", "deseq-rlog", "logcpm"
#' @return An MLSeq model that can be passed into use_model. If train-test-split is true, then this outputs a named list with components:
#' @return \item{model}{The trained model}
#' @return \item{confusion_matrix}{The confusion matrix from the tts}
#' @export
generate_model <- function(
  data,
  class,
  positive = "T",
  use_svmRadial = FALSE,
  voom_ctrl = MLSeq::voomControl(
    method = "repeatedcv",
    number = 5,
    repeats = 2,
    tuneLength = 7
  ),
  svm_ctrl = caret::trainControl(
    method = "repeatedcv",
    number = 2,
    repeats = 2,
    classProbs = TRUE
  ),
  train_test_split = FALSE,
  negative = "F",
  normalize = "deseq",
  preProcessing = "deseq-vst"
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
      use_svmRadial = use_svmRadial,
      voom_ctrl = voom_ctrl,
      svm_ctrl = svm_ctrl,
      train_test_split = FALSE,
      normalize = normalize,
      preProcessing = preProcessing
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

  if (use_svmRadial) {
    model <- MLSeq::classify(
      data = data.S4,
      method = "svmRadial",
      ref = positive,
      control = svm_ctrl,
      preProcessing = preProcessing
    )
  } else {
    model <- MLSeq::classify(
      data = data.S4,
      method = "voomNSC",
      ref = positive,
      control = voom_ctrl,
      normalize = normalize
    )
  }

  return(model)
}
