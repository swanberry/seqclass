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
#' @param positive The string meaning positive for the class; default "T"
#' @param ctrl A control function to modify the model; by default repeatedcv, 2 repeats x 5, tune Length 7
#' @return An MLSeq model that can be passed into use_model
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
  )
) {
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
