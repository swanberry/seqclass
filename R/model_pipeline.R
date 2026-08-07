#' Use an existing model
#'
#' @param data A data frame or count matrix (Features as rows, samples as columns)
#' @param model A model generated from generate_model
#' @param features The features output from generate_model (ls$features)
#' @param yield_probs (Only for svmModels) show probabilities of "TRUE" instead of binary classifications
#' @param add_one Boolean: Adds one to the data to avoid too many zeroes error; helpful if low feature count
#' @return A factor of predictions
#' @export
use_model <- function(
  data,
  model,
  features = "all",
  add_one = FALSE,
  yield_probs = FALSE
) {
  if (!identical(features, "all")) {
    data <- data[features, ]
  }
  if (add_one) {
    matrix <- as.matrix(data + 1)
  } else {
    matrix <- as.matrix(data)
  }

  rownames(matrix) <- rownames(data)
  data <- matrix
  dummy <- S4Vectors::DataFrame(colnames(data))
  data.S4 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data,
    colData = dummy,
    design = ~1
  )

  if (yield_probs) {
    return(MLSeq::predict(model, data.S4, type = "prob"))
  }
  return(MLSeq::predict(model, data.S4))
}

#' Generate a model using MLSeq
#'
#' @param data A data frame or count matrix (Features as rows, samples as columns)
#' @param class A binary vector of strings containing the condition factor
#' @param features Integer indicating how many of the best features to use-- also accepts arguments "all" and "tenthOfN"
#' @param add_one Boolean: Adds one to the data to avoid the too many zeroes error; helpful if using not many features.
#' @param positive The string meaning positive for the class. Note this cannot be "TRUE".
#' @param use_svmRadial Determines whether or not an svmRadial model is used instead of voomNSC (needed for roc dataframe)
#' @param voom_ctrl A control function to modify the model for voomNSC models
#' @param svm_ctrl A control function to modify the model for SVM models (use only if use_svmRadial is TRUE)
#' @param train_test_split Determine whether or not to split (70:30) into train and test data (tts)
#' @param roc_curve_samples For tts-- how many thresholds to test (default 1000) for the ROC df
#' @param normalize For voom-based clasifiers. Character string indicating type of normalization-- 'deseq', 'tmm', or 'none'
#' @param negative Required if train_test_split is true-- the string meaning negative for the class; default "F"
#' @param preProcessing For caret-based classifiers. Name of preprocessing method-- eg "deseq-vst", "deseq-rlog", "logcpm"
#' @return An MLSeq model that can be passed into use_model. Returns a named list with components:
#' @return \item{model}{The trained model}
#' @return \item{features}{The features selected-- make sure to pass this into use_model if you selected features}
#' @return \item{confusion_matrix}{The confusion matrix from the tts-- only outputted if tts is TRUE}
#' @return \item{roc_df}{A dataframe from which an roc curve can be easily plotted-- only outputted if using an SVM model}
#' @export
generate_model <- function(
  data,
  class,
  features = "all",
  positive = "T",
  add_one = FALSE,
  use_svmRadial = FALSE,
  voom_ctrl = MLSeq::voomControl(
    method = "repeatedcv",
    number = 3,
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
  roc_curve_samples = 1000,
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
      features = features,
      positive = positive,
      add_one = add_one,
      use_svmRadial = use_svmRadial,
      voom_ctrl = voom_ctrl,
      svm_ctrl = svm_ctrl,
      train_test_split = FALSE,
      normalize = normalize,
      preProcessing = preProcessing
    )

    # separate probs and predictions logic
    predictions <- probs <- use_model(
      data = data_ts,
      model = model$model,
      features = model$features,
      add_one = add_one,
      yield_probs = FALSE
    )
    class_ts <- factor(class_ts, levels = c(negative, positive))
    confusion_matrix <- caret::confusionMatrix(
      data = predictions,
      reference = class_ts,
      positive = positive
    )

    # probs only makes sense if we use svmRadial
    if (use_svmRadial) {
      df <- data.frame(
        true_pos_rate = numeric(),
        false_pos_rate = numeric(),
        threshold = numeric()
      )
      class <- as.character(class_ts)
      positives <- length(class[class == positive])
      negatives <- length(class[class == negative])
      probs <- use_model(
        data = data_ts,
        model = model$model,
        features = model$features,
        add_one = add_one,
        yield_probs = TRUE
      )[, positive] # so is just a list of probabilities

      for (i in 1:roc_curve_samples) {
        tpc <- 0
        fpc <- 0

        threshold <- i / roc_curve_samples
        boolvec <- ifelse(probs >= threshold, positive, negative)
        for (j in seq_along(boolvec)) {
          prediction <- boolvec[j]
          actual <- class[j]
          if (prediction == positive && actual == positive) {
            # if prediction is true and actual is true (true positive)
            tpc <- tpc + 1
          }
          if (prediction == positive && actual == negative) {
            # if prediction is true and actual is false (false positive)
            fpc <- fpc + 1
          }
        }
        tpr <- tpc / positives
        fpr <- fpc / negatives
        new_row <- data.frame(
          true_pos_rate = tpr,
          false_pos_rate = fpr,
          threshold = threshold
        )
        df <- rbind(df, new_row)
      }

      return(list(
        "model" = model,
        "features" = features,
        "confusion_matrix" = confusion_matrix,
        "roc_df" = df
      ))
    }
    return(list(
      "model" = model,
      "features" = features,
      "confusion_matrix" = confusion_matrix
    ))
  }

  class <- S4Vectors::DataFrame(
    condition = factor(class)
  )

  if (add_one) {
    data <- as.matrix(data + 1)
  } else {
    data <- as.matrix(data)
  }

  print("Converting to DESeq2 DS-- this may take a while")
  print("If this takes a very long time, try passing an argument into features")
  # convert to DESeq2 datasets
  data.S4 <- DESeq2::DESeqDataSetFromMatrix(
    countData = data,
    colData = class,
    design = ~condition
  )

  # handle feature use
  if (features != "all") {
    if (features == "tenthOfN") {
      features <- ceiling(ncol(data) * 0.1)
    }
    dds <- DESeq2::DESeq(data.S4)
    res <- DESeq2::results(dds)
    sres <- res[order(res$padj), ]
    features <- rownames(sres[1:features, ])
    data.S4 <- data.S4[features]
  }

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
  return(list("model" = model, "features" = features))
}
