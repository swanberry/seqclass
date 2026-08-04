# seqclass

`seqclass` provides lightweight, streamlined wrapper functions for training and predicting with RNA-seq classification models using `MLSeq` and `DESeq2`. It removes the tedious boilerplate of coercing data frames, formatting S4 datasets, and aligning factor levels.

## installation

you can install the development version of `seqclass` from GitHub:

```r
# install.packages("remotes")
remotes::install_github("qryce-01/seqclass")
```

## usage

```r
library(seqclass)

# train a voomNSC model using your count matrix and class labels
model <- generate_model(
  data = my_counts, 
  class = my_classes
)

# predict on new, similar data
predictions <- use_model(my_new_counts, model)
print(predictions)
```
