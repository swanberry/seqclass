#* Use of explicit features for brian-cross-predict-ex.r
# install.packages("remotes")
# remotes::install_github("swanberry/seqclass")
library(seqclass)
library(data.table)
library(dplyr)
library(caret)

if (!file.exists("GSE205450_counts.table.txt.gz")) {
  print("Data missing in directory-- downloading")
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205450/suppl/GSE205450_counts.table.txt.gz",
    destfile = "GSE205450_counts.table.txt.gz",
    mode = "wb"
  )
  print("File downloaded")
}

#! WARN The use of rownames is mandatory-- see below for how to deal with this
df <- as.data.frame(fread("GSE205450_counts.table.txt.gz")) |>
  filter(!is.na(Gene_symbol))
genes <- df[[1]]
df <- df[-1]
rownames(df) <- genes
df_cau <- df |>
  select(matches(".*CAU.*"))
df_put <- df |>
  select(matches(".*PUT.*"))

vec <- grepl(".*PD.*", colnames(df_cau))
class <- c("F", "T")[vec + 1]

# see how to use features here
model <- generate_model(data = df_cau, class = class, features = 100)
predictions <- use_model(
  data = df_put,
  model = model$model,
  features = model$features
)

# generate a confusionMatrix-- note the use of features here
ans <- grepl(".*PD.*", colnames(df_put[model$features, ]))
ans <- c("F", "T")[ans + 1]
ans <- factor(ans, levels = c("F", "T"))
confusionMatrix(data = predictions, reference = ans, positive = "T")
