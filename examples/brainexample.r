#* train a Parkinson's classifier on CAU samples and test it on PUT samples
#! WARNING the model may take ~2-5mins to generate

# install.packages("remotes")
# remotes::install_github("qryce-01/seqclass")
library(seqclass)
library(data.table)
library(dplyr)
library(caret)

# download the required dataset
if (!file.exists("GSE205450_counts.table.txt.gz")) {
  print("Data missing in directory-- downloading")
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205450/suppl/GSE205450_counts.table.txt.gz",
    destfile = "GSE205450_counts.table.txt.gz",
    mode = "wb"
  )
  print("File downloaded")
}

# read the dataset and select caudate and putamen separately
df <- as.data.frame(fread("GSE205450_counts.table.txt.gz"))[-1]
df_cau <- df |>
  select(matches(".*CAU.*"))
df_put <- df |>
  select(matches(".*PUT.*"))

# true or false if a column has parkinson's
vec <- grepl(".*PD.*", colnames(df_cau))
# convert to the "T"/"F" default using some clever indexing
class <- c("F", "T")[vec + 1]

model <- generate_model(data = df_cau, class = class)
predictions <- use_model(data = df_put, model = model)

# generate a confusionMatrix 
ans <- grepl(".*PD.*", colnames(df_put))
ans <- c("F", "T")[ans + 1]
ans <- factor(ans, levels = c("F", "T"))
confusionMatrix(data = predictions, reference = ans, positive = "T")
