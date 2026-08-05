#* INFO Show use of Train-Test-Split for the same brain dataset as brain-cross-predict-ex.r
# install.packages("remotes")
# remotes::install_github("swanberry/seqclass")
library(seqclass)
library(data.table)
library(dplyr)

# process is very similar to brain-cross-predict-ex.r
if (!file.exists("GSE205450_counts.table.txt.gz")) {
  print("Data missing in directory-- downloading")
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205450/suppl/GSE205450_counts.table.txt.gz",
    destfile = "GSE205450_counts.table.txt.gz",
    mode = "wb"
  )
  print("File downloaded")
}

data <- as.data.frame(fread("GSE205450_counts.table.txt.gz"))[-1]

vec <- grepl(".*PD.*", colnames(data))
class <- c("F", "T")[vec + 1]

# highlight the use of tts and show the confusion matrix output
ls <- generate_model(data = data, class = class, train_test_split = TRUE)
ls$confusion_matrix
