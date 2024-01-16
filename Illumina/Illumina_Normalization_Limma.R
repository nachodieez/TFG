library(limma)
maqc <- read.ilmn(files="GSE119842_01_gex_matrix_non-normalized.txt", probeid="ID_REF", other.columns=c("Detection"))
maqc <- read.ilmn(files="GSE120991_non-normalized_data.txt", probeid = "ID_REF", expr = "IB")
# neqc functions performs normexp background correction using negative controls, then quantile normalizes and finally log2 transforms
# The data file has been edited and no longer has standard Illumina column headings, so you have to tell read.ilmn how to recognize the probeid column and the expression columns. 
maqc.norm <- neqc(maqc)
options("max.print" = 1E9)
options("width" = 10000)
sink("c.txt")
print(as.matrix(maqc.norm))
sink()

#maqc$other$Detection de texto a numeritos
colnames <- maqc$other$Detection %>% colnames()
ncol <- maqc$other$Detection %>% ncol
nrow <- maqc$other$Detection %>% nrow
maqc$other$Detection <- mapply(maqc$other$Detection, FUN = as.numeric)
maqc$other$Detection <- matrix(data = maqc$other$Detection, ncol = ncol, nrow = nrow)
colnames(maqc$other$Detection) <- colnames
maqc$other$Detection <- ifelse(is.na(maqc$other$Detection), 0, maqc$other$Detection)
