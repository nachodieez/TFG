# Read the Target file (experimental description) 
library(limma)
targets <- readTargets()

# Reading the Images
RG  <- read.maimages(targets$FileName, source = "agilent", green.only = TRUE)
RGb <- backgroundCorrect(RG, method = "normexp")
## Within for two colours
# MA <- normalizeWithinArrays(RGb, method="loess")
##
Qnorm <- normalizeBetweenArrays(RGb, method = "quantile")


