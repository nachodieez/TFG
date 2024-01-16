library(limma)
targets <- readTargets("targets.txt")
RG <- read.maimages(targets, source = "agilent")
RG <- backgroundCorrect(RG, method = "normexp", offset=16)
MA <- normalizeWithinArrays(RG, method = "loess")
options("max.print"=1E9)
options("width"=10000)

# Normalized matrix
sink("a_GSE58791_JQ1_Expression_MA.txt")
print(MA$A)
sink()

# Genes list
sink("b_GSE58791_JQ1_Genes_MA.txt")
print(MA$genes)
sink()

## A:numeric matrix containing the A-values (average log-2 expression values)

