library(readxl)

setwd("C:/Users/nacho/OneDrive/Universidad/Cuarto/TFG/Meta-analisis")
significativos <- read_excel("significativos.xlsx", 
                             skip = 1)

experimentos <- as.factor(c(1,1,2,3,3,4,5,6,6,7,7,8,9,9,10,10,11,12,12,12,13,13,13,13,14,14,14,15,
                          16,17,17,18,19,
                          20,21,22,23,23,23,23,23,23,24,25,26,26,27,27,
                          28,28,29,30,30,31,32,33,33))

individuos <- data.frame(experimento = colnames(significativos), grupo = experimentos)
individuos

genes_sig <- list()
for(i in levels(experimentos)) {
  a <- significativos[,experimentos == i]
  #a <- na.omit(a)
  b <- c()
  for(j in 1:length(a)){
    b <- c(b, a[[j]])
  }
  b <- na.omit(b)
  b <- unique(b)
  genes_sig[as.numeric(i)] <- list(b)
}

total_genes <- c()
for(i in 1:length(genes_sig)){
  total_genes <- c(total_genes, genes_sig[[i]])
}
total_genes <- unique(total_genes)

repeticiones <- data.frame(n = rep(0,length(total_genes)))
rownames(repeticiones) <- total_genes


for(i in 1:length(genes_sig)){
  for(j in 1:length(genes_sig[[i]])){
    repeticiones[genes_sig[[i]][j],] <- repeticiones[genes_sig[[i]][j],] +1
  }
}

sum(repeticiones>1)

ngenes_repetidos <- c()
for (i in 1:length(levels(experimentos))) {
  ngenes_repetidos <- c(ngenes_repetidos, sum(repeticiones>=i))
}

rownames(repeticiones)[repeticiones>19]
rownames(repeticiones)[repeticiones>27]

data.frame(repeticiones = ngenes_repetidos)
write.csv2(data.frame(gen = rownames(repeticiones)[repeticiones>1]), "genecillos.csv")


nombres <- c()
for(i in 1:length(genes_sig)){
  nombres <- c(nombres, unlist(significativos[experimentos==i][1,1]))
}


x <- c()
for(i in 1:length(genes_sig)){
  x <- c(x, length(genes_sig[[i]]))
}
sig <- data.frame(nombres = nombres,genes_sig = x)





