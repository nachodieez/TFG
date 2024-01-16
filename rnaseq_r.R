library(edgeR)
library(tidyverse)

options("max.print"=1E9)

counts <- read.table("C:/Users/nacho/OneDrive/Universidad/Cuarto/TFG/Meta-analisis/RNA-seq/GSE104714/counts.txt", quote="\"", comment.char="")
#counts <- counts[,c()]

groups <- as.factor(rep(c(), each=2))

people <- data.frame(individuo = colnames(counts), grupo = groups); people

d0 <- DGEList(counts, group = groups)
d0 <- calcNormFactors(d0)
d0

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

plotMDS(d0, col = as.numeric(groups))

mm <- model.matrix(~0 + groups)
y <- voom(d0, mm, plot = T)

fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupsControl - groupsDex, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

length(which(top.table$adj.P.Val < 0.05))

nombres <- rownames(top.table[top.table$adj.P.Val<0.05 ,])
length(nombres)

index_1 <- 1 #index del primer grupo del contraste el levels(groups)
index_2 <- 2 #index del segundo grupo del contraste
if(length(nombres)!=0){
  genes <- 2^y$E[nombres,]
  medias_A <- c()
  medias_B <- c()
  desviaciones_A <- c()
  desviaciones_B <- c()
  for(gen in nombres){
    medias_A <- c(medias_A,
                  genes[gen,groups==levels(groups)[index_1]] %>% as.numeric %>% mean)
    medias_B <- c(medias_B,
                  genes[gen,groups==levels(groups)[index_2]] %>% as.numeric %>% mean)
    desviaciones_A <- c(desviaciones_A,
                        genes[gen,groups == levels(groups)[index_1]] %>% as.numeric %>% sd)
    desviaciones_B <- c(desviaciones_B,
                        genes[gen,groups == levels(groups)[index_2]] %>% as.numeric %>% sd)
  }
  
  df <- data.frame(A = medias_A,
                   B = desviaciones_A,
                   C = medias_B,
                   D = desviaciones_B,
                   E = rep(plyr::count(groups)[index_1,2],length(nombres)),
                   F = rep(plyr::count(groups)[index_2,2],length(nombres)))
  
  rownames(df) <- nombres
  
  colnames(df) <- c(       paste0("Media ", levels(groups))[index_1],
                           paste0("SD ", levels(groups)[index_1]),
                           paste0("Media ", levels(groups))[index_2],
                           paste0("SD ", levels(groups)[index_2]),
                           paste0("n ", levels(groups)[index_1]),
                           paste0("n ", levels(groups)[index_2])
  )
  
  write.csv2(df, file ="significativos.csv")
  
}
if(T){
  
  nombres <- rownames(y$E)
  
  genes <- 2^y$E[nombres,]
  medias_A <- c()
  medias_B <- c()
  desviaciones_A <- c()
  desviaciones_B <- c()
  for(gen in nombres){
    medias_A <- c(medias_A,
                  genes[gen,groups==levels(groups)[index_1]] %>% as.numeric %>% mean)
    medias_B <- c(medias_B,
                  genes[gen,groups==levels(groups)[index_2]] %>% as.numeric %>% mean)
    desviaciones_A <- c(desviaciones_A,
                        genes[gen,groups == levels(groups)[index_1]] %>% as.numeric %>% sd)
    desviaciones_B <- c(desviaciones_B,
                        genes[gen,groups == levels(groups)[index_2]] %>% as.numeric %>% sd)
  }
  
  df <- data.frame(A = medias_A,
                   B = desviaciones_A,
                   C = medias_B,
                   D = desviaciones_B,
                   E = rep(plyr::count(groups)[index_1,2],length(nombres)),
                   F = rep(plyr::count(groups)[index_2,2],length(nombres)))
  
  rownames(df) <- nombres
  
  colnames(df) <- c(       paste0("Media ", levels(groups))[index_1],
                           paste0("SD ", levels(groups)[index_1]),
                           paste0("Media ", levels(groups))[index_2],
                           paste0("SD ", levels(groups)[index_2]),
                           paste0("n ", levels(groups)[index_1]),
                           paste0("n ", levels(groups)[index_2])
  )
  
  write.csv2(df, file ="genes.csv")
}