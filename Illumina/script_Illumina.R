library(GEOquery)
library(siggenes)
library(car)
library(ggpubr)
library(ggfortify)
library(limma)
library(Glimma)
library(edgeR)
library(gplots)
library(dplyr)

codigo <- ""
datos0 <- read.csv("c2.txt", sep="")
head(datos0)

grupos0 <- as.factor()

individuos <- data.frame(individuo = colnames(datos0), grupo = grupos0)
individuos

design <- model.matrix(~0+grupos)
colnames(design) <- levels(grupos)

contr.matrix <- makeContrasts(
  AvsB = A - B,
  levels = grupos
)

fit <- lmFit(datos, design)
fit <- contrasts.fit(fit, contrasts = contr.matrix)
summary(decideTests(fit))
AvsB <- topTreat(eBayes(fit), n=Inf)
genes <- datos0$gen
nombs <- data.frame(sonda = rownames(datos), gen = genes)
genes <- c()
for (sonda in rownames(AvsB)) {
  genes <- c(genes, nombs$gen[nombs$sonda==sonda])
}
AvsB$gen <- genes

tabla <- c()
for (nombre in base::unique(AvsB$gen)){
  g <- AvsB[AvsB$gen==nombre,]
  AvsB <- AvsB[!AvsB$gen==nombre,]
  g <- g[which.max(g$adj.P.Val),]
  tabla <- rbind(tabla,g)
}

AvsB <- tabla
datos <- datos[rownames(AvsB),]
rownames(AvsB) <- AvsB$gen
rownames(datos) <- AvsB$gen
AvsB <- subset(AvsB,select = -gen)
AvsB <- AvsB[order(tabla$adj.P.Val),]
head(AvsB)
contraste <- deparse(substitute(AvsB))


head(AvsB)

nombres <- rownames(AvsB[AvsB$adj.P.Val<0.05 ,])
length(nombres)

#write.table(file = paste0(contraste,"_significativos.txt"), x = nombres, sep = "\n", row.names = F, col.names = F)

index_1 <- 1 #index del primer grupo del contraste el levels(groups)
index_2 <- 2 #index del segundo grupo del contraste
if(length(nombres)!=0){
  genes <- 2^datos[nombres,]
  medias_A <- c()
  medias_B <- c()
  desviaciones_A <- c()
  desviaciones_B <- c()
  for(gen in nombres){
    medias_A <- c(medias_A,
                  genes[gen,grupos==levels(grupos)[index_1]] %>% as.numeric %>% mean)
    medias_B <- c(medias_B,
                  genes[gen,grupos==levels(grupos)[index_2]] %>% as.numeric %>% mean)
    desviaciones_A <- c(desviaciones_A,
                        genes[gen,grupos == levels(grupos)[index_1]] %>% as.numeric %>% sd)
    desviaciones_B <- c(desviaciones_B,
                        genes[gen,grupos == levels(grupos)[index_2]] %>% as.numeric %>% sd)
  }
  
  df <- data.frame(A = medias_A,
                   B = desviaciones_A,
                   C = medias_B,
                   D = desviaciones_B,
                   E = rep(plyr::count(grupos)[index_1,2],length(nombres)),
                   F = rep(plyr::count(grupos)[index_2,2],length(nombres)))
  
  rownames(df) <- nombres
  
  colnames(df) <- c(       paste0("Media ", levels(grupos))[index_1],
                           paste0("SD ", levels(grupos)[index_1]),
                           paste0("Media ", levels(grupos))[index_2],
                           paste0("SD ", levels(grupos)[index_2]),
                           paste0("n ", levels(grupos)[index_1]),
                           paste0("n ", levels(grupos)[index_2])
  )
  
  write.csv2(df, file =paste0(contraste,"_significativos.csv"))
  
}
if(T){
  
  nombres <- rownames(datos)
  
  genes <- 2^datos[nombres,]
  medias_A <- c()
  medias_B <- c()
  desviaciones_A <- c()
  desviaciones_B <- c()
  for(gen in nombres){
    medias_A <- c(medias_A,
                  genes[gen,grupos==levels(grupos)[index_1]] %>% as.numeric %>% mean)
    medias_B <- c(medias_B,
                  genes[gen,grupos==levels(grupos)[index_2]] %>% as.numeric %>% mean)
    desviaciones_A <- c(desviaciones_A,
                        genes[gen,grupos == levels(grupos)[index_1]] %>% as.numeric %>% sd)
    desviaciones_B <- c(desviaciones_B,
                        genes[gen,grupos == levels(grupos)[index_2]] %>% as.numeric %>% sd)
  }
  
  df <- data.frame(A = medias_A,
                   B = desviaciones_A,
                   C = medias_B,
                   D = desviaciones_B,
                   E = rep(plyr::count(grupos)[index_1,2],length(nombres)),
                   F = rep(plyr::count(grupos)[index_2,2],length(nombres)))
  
  rownames(df) <- nombres
  
  colnames(df) <- c(       paste0("Media ", levels(grupos))[index_1],
                           paste0("SD ", levels(grupos)[index_1]),
                           paste0("Media ", levels(grupos))[index_2],
                           paste0("SD ", levels(grupos)[index_2]),
                           paste0("n ", levels(grupos)[index_1]),
                           paste0("n ", levels(grupos)[index_2])
  )
  
  write.csv2(df, file =paste0(contraste,"_genes.csv"))
}

