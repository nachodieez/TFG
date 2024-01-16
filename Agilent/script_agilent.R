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
head(datos0)
datos00 <- df
datos0 <- datos00[,-c(1,ncol(datos00))]
grupos0 <- as.factor(c(rep(c("Control_WTGR", "Dex_WTGR"),each = 3), rep(c("Control_AGR", "Dex_AGR"),each = 3)))

individuos <- data.frame(individuo = colnames(datos0), grupo = grupos0)
individuos

datos <- datos0[, grupos0 == "Control_WTGR" | grupos0 == "Dex_WTGR"]
grupos <- droplevels(grupos0[grupos0 == "Control_WTGR" | grupos0 == "Dex_WTGR"])

design <- model.matrix(~0+grupos)
colnames(design) <- levels(grupos)

contr.matrix <- makeContrasts(
  ControlvsDex = Control - Dex,
  levels = grupos
)

fit <- lmFit(datos, design)
fit <- contrasts.fit(fit, contrasts = contr.matrix)
summary(decideTests(fit))
ControlvsDex <- topTreat(eBayes(fit), n=Inf)
genes <- genes$Gene.stable.ID
nombs <- data.frame(sonda = rownames(datos), gen = genes)
genes <- c()
for (sonda in rownames(ControlvsDex)) {
  genes <- c(genes, nombs$gen[nombs$sonda==sonda])
}
ControlvsDex$gen <- genes

tabla <- c()
for (nombre in base::unique(ControlvsDex$gen)){
  g <- ControlvsDex[ControlvsDex$gen==nombre,]
  ControlvsDex <- ControlvsDex[!ControlvsDex$gen==nombre,]
  g <- g[which.min(g$adj.P.Val),]
  tabla <- rbind(tabla,g)
}

ControlvsDex <- tabla
datos <- datos[rownames(ControlvsDex),]
rownames(ControlvsDex) <- ControlvsDex$gen
rownames(datos) <- ControlvsDex$gen
ControlvsDex <- subset(ControlvsDex,select = -gen)
ControlvsDex <- ControlvsDex[order(tabla$adj.P.Val),]
head(ControlvsDex)
contraste <- deparse(substitute(ControlvsDex))


head(ControlvsDex)

nombres <- rownames(ControlvsDex[ControlvsDex$adj.P.Val<0.05 ,])
length(nombres)

#write.table(file = paste0(contraste,"_significativos.txt"), x = nombres, sep = "\n", row.names = F, col.names = F)

index_1 <- 1 #index del primer grupo del contraste el levels(groups)
index_2 <- 2 #index del segundo grupo del contraste
if(length(nombres)!=0){
  genes <- datos[nombres,]
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
  
  genes <- datos[nombres,]
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
