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

options("max.print"=1E9)

#codigo <- ""
datos0 <- read.csv(paste0(codigo,".txt"), sep="")
head(datos0)

grupos0 <- as.factor()

individuos <- data.frame(individuo = colnames(datos0), grupo = grupos0)
individuos

pca <- prcomp(t(datos), scale = T, rank = 2)
summary(pca)

plot(pca, col = "blue", main = "Varianza explicada por cada una de las componentes")

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
contraste <- deparse(substitute(AvsB))

png(file = paste0("MDplot_",contraste,".png"))
plotMD(fit, column = 1, status = decideTests(fit),
       main = paste0("MD Plot ",codigo, " (",contraste, ")"))
dev.off()

png(file = paste0("PCA_",contraste,".png"))
autoplot(pca, data = data.frame(t(datos), grupos = grupos), colour = 'grupos') + 
  ggtitle(paste0("PCA de las muestras del experimento ",codigo, " (",contraste,")"))
dev.off()

sink(paste0(contraste,"_topTreat.txt"))
print(AvsB)
sink()


png(file = paste0("heatmap_",contraste,".png"))
pheatmap::pheatmap(as.matrix(datos[rownames(AvsB)[1:30],]))
dev.off()
dev.off()
png(file = paste0("heatmap_",contraste,".png"))
pheatmap::pheatmap(as.matrix(datos[rownames(AvsB)[1:30],]))
dev.off()

head(AvsB)

nombres <- rownames(AvsB[AvsB$adj.P.Val<0.05 ,])
length(nombres)

write.table(file = paste0(contraste,"_significativos.txt"), x = nombres, sep = "\n", row.names = F, col.names = F)

if(length(nombres)!=0){
  genes <- 2^datos[nombres,]
  medias_A <- c()
  medias_B <- c()
  desviaciones_A <- c()
  desviaciones_B <- c()
  for(gen in nombres){
    medias_A <- c(medias_A,
                  genes[gen,grupos==levels(grupos)[1]] %>% as.numeric %>% mean)
    medias_B <- c(medias_B,
                  genes[gen,grupos==levels(grupos)[2]] %>% as.numeric %>% mean)
    desviaciones_A <- c(desviaciones_A,
                        genes[gen,grupos == levels(grupos)[1]] %>% as.numeric %>% sd)
    desviaciones_B <- c(desviaciones_B,
                        genes[gen,grupos == levels(grupos)[2]] %>% as.numeric %>% sd)
  }
  
  df <- data.frame(A = medias_A,
                   B = medias_B,
                   C = desviaciones_A,
                   D = desviaciones_B,
                   E = rep(plyr::count(grupos)[1,2],length(nombres)),
                   F = rep(plyr::count(grupos)[2,2],length(nombres)))
  
  rownames(df) <- nombres
  
  colnames(df) <- c(       paste0("Media ", levels(grupos))[1],
                           paste0("Media ", levels(grupos))[2],
                           paste0("SD ", levels(grupos)[1]),
                           paste0("SD ", levels(grupos)[2]),
                           paste0("n ", levels(grupos)[1]),
                           paste0("n ", levels(grupos)[2])
  )
  
  sink(paste0(contraste,"_significativos_MediasYDesviaciones.txt"))
  print(df)
  sink()
}
