library(hash)
library(readxl)

options("max.print"=1E9)
options(digits=4)


gpl <- read_excel("C:/Users/nacho/Desktop/GPL_illumina/GPL6104-11576.xlsx")
gpl <- read_excel("C:/Users/nacho/Desktop/GPL_illumina/GPL6883-11606.xlsx")
gpl <- read_excel("C:/Users/nacho/Desktop/GPL_illumina/GPL6947-13512.xlsx")
gpl <- read_excel("C:/Users/nacho/Desktop/GPL_illumina/GPL10558-50081.xlsx")

datos <- read.table("c.txt", quote="\"", comment.char="")
gpl$`Gene stable ID`[is.na(gpl$`Gene stable ID`)] <- 0

hm <- hash(keys=gpl$ID, values=gpl$`Gene stable ID`)

rnames <- rownames(datos)
for(i in 1:nrow(datos)){
  rnames[i] <- if(is.null(hm[[rnames[i]]])){0}else{hm[[rnames[i]]]}
}
datos$gen <- rnames
datos <- datos[datos$gen!="0",]
sink("c2.txt")
print(datos)
sink()

# Por si hay sondas repetidas
c <- read.table("C:/Users/nacho/Desktop/Illumina/GSE46448_/c.txt", row.names=NULL, quote="\"", comment.char="")
n_occur <- data.frame(table(datos$row.names))
datos <- datos[!(n_occur$Freq > 1),]
rownames(datos) <- datos$row.names
datos[,-1]

genes <- datos0$gen
nombs <- data.frame(sonda = rownames(datos), gen = genes)
genes <- c()
for (sonda in rownames(PacientevsControl)) {
  genes <- c(genes, nombs$gen[nombs$sonda==sonda])
}
PacientevsControl$gen <- genes

tabla <- c()
for (nombre in base::unique(PacientevsControl$gen)){
  g <- PacientevsControl[PacientevsControl$gen==nombre,]
  PacientevsControl <- PacientevsControl[!PacientevsControl$gen==nombre,]
  g <- g[which.max(g$adj.P.Val),]
  tabla <- rbind(tabla,g)
}

PacientevsControl <- tabla
datos <- datos[rownames(PacientevsControl),]
rownames(PacientevsControl) <- PacientevsControl$gen
rownames(datos) <- PacientevsControl$gen
PacientevsControl <- subset(PacientevsControl,select = -gen)
PacientevsControl <- PacientevsControl[order(tabla$adj.P.Val),]
head(PacientevsControl)
