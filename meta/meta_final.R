library(metafor)
library(readxl)

patologia <- c(3,5,5,1,5,0,6,5,5,5,5,3,5,3,6,3,
               2,0,0,2,6,6,0,
               0,4,4,3,3,5,
               0,1,0,1,1,5,5,3)

contrastes_por_grupo <- c(3,1,2,1,2,3,2,1,2,2,1,3,4,1,2,1,
                          1,2,1,3,1,1,1,
                          2,1,2,1,1,2,
                          1,1,2,5,1,1,2,2)

grupos <- as.factor(rep(patologia,contrastes_por_grupo*6))
levels(grupos) <- c(0,"pulm?n","pr?stata","mama","hueso","hemato_ALL","hemato_mieloma")


tabla_completa <- read_excel("C:/Users/nacho/OneDrive/Universidad/Cuarto/TFG/Meta-analisis/estudios_completo_final.xlsx", 
                             sheet = "completo")
tabla_completa <- na.omit(tabla_completa)

tabla_completa <- as.data.frame(tabla_completa)
rownames(tabla_completa) <- tabla_completa$...1
tabla_completa <- tabla_completa[,-1]
tabla_completa <- na.omit(tabla_completa)

#pulm?n

tabla <- as.data.frame(tabla_completa)[,grupos == "pulm?n"]
index_xm2 <- seq(1,dim(tabla)[2]-5,6)
index_s2 <- seq(2,dim(tabla)[2]-4,6)
index_xm1 <- seq(3,dim(tabla)[2]-3,6)
index_s1 <- seq(4,dim(tabla)[2]-2,6)
index_n2 <- seq(5,dim(tabla)[2]-1,6)
index_n1 <- seq(6,dim(tabla)[2],6)

options("max.print"=1E9)
options("width"=10000)

res.pulmon <- NULL
for(i in 1:nrow(tabla)){
  Xm2 <- tabla[i,index_xm2]
  Xm1 <- tabla[i,index_xm1]
  S2 <- tabla[i,index_s2]
  S1 <- tabla[i,index_s1]
  n2 <- tabla[i,index_n2]
  n1 <- tabla[i,index_n1]
  ##Create matrix
  matrix_1 <- matrix(c(Xm1,Xm2,S1,S2,n1,n2),nrow=length(n1),dimnames=list(1:(dim(tabla)[2]/6),c("Xm1","Xm2","S1","S2","n1","n2")))
  datos_1 <- escalc(measure="ROM",m1i= as.numeric(Xm1), m2i= as.numeric(Xm2), sd1i= as.numeric(S1), sd2i= as.numeric(S2), n1i= as.numeric(n1), n2i= as.numeric(n2), var.names=c("LRR","LRR_var"),vtype="LS", data= matrix_1)
  res.1 <- rma.uni(measure="GEN",yi=LRR,vi=LRR_var,method="DL",data=datos_1,vtype="LS")
  res.pulmon <-rbind(res.pulmon,res.1)
}
res.pulmon <- as.data.frame(res.pulmon)
rownames(res.pulmon) <- row.names(tabla_completa)

sink("meta_pulmon.txt")
print(res.pulmon)
sink()

p_ajustados_pulmon <- p.adjust((unlist(res.pulmon$pval)), method = "BH")
sum(p_ajustados_pulmon<0.05)
rownames(tabla)[p_ajustados_pulmon<0.05]


#prostata

tabla <- as.data.frame(tabla_completa)[,grupos == "pr?stata"]
index_xm2 <- seq(1,dim(tabla)[2]-5,6)
index_s2 <- seq(2,dim(tabla)[2]-4,6)
index_xm1 <- seq(3,dim(tabla)[2]-3,6)
index_s1 <- seq(4,dim(tabla)[2]-2,6)
index_n2 <- seq(5,dim(tabla)[2]-1,6)
index_n1 <- seq(6,dim(tabla)[2],6)


options("max.print"=1E9)
options("width"=10000)

res.prostata <- NULL
for(i in 1:nrow(tabla)){
  Xm2 <- tabla[i,index_xm2]
  Xm1 <- tabla[i,index_xm1]
  S2 <- tabla[i,index_s2]
  S1 <- tabla[i,index_s1]
  n2 <- tabla[i,index_n2]
  n1 <- tabla[i,index_n1]
  ##Create matrix
  matrix_1 <- matrix(c(Xm1,Xm2,S1,S2,n1,n2),nrow=length(n1),dimnames=list(1:(dim(tabla)[2]/6),c("Xm1","Xm2","S1","S2","n1","n2")))
  datos_1 <- escalc(measure="ROM",m1i= as.numeric(Xm1), m2i= as.numeric(Xm2), sd1i= as.numeric(S1), sd2i= as.numeric(S2), n1i= as.numeric(n1), n2i= as.numeric(n2), var.names=c("LRR","LRR_var"),vtype="LS", data= matrix_1)
  res.1 <- rma.uni(measure="GEN",yi=LRR,vi=LRR_var,method="DL",data=datos_1,vtype="LS")
  res.prostata <-rbind(res.prostata,res.1)
}
res.prostata <- as.data.frame(res.prostata)
rownames(res.prostata) <- row.names(tabla_completa)

sink("meta_pr?stata.txt")
print(res.prostata)
sink()

p_ajustados_prostata <- p.adjust((unlist(res.prostata$pval)), method = "BH")
sum(p_ajustados_prostata<0.05)


#mama

tabla <- as.data.frame(tabla_completa)[,grupos == "mama"]
index_xm2 <- seq(1,dim(tabla)[2]-5,6)
index_s2 <- seq(2,dim(tabla)[2]-4,6)
index_xm1 <- seq(3,dim(tabla)[2]-3,6)
index_s1 <- seq(4,dim(tabla)[2]-2,6)
index_n2 <- seq(5,dim(tabla)[2]-1,6)
index_n1 <- seq(6,dim(tabla)[2],6)


options("max.print"=1E9)
options("width"=10000)

res.mama <- NULL
for(i in 1:nrow(tabla)){
  Xm2 <- tabla[i,index_xm2]
  Xm1 <- tabla[i,index_xm1]
  S2 <- tabla[i,index_s2]
  S1 <- tabla[i,index_s1]
  n2 <- tabla[i,index_n2]
  n1 <- tabla[i,index_n1]
  ##Create matrix
  matrix_1 <- matrix(c(Xm1,Xm2,S1,S2,n1,n2),nrow=length(n1),dimnames=list(1:(dim(tabla)[2]/6),c("Xm1","Xm2","S1","S2","n1","n2")))
  datos_1 <- escalc(measure="ROM",m1i= as.numeric(Xm1), m2i= as.numeric(Xm2), sd1i= as.numeric(S1), sd2i= as.numeric(S2), n1i= as.numeric(n1), n2i= as.numeric(n2), var.names=c("LRR","LRR_var"),vtype="LS", data= matrix_1)
  res.1 <- rma.uni(measure="GEN",yi=LRR,vi=LRR_var,method="DL",data=datos_1,vtype="LS")
  res.mama <-rbind(res.mama,res.1)
}
res.mama <- as.data.frame(res.mama)
rownames(res.mama) <- row.names(tabla)

sink("meta_mama.txt")
print(res.mama)
sink()

p_ajustados_mama <- p.adjust((unlist(res.mama$pval)), method = "BH")
sum(p_ajustados_mama<0.05)


#hueso

tabla <- as.data.frame(tabla_completa)[,grupos == "hueso"]
index_xm2 <- seq(1,dim(tabla)[2]-5,6)
index_s2 <- seq(2,dim(tabla)[2]-4,6)
index_xm1 <- seq(3,dim(tabla)[2]-3,6)
index_s1 <- seq(4,dim(tabla)[2]-2,6)
index_n2 <- seq(5,dim(tabla)[2]-1,6)
index_n1 <- seq(6,dim(tabla)[2],6)

options("max.print"=1E9)
options("width"=10000)

res.hueso <- NULL
for(i in 1:nrow(tabla)){
  Xm2 <- tabla[i,index_xm2]
  Xm1 <- tabla[i,index_xm1]
  S2 <- tabla[i,index_s2]
  S1 <- tabla[i,index_s1]
  n2 <- tabla[i,index_n2]
  n1 <- tabla[i,index_n1]
  ##Create matrix
  matrix_1 <- matrix(c(Xm1,Xm2,S1,S2,n1,n2),nrow=length(n1),dimnames=list(1:(dim(tabla)[2]/6),c("Xm1","Xm2","S1","S2","n1","n2")))
  datos_1 <- escalc(measure="ROM",m1i= as.numeric(Xm1), m2i= as.numeric(Xm2), sd1i= as.numeric(S1), sd2i= as.numeric(S2), n1i= as.numeric(n1), n2i= as.numeric(n2), var.names=c("LRR","LRR_var"),vtype="LS", data= matrix_1)
  res.1 <- rma.uni(measure="GEN",yi=LRR,vi=LRR_var,method="DL",data=datos_1,vtype="LS")
  res.hueso <-rbind(res.hueso,res.1)
}
res.hueso <- as.data.frame(res.hueso)
rownames(res.hueso) <- row.names(tabla_completa)

sink("meta_hueso.txt")
print(res.hueso)
sink()

p_ajustados_hueso <- p.adjust((unlist(res.hueso$pval)), method = "BH")
sum(p_ajustados_hueso<0.05)
rownames(tabla)[p_ajustados_hueso<0.05]

#hemato_ALL

tabla <- as.data.frame(tabla_completa)[,grupos == "hemato_ALL"]
index_xm2 <- seq(1,dim(tabla)[2]-5,6)
index_s2 <- seq(2,dim(tabla)[2]-4,6)
index_xm1 <- seq(3,dim(tabla)[2]-3,6)
index_s1 <- seq(4,dim(tabla)[2]-2,6)
index_n2 <- seq(5,dim(tabla)[2]-1,6)
index_n1 <- seq(6,dim(tabla)[2],6)


options("max.print"=1E9)
options("width"=10000)

res.hemato_ALL <- NULL
for(i in 1:nrow(tabla)){
  Xm2 <- tabla[i,index_xm2]
  Xm1 <- tabla[i,index_xm1]
  S2<- tabla[i,index_s2]
  S1<- tabla[i,index_s1]
  n2 <- tabla[i,index_n2]
  n1 <- tabla[i,index_n1]
  ##Create matrix
  matrix_1 <- matrix(c(Xm1,Xm2,S1,S2,n1,n2),nrow=length(n1),dimnames=list(1:(dim(tabla)[2]/6),c("Xm1","Xm2","S1","S2","n1","n2")))
  datos_1 <- escalc(measure="ROM",m1i= as.numeric(Xm1), m2i= as.numeric(Xm2), sd1i= as.numeric(S1), sd2i= as.numeric(S2), n1i= as.numeric(n1), n2i= as.numeric(n2), var.names=c("LRR","LRR_var"),vtype="LS", data= matrix_1)
  res.1 <- rma.uni(measure="GEN",yi=LRR,vi=LRR_var,method="DL",data=datos_1,vtype="LS")
  res.hemato_ALL <-rbind(res.hemato_ALL,res.1)
}
res.hemato_ALL <- as.data.frame(res.hemato_ALL)
rownames(res.hemato_ALL) <- row.names(tabla_completa)

sink("meta_hemato_ALL.txt")
print(res.hemato_ALL)
sink()

p_ajustados_hemato_ALL <- p.adjust((unlist(res.hemato_ALL$pval)), method = "BH")
sum(p_ajustados_hemato_ALL<0.05)

#hemato_mieloma

tabla <- as.data.frame(tabla_completa)[,grupos == "hemato_mieloma"]
index_xm2 <- seq(1,dim(tabla)[2]-5,6)
index_s2 <- seq(2,dim(tabla)[2]-4,6)
index_xm1 <- seq(3,dim(tabla)[2]-3,6)
index_s1 <- seq(4,dim(tabla)[2]-2,6)
index_n2 <- seq(5,dim(tabla)[2]-1,6)
index_n1 <- seq(6,dim(tabla)[2],6)


options("max.print"=1E9)
options("width"=10000)

res.hemato_mieloma <- NULL
for(i in 1:nrow(tabla)){
  Xm2 <- tabla[i,index_xm2]
  Xm1 <- tabla[i,index_xm1]
  S2<- tabla[i,index_s2]
  S1<- tabla[i,index_s1]
  n2 <- tabla[i,index_n2]
  n1 <- tabla[i,index_n1]
  ##Create matrix
  matrix_1 <- matrix(c(Xm1,Xm2,S1,S2,n1,n2),nrow=length(n1),dimnames=list(1:(dim(tabla)[2]/6),c("Xm1","Xm2","S1","S2","n1","n2")))
  datos_1 <- escalc(measure="ROM",m1i= as.numeric(Xm1), m2i= as.numeric(Xm2), sd1i= as.numeric(S1), sd2i= as.numeric(S2), n1i= as.numeric(n1), n2i= as.numeric(n2), var.names=c("LRR","LRR_var"),vtype="LS", data= matrix_1)
  res.1 <- rma.uni(measure="GEN",yi=LRR,vi=LRR_var,method="DL",data=datos_1,vtype="LS")
  res.hemato_mieloma <-rbind(res.hemato_mieloma,res.1)
}
res.hemato_mieloma <- as.data.frame(res.hemato_mieloma)
rownames(res.hemato_mieloma) <- row.names(tabla_completa)

sink("meta_hemato_mieloma.txt")
print(res.hemato_mieloma)
sink()

p_ajustados_hemato_mieloma <- p.adjust((unlist(res.hemato_mieloma$pval)), method = "BH")
sum(p_ajustados_hemato_mieloma<0.05)

#GLOBAL

patologia <- c(3,5,5,1,5,0,6,5,5,5,5,3,5,3,6,3,
               2,0,0,2,6,6,0,
               0,4,4,3,3,5,
               0,1,0,1,1,5,5,3)

contrastes_por_grupo <- c(3,1,2,1,2,3,2,1,2,2,1,3,4,1,2,1,
                          1,2,1,3,1,1,1,
                          2,1,2,1,1,2,
                          1,1,2,5,1,1,2,2)

grupos <- as.factor(rep(patologia,contrastes_por_grupo*6))
levels(grupos) <- c(0,"pulm?n","pr?stata","mama","hueso","hemato_ALL","hemato_mieloma")


tabla_completa <- read_excel("C:/Users/nacho/OneDrive/Universidad/Cuarto/TFG/Meta-analisis/estudios_completo_final.xlsx", 
                             sheet = "completo")
tabla_completa <- na.omit(tabla_completa)

tabla_completa <- as.data.frame(tabla_completa)
rownames(tabla_completa) <- tabla_completa$...1
tabla_completa <- tabla_completa[,-1]
tabla_completa <- na.omit(tabla_completa)

tabla <- as.data.frame(tabla_completa)
index_xm2 <- seq(1,dim(tabla)[2]-5,6)
index_s2 <- seq(2,dim(tabla)[2]-4,6)
index_xm1 <- seq(3,dim(tabla)[2]-3,6)
index_s1 <- seq(4,dim(tabla)[2]-2,6)
index_n2 <- seq(5,dim(tabla)[2]-1,6)
index_n1 <- seq(6,dim(tabla)[2],6)


options("max.print"=1E9)
options("width"=10000)


res.GLOBAL <- NULL
for(i in 1:nrow(tabla)){
  Xm2 <- tabla[i,index_xm2]
  Xm1 <- tabla[i,index_xm1]
  S2<- tabla[i,index_s2]
  S1<- tabla[i,index_s1]
  n2 <- tabla[i,index_n2]
  n1 <- tabla[i,index_n1]
  ##Create matrix
  matrix_1 <- matrix(c(Xm1,Xm2,S1,S2,n1,n2),nrow=length(n1),dimnames=list(1:(dim(tabla)[2]/6),c("Xm1","Xm2","S1","S2","n1","n2")))
  datos_1 <- escalc(measure="ROM",m1i= as.numeric(Xm1), m2i= as.numeric(Xm2), sd1i= as.numeric(S1), sd2i= as.numeric(S2), n1i= as.numeric(n1), n2i= as.numeric(n2), var.names=c("LRR","LRR_var"),vtype="LS", data= matrix_1)
  res.1 <- rma.uni(measure="GEN",yi=LRR,vi=LRR_var,method="DL",data=datos_1,vtype="LS")
  res.GLOBAL <-rbind(res.GLOBAL,res.1)
}
res.GLOBAL <- as.data.frame(res.GLOBAL)
rownames(res.GLOBAL) <- row.names(tabla_completa)

sink("meta_GLOBAL.txt")
print(res.GLOBAL)
sink()

p_ajustados_GLOBAL <- p.adjust((unlist(res.GLOBAL$pval)), method = "BH")
length(p_ajustados_GLOBAL)
sum(p_ajustados_GLOBAL<0.05)
max(p_ajustados_GLOBAL)
rownames(tabla)[p_ajustados_GLOBAL>0.05]

efecto <- unlist(res.GLOBAL$b)[p_ajustados_GLOBAL<0.05]
z <- unlist(res.GLOBAL$zval)[p_ajustados_GLOBAL<0.05]
genes_con_significacion <- data.frame(ids = rownames(tabla)[p_ajustados_GLOBAL<0.05],
                                      scores = efecto)
writexl::write_xlsx(genes_con_significacion, "genes_para_webgestalt.xlsx",
                    col_names = F)

write.table(x = genes_con_significacion, 
            file = "genes_para_webgestalt.rnk", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

# -----------------------------------------------------------

i <- 4

Xm2 <- tabla[i,index_xm2]
Xm1 <- tabla[i,index_xm1]
S2<- tabla[i,index_s2]
S1<- tabla[i,index_s1]
n2 <- tabla[i,index_n2]
n1 <- tabla[i,index_n1]

datos_G <- escalc(measure="ROM",m1i= Xm1, m2i= Xm2, sd1i= S1, sd2i= S2, n1i= n1, n2i= n2, data= datos, var.names=c("LRR","LRR_var"),vtype="LS")
res.G <- rma.uni(measure="GEN",yi=LRR,vi=LRR_var,method="DL",data=datos_G,vtype="LS",slab=paste(Author, Year, Cell.line, Dose, Time, sep=", "))

for(i in 1:length(res.G)){
  print(res.G[i])
}

tmn_efecto <- exp(res.GLOBAL$b)


df_n_genes_con_significacion <- data.frame(n_genes = c(sum(p_ajustados_pulmon<0.05),
                                                       sum(p_ajustados_prostata<0.05),
                                                       sum(p_ajustados_mama<0.05),
                                                       sum(p_ajustados_hueso<0.05),
                                                       sum(p_ajustados_hemato_mieloma<0.05),
                                                       sum(p_ajustados_hemato_ALL<0.05),
                                                       sum(p_ajustados_GLOBAL<0.05)),
                                           row.names = c("Pulm?n",
                                                         "Pr?stata",
                                                         "Mama",
                                                         "?seo",
                                                         "Hemato_mieloma",
                                                         "Hemato_ALL",
                                                         "GLOBAL"))

write.csv2(df_n_genes_con_significacion, "n_genes_significacion.csv")

