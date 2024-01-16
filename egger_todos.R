library(tidyverse)
library(readxl)
library(meta)
library(metafor)
library(metaviz)
library(grid)


GEO_ID <- c("GSE4917", "GSE7067", "GSE11336","GSE17307","GSE22152", "GSE22779", "GSE30644", "GSE33562", "GSE39338", "GSE39339", "GSE48680", "GSE53405", "GSE55878", "GSE79761", "GSE94341", "GSE139870",
            "GSE39654", "GSE48328", "GSE54608", "GSE56172", "GSE59080", "GSE64704", "GSE119842",
            "GSE16643", "GSE26857", "GSE28912", "GSE30592", "GSE42619", "GSE135130",
            "GSE52778", "GSE79432", "GSE80651", "GSE104714", "GSE104908", "GSE117796", "GSE137893", "GSE152201")

Author <- c("Wu","Real","Rainer","Muzikar","Carlet","Carlet","Rickles","Samon","Chen","Chen","Aneichyk","Lauriola","Bindreither","West","Hernández-García","Moore",
            "Sahu","Paakinaho","Backman","Morhayim","Burwick","Thomas","Provençal",
            "Nehmé","Jewell","Galliher-Beckley","Burd","Kittler","Diaz-Jimenez",
            "Himes","Vockley","Jiang","McDOwell","Li","Poulard","Meyer","Conway")

Year <- c("2006","2009","2009","2009","2010","2010","2012","2012","2013","2013","2013","2014","2014","2016","2017","2020",
          "2013","2013","2014","2015","2014","2015","2019",
          "2009","2011","2011","2012","2013","2020",
          "2014","2016","2016","2018","2019","2018","2019","2020")

contrastes_por_grupo <- c(3,1,2,1,2,3,2,1,2,2,1,3,4,1,2,1,
                          1,2,1,3,1,1,1,
                          2,1,2,1,1,2,
                          1,1,2,5,1,1,2,2)

GEO_ID <- as.factor(rep(GEO_ID, contrastes_por_grupo))
Author <- as.factor(rep(Author, contrastes_por_grupo))
Year <- as.factor(rep(Year, contrastes_por_grupo))

tabla_completa <- read_excel("C:/Users/nacho/OneDrive/Universidad/Cuarto/TFG/Meta-analisis/estudios_completo_final.xlsx", 
                             sheet = "completo")
tabla <- as.data.frame(tabla_completa)
rownames(tabla) <- tabla$...1
tabla <- tabla[,-1]
tabla <- na.omit(tabla)

index_xm2 <- seq(1,dim(tabla)[2]-5,6)
index_s2 <- seq(2,dim(tabla)[2]-4,6)
index_xm1 <- seq(3,dim(tabla)[2]-3,6)
index_s1 <- seq(4,dim(tabla)[2]-2,6)
index_n2 <- seq(5,dim(tabla)[2]-1,6)
index_n1 <- seq(6,dim(tabla)[2],6)


pvalores <- c()


for(i in 1:nrow(tabla)){
  
  Xm2 <- t(tabla)[index_xm2,i]
  Xm1 <- t(tabla)[index_xm1,i]
  S2<- t(tabla)[index_s2,i]
  S1<- t(tabla)[index_s1,i]
  n2 <- t(tabla)[index_n2,i]
  n1 <- t(tabla)[index_n1,i]
  
  
  patologia <- c(3,5,5,1,5,0,6,5,5,5,5,3,5,3,6,3,
                 2,0,0,2,6,6,0,
                 0,4,4,3,3,5,
                 0,1,0,1,1,5,5,3)
  contrastes_por_grupo <- c(3,1,2,1,2,3,2,1,2,2,1,3,4,1,2,1,
                            1,2,1,3,1,1,1,
                            2,1,2,1,1,2,
                            1,1,2,5,1,1,2,2)
  tejido <- as.factor(rep(patologia,contrastes_por_grupo))
  levels(tejido) <- c(0,"pulmón","próstata","mama","hueso","hemato_ALL","hemato_mieloma")
  
  df <- data.frame(Author = Author,
                   Year = Year,
                   GEO_ID = GEO_ID,
                   Tejido = tejido,
                   Xm1 = Xm1,
                   S1 = S1,
                   Xm2 = Xm2,
                   S2 = S2,
                   n1 = n1,
                   n2 = n2)
  df <- df[order(df$Year),]
  #write.csv2(df,"df.csv")
  
  
  
  m.cont <- metacont(n.e = n1,
                     mean.e = Xm1,
                     sd.e = S1,
                     n.c = n2,
                     mean.c = Xm2,
                     sd.c = S2,
                     studlab = Author,
                     data = df,
                     sm = "ROM",
                     fixed = F,
                     random = T,
                     method.tau = "DL",
                     hakn = F,
                     title = "meta")
  pvalores <- c(pvalores, metabias(m.cont, method.bias = "linreg")$pval)
}

