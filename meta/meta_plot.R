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

Author <- c("Wu","Real","Rainer","Muzikar","Carlet","Carlet","Rickles","Samon","Chen","Chen","Aneichyk","Lauriola","Bindreither","West","Hern?ndez-Garc?a","Moore",
            "Sahu","Paakinaho","Backman","Morhayim","Burwick","Thomas","Proven?al",
            "Nehm?","Jewell","Galliher-Beckley","Burd","Kittler","Diaz-Jimenez",
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

tabla[i,]

index_xm2 <- seq(1,dim(tabla)[2]-5,6)
index_s2 <- seq(2,dim(tabla)[2]-4,6)
index_xm1 <- seq(3,dim(tabla)[2]-3,6)
index_s1 <- seq(4,dim(tabla)[2]-2,6)
index_n2 <- seq(5,dim(tabla)[2]-1,6)
index_n1 <- seq(6,dim(tabla)[2],6)


which((p_ajustados_pulmon>0.05)==T)
i <- 2350
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
levels(tejido) <- c(0,"pulm?n","pr?stata","mama","hueso","hemato_ALL","hemato_mieloma")

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
#write.csv2(df,"df.csv")

m.cont <- metacont(n.e = n1,
                   mean.e = Xm1,
                   sd.e = S1,
                   n.c = n2,
                   mean.c = Xm2,
                   sd.c = S2,
                   studlab = Author,
                   data = df[df$Tejido=="pulm?n",],
                   sm = "ROM",
                   fixed = F,
                   random = T,
                   method.tau = "DL",
                   hakn = F,
                   title = "minimini")

#m.cont <- update.meta(m.cont, prediction = TRUE)

dev.new(width=9, height=7.5)
par(mar=c(4,4,1,2))
plot.new()
forest.meta(m.cont, 
            sortvar = TE,
            print.tau2 = FALSE,
            leftcols = c("studlab","Year" ,"GEO_ID"),
            leftlabs = c("Autor","A?o", "GEO_ID"),
            rightlabs = c("FC", "95% IC", "Peso"),
            text.random = "Global",
            text.predict = "intervalo",
            test.overall.random = T,
            label.test.overall.random = "Modelo efectos aleatorios: ",
            print.Q = T,
            print.pval.Q = T,
            print.tau = T,
            smlab = "Fold Change (t/c)",
            label.e = "e",
            label.c = "c",
            hetlab = "Heterogeneidad: ")
grid.text("Gen GFI1 y dexametasona (pulm?n)", .355, .95, gp=gpar(cex=1))


# Funnel plot
funnel.meta(m.cont,
            xlim = c(0.004, 4),
            ylim = c(5,0),
            studlab = TRUE)
metabias(m.cont, method.bias = "linreg")

tf <- trimfill(m.cont)

funnel.meta(tf,
            xlim = c(0.005, 2.5),
            ylim = c(3,0),
            studlab = TRUE)


metabias(m.cont, method.bias = "linreg")
y = m.cont$TE/m.cont$seTE
x = (1/m.cont$seTE) 
model <- lm(y ~ x, data = m.cont$data)


a <- data.frame(TE = 2^m.cont$TE, seTE = m.cont$seTE)
viz_funnel(a)

# Define fill colors for contour
col.contour = c("gray75", "gray85", "gray95")
col.contour <-  RColorBrewer::brewer.pal(3, "Blues")

# Generate funnel plot (we do not include study labels here)
funnel.meta(m.cont, xlim = c(0.01, 2),
            contour = c(0.9, 0.95, 0.99),
            col.contour = col.contour,
            xlab = "Ratio of Means",
            ylab = "Error estandar")

# Legend
legend(x = 0.8, y = 0.01, 
       legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
       fill = col.contour)

# Title
title("Contour-Enhanced Funnel Plot")


# --------------------- pr?stata -------------------------------
GEO_ID <- c("GSE4917", "GSE7067", "GSE11336","GSE17307","GSE22152", "GSE22779", "GSE30644", "GSE33562", "GSE39338", "GSE39339", "GSE48680", "GSE53405", "GSE55878", "GSE79761", "GSE94341", "GSE139870",
            "GSE39654", "GSE48328", "GSE54608", "GSE56172", "GSE59080", "GSE64704", "GSE119842",
            "GSE16643", "GSE26857", "GSE28912", "GSE30592", "GSE42619", "GSE135130",
            "GSE52778", "GSE79432", "GSE80651", "GSE104714", "GSE104908", "GSE117796", "GSE137893", "GSE152201")

Author <- c("Wu","Real","Rainer","Muzikar","Carlet","Carlet","Rickles","Samon","Chen","Chen","Aneichyk","Lauriola","Bindreither","West","Hern?ndez-Garc?a","Moore",
            "Sahu","Paakinaho","Backman","Morhayim","Burwick","Thomas","Proven?al",
            "Nehm?","Jewell","Galliher-Beckley","Burd","Kittler","Diaz-Jimenez",
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


which((p_ajustados_pulmon>0.05)==T)
i <- 2350
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
levels(tejido) <- c(0,"pulm?n","pr?stata","mama","hueso","hemato_ALL","hemato_mieloma")

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

#m.cont <- update.meta(m.cont, prediction = TRUE)

dev.new(width=9, height=7.5)
par(mar=c(4,4,1,2))
plot.new()
forest.meta(m.cont, 
            sortvar = data$Year,
            print.tau2 = FALSE,
            leftcols = c("studlab","Year" ,"GEO_ID"),
            leftlabs = c("Autor","A?o", "GEO_ID"),
            rightlabs = c("FC", "95% IC", "Peso"),
            text.random = "Global",
            text.predict = "intervalo",
            test.overall.random = T,
            label.test.overall.random = " ",
            print.Q = T,
            print.pval.Q = T,
            print.tau = T,
            smlab = "Fold Change (t/c)",
            label.e = "e",
            label.c = "c",
            hetlab = "")
grid.text("Gen GFI1 y dexametasona (pulm?n)", .355, .95, gp=gpar(cex=1))





# --------------------- mama -------------------------------
GEO_ID <- c("GSE4917", "GSE7067", "GSE11336","GSE17307","GSE22152", "GSE22779", "GSE30644", "GSE33562", "GSE39338", "GSE39339", "GSE48680", "GSE53405", "GSE55878", "GSE79761", "GSE94341", "GSE139870",
            "GSE39654", "GSE48328", "GSE54608", "GSE56172", "GSE59080", "GSE64704", "GSE119842",
            "GSE16643", "GSE26857", "GSE28912", "GSE30592", "GSE42619", "GSE135130",
            "GSE52778", "GSE79432", "GSE80651", "GSE104714", "GSE104908", "GSE117796", "GSE137893", "GSE152201")

Author <- c("Wu","Real","Rainer","Muzikar","Carlet","Carlet","Rickles","Samon","Chen","Chen","Aneichyk","Lauriola","Bindreither","West","Hern?ndez-Garc?a","Moore",
            "Sahu","Paakinaho","Backman","Morhayim","Burwick","Thomas","Proven?al",
            "Nehm?","Jewell","Galliher-Beckley","Burd","Kittler","Diaz-Jimenez",
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
index_xm1 <- seq(2,dim(tabla)[2]-4,6)
index_s2<- seq(3,dim(tabla)[2]-3,6)
index_s1 <- seq(4,dim(tabla)[2]-2,6)
index_n2 <- seq(5,dim(tabla)[2]-1,6)
index_n1 <- seq(6,dim(tabla)[2],6)


which((p_ajustados_pulmon>0.05)==T)
i <- 2350
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
levels(tejido) <- c(0,"pulm?n","pr?stata","mama","hueso","hemato_ALL","hemato_mieloma")

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
#write.csv2(df,"df.csv")




m.cont <- metacont(n.e = n1,
                   mean.e = Xm1,
                   sd.e = S1,
                   n.c = n2,
                   mean.c = Xm2,
                   sd.c = S2,
                   studlab = Author,
                   data = df[df$Tejido=="mama",],
                   sm = "ROM",
                   fixed = F,
                   random = T,
                   method.tau = "DL",
                   hakn = F,
                   title = "meta mama")

#m.cont <- update.meta(m.cont, prediction = TRUE)


forest.meta(m.cont, 
            sortvar = TE,
            print.tau2 = FALSE,
            leftcols = c("studlab","Year" ,"GEO_ID"),
            leftlabs = c("Autor","A?o", "GEO_ID"),
            rightlabs = c("FC", "95% IC", "Peso"),
            text.random = "Global",
            text.predict = "intervalo",
            test.overall.random = T,
            label.test.overall.random = "Modelo efectos aleatorios: ",
            print.Q = T,
            print.pval.Q = T,
            print.tau = T,
            smlab = "Fold Change (t/c)",
            label.e = "e",
            label.c = "c",
            hetlab = "Heterogeneidad: ")
grid.text("Gen GFI1 y dexametasona (pulm?n)", .355, .95, gp=gpar(cex=1))




# --------------------- hemato_mieloma -------------------------------
GEO_ID <- c("GSE4917", "GSE7067", "GSE11336","GSE17307","GSE22152", "GSE22779", "GSE30644", "GSE33562", "GSE39338", "GSE39339", "GSE48680", "GSE53405", "GSE55878", "GSE79761", "GSE94341", "GSE139870",
            "GSE39654", "GSE48328", "GSE54608", "GSE56172", "GSE59080", "GSE64704", "GSE119842",
            "GSE16643", "GSE26857", "GSE28912", "GSE30592", "GSE42619", "GSE135130",
            "GSE52778", "GSE79432", "GSE80651", "GSE104714", "GSE104908", "GSE117796", "GSE137893", "GSE152201")

Author <- c("Wu","Real","Rainer","Muzikar","Carlet","Carlet","Rickles","Samon","Chen","Chen","Aneichyk","Lauriola","Bindreither","West","Hern?ndez-Garc?a","Moore",
            "Sahu","Paakinaho","Backman","Morhayim","Burwick","Thomas","Proven?al",
            "Nehm?","Jewell","Galliher-Beckley","Burd","Kittler","Diaz-Jimenez",
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


which((p_ajustados_pulmon>0.05)==T)
i <- 2350
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
levels(tejido) <- c(0,"pulm?n","pr?stata","mama","hueso","hemato_ALL","hemato_mieloma")

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

m.cont <- metacont(n.e = n1,
                   mean.e = Xm1,
                   sd.e = S1,
                   n.c = n2,
                   mean.c = Xm2,
                   sd.c = S2,
                   studlab = Author,
                   data = df[df$Tejido=="hemato_mieloma",],
                   sm = "ROM",
                   fixed = F,
                   random = T,
                   method.tau = "DL",
                   hakn = F,
                   title = "meta mama")

#m.cont <- update.meta(m.cont, prediction = TRUE)


forest.meta(m.cont, 
            sortvar = data$Year,
            print.tau2 = FALSE,
            leftcols = c("studlab","Year" ,"GEO_ID"),
            leftlabs = c("Autor","A?o", "GEO_ID"),
            rightlabs = c("FC", "95% IC", "Peso"),
            text.random = "Global",
            text.predict = "intervalo",
            test.overall.random = T,
            label.test.overall.random = "Modelo efectos aleatorios: ",
            print.Q = T,
            print.pval.Q = T,
            print.tau = T,
            smlab = "Fold Change (t/c)",
            label.e = "e",
            label.c = "c",
            hetlab = "Heterogeneidad: ")
grid.text("Gen GFI1 y dexametasona (pulm?n)", .355, .95, gp=gpar(cex=1))




# --------------------- hemato_ALL -------------------------------
GEO_ID <- c("GSE4917", "GSE7067", "GSE11336","GSE17307","GSE22152", "GSE22779", "GSE30644", "GSE33562", "GSE39338", "GSE39339", "GSE48680", "GSE53405", "GSE55878", "GSE79761", "GSE94341", "GSE139870",
            "GSE39654", "GSE48328", "GSE54608", "GSE56172", "GSE59080", "GSE64704", "GSE119842",
            "GSE16643", "GSE26857", "GSE28912", "GSE30592", "GSE42619", "GSE135130",
            "GSE52778", "GSE79432", "GSE80651", "GSE104714", "GSE104908", "GSE117796", "GSE137893", "GSE152201")

Author <- c("Wu","Real","Rainer","Muzikar","Carlet","Carlet","Rickles","Samon","Chen","Chen","Aneichyk","Lauriola","Bindreither","West","Hern?ndez-Garc?a","Moore",
            "Sahu","Paakinaho","Backman","Morhayim","Burwick","Thomas","Proven?al",
            "Nehm?","Jewell","Galliher-Beckley","Burd","Kittler","Diaz-Jimenez",
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


which((p_ajustados_pulmon>0.05)==T)
i <- 2350
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
levels(tejido) <- c(0,"pulm?n","pr?stata","mama","hueso","hemato_ALL","hemato_mieloma")

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
                   data = df[df$Tejido=="hemato_ALL",],
                   sm = "ROM",
                   fixed = F,
                   random = T,
                   method.tau = "DL",
                   hakn = F,
                   title = "meta mama")

#m.cont <- update.meta(m.cont, prediction = TRUE)


forest.meta(m.cont, 
            sortvar = data$Year,
            print.tau2 = FALSE,
            leftcols = c("studlab","Year" ,"GEO_ID"),
            leftlabs = c("Autor","A?o", "GEO_ID"),
            rightlabs = c("FC", "95% IC", "Peso"),
            text.random = "Global",
            text.predict = "intervalo",
            test.overall.random = T,
            label.test.overall.random = "",
            print.Q = T,
            print.pval.Q = T,
            print.tau = T,
            smlab = "Fold Change (t/c)",
            label.e = "e",
            label.c = "c",
            hetlab = "")
grid.text("Gen GFI1 y dexametasona (pulm?n)", .355, .95, gp=gpar(cex=1))


# -------------------------------------------------------
# --------------------- GLOBAL -------------------------------
GEO_ID <- c("GSE4917", "GSE7067", "GSE11336","GSE17307","GSE22152", "GSE22779", "GSE30644", "GSE33562", "GSE39338", "GSE39339", "GSE48680", "GSE53405", "GSE55878", "GSE79761", "GSE94341", "GSE139870",
            "GSE39654", "GSE48328", "GSE54608", "GSE56172", "GSE59080", "GSE64704", "GSE119842",
            "GSE16643", "GSE26857", "GSE28912", "GSE30592", "GSE42619", "GSE135130",
            "GSE52778", "GSE79432", "GSE80651", "GSE104714", "GSE104908", "GSE117796", "GSE137893", "GSE152201")

Author <- c("Wu","Real","Rainer","Muzikar","Carlet","Carlet","Rickles","Samon","Chen","Chen","Aneichyk","Lauriola","Bindreither","West","Hern?ndez-Garc?a","Moore",
            "Sahu","Paakinaho","Backman","Morhayim","Burwick","Thomas","Proven?al",
            "Nehm?","Jewell","Galliher-Beckley","Burd","Kittler","Diaz-Jimenez",
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


which(rownames(tabla)=="ENSG00000171848") # gen importante en GSEA
i <- 2857
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
levels(tejido) <- c(0,"pulm?n","pr?stata","mama","hueso","hemato_ALL","hemato_mieloma")

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

#m.cont <- update.meta(m.cont, prediction = TRUE)


forest.meta(m.cont, 
            sortvar = data$Year,
            print.tau2 = FALSE,
            leftcols = c("studlab","Year" ,"GEO_ID"),
            leftlabs = c("Autor","A?o", "GEO_ID"),
            rightlabs = c("FC", "95% IC", "Peso"),
            text.random = "Global",
            text.predict = "intervalo",
            test.overall.random = T,
            label.test.overall.random = "",
            print.Q = T,
            print.pval.Q = T,
            print.tau = T,
            smlab = "Fold Change (t/c)",
            label.e = "e",
            label.c = "c",
            hetlab = "")
grid.text("Gen GFI1 y dexametasona (pulm?n)", .355, .95, gp=gpar(cex=1))


# ------------------------ sesgo de publicaci?n ---------------------------------
col.contour = c("gray75", "gray85", "gray95")

# Generate funnel plot (we do not include study labels here)
funnel.meta(m.cont, xlim = c(0.5, 1.5),
            ylim = c(0.5,0),
            contour = c(0.9, 0.95, 0.99),
            col.contour = col.contour,
            xlab = "Ratio of Means",
            ylab = "Error est?ndar")

# Add a legend
legend(x = 1.25, y = 0.01, 
       legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
       fill = col.contour)

# Add a title
title("Contour-Enhanced Funnel Plot para el metaan?lisis global del gen RRM2")


# egger
metabias(m.cont, method.bias = "linreg")
y = m.cont$TE/m.cont$seTE
x = (1/m.cont$seTE) 
model <- lm(y ~ x, data = m.cont$data)

df<-data.frame(x,y)
library(ggplot2)
ggplot(df,aes(x,y))+
  geom_point()+stat_smooth(method="lm",se=F)+annotate("text",x=3,y=1,label=(paste0("Intercepto==",coef(lm(df$y~df$x))[1])),parse=TRUE) +
  theme_minimal() + xlab("1/seTE") + ylab("TE/seTE")
tf <- trimfill(m.cont)

# Analyze with outliers removed
tf.no.out <- trimfill(update(m.cont, 
                             subset = -c(3, 16)))



# Define fill colors for contour
contour <- c(0.9, 0.95, 0.99)
col.contour <- c("gray75", "gray85", "gray95")
ld <- c("p < 0.1", "p < 0.05", "p < 0.01")

# Use 'par' to create two plots in one row (row, columns)
par(mfrow=c(1,2))

# Contour-enhanced funnel plot (full data)
funnel.meta(m.cont, 
            xlim = c(0.005, 2), contour = contour,
            col.contour = col.contour,
            ylim = c(4.5,0))
legend(x = 1, y = 0.01, 
       legend = ld, fill = col.contour)
title("Funnel Plot")

# Contour-enhanced funnel plot (outliers removed)
funnel.meta(tf, 
            xlim = c(0.005, 2), contour = contour,
            col.contour = col.contour,
            ylim = c(4.5,0))
legend(x = 1, y = 0.01, 
       legend = ld, fill = col.contour)
title("Funnel Plot (Trim and Fill)")




funnel.meta(m.cont, xlim = c(0.01, 2),
            ylim = c(1.4,0),
            contour = c(0.9, 0.95, 0.99),
            col.contour = col.contour,
            xlab = "Ratio of Means",
            ylab = "Error est?ndar")

