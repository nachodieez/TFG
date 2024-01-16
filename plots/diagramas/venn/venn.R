library(ggvenn)
library(venn)
x <- list(
  pulmón = rownames(tabla)[p_ajustados_pulmon<0.05],
  próstata = rownames(tabla)[p_ajustados_prostata<0.05],
  mama = rownames(tabla)[p_ajustados_mama<0.05],
  óseo = rownames(tabla)[p_ajustados_hueso<0.05],
  hemato_mieloma = rownames(tabla)[p_ajustados_hemato_mieloma<0.05],
  hemato_ALL = rownames(tabla)[p_ajustados_hemato_ALL<0.05]
)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
draw.quintuple.venn(pulmón, próstata, mama, óseo, hemato_mieloma)

pulmón = rownames(tabla)[p_ajustados_pulmon<0.05]
próstata = rownames(tabla)[p_ajustados_prostata<0.05]
mama = rownames(tabla)[p_ajustados_mama<0.05]
óseo = rownames(tabla)[p_ajustados_hueso<0.05]
hemato_mieloma = rownames(tabla)[p_ajustados_hemato_mieloma<0.05]
hemato_ALL = rownames(tabla)[p_ajustados_hemato_ALL<0.05]

write.csv2(rownames(tabla)[p_ajustados_hemato_ALL<0.05],"hemato_ALL")
venn(x)
