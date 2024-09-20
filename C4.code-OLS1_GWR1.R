#Universidade de Coimbra - master thesis - Ordinary Least Square Regression Analysis and Geographically Weighted Regression Analysis


library(ggthemes)
library(ggspatial)
library(gstat)
library(grid)
library(GWmodel)
library(remotes)
#library(rgdal) retired
#library(rgeos) retired
library(RColorBrewer)
library(reshape2)
library(sf)
library(spdep)
library(sp)
library(spgwr)
library(spData)
#install.packages('spDataLarge', repos='https://nowosad.github.io/drat/',type='source') #if not already installed
library(spDataLarge)
library(abind)
library(stars)
library(tidyverse)
library(ggplot2)
library(ggpubr)
#remotes::install_github('r-tmap/tmap')
library(tmap)
#library(tmaptools) retired
library(terra)
library(lmtest)
library(sandwich)


#change the presentation of decimal numbers to 4 and avoid scientific notation
options(prompt="R> ", digits=4, scipen=999)


#set working directory (=> wd) and load data
setwd("C:/Users/Antonia/Desktop/Documents/R/master/")


#load data
variable1 <- read.csv("health determinants data.csv")
variable2 <- read.csv("health outcomes data CID R.csv")




pca2.pc1_data <- cbind(variable2$diffwalk.pop,
                       variable2$diffhyg.pop,
                       variable2$diffmemo.pop,
                       variable2$diffunder.pop,
                       variable2$circ.pop,
                       variable2$card.pop,
                       variable2$diffhear.pop,
                       variable2$diffsee.pop,
                       variable2$resp.pop,
                       variable2$tummal.pop)


pca2.pc1_data <- scale(pca2.pc1_data)
pca2.pc1_data <- as.data.frame(pca2.pc1_data)
pca2.pc1_score <- rowSums(pca2.pc1_data)/10
pca2.pc1_score <- as.data.frame(pca2.pc1_score)


#checking on missing values in var1
pca2.pc1_score[pca2.pc1_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca2.pc1_score)) {
  pca2.pc1_score[ , i][is.na(pca2.pc1_score[ , i])] <- mean(pca2.pc1_score[ , i], na.rm = TRUE)
}




pca2.pc2_data <- cbind(variable2$suic.pop,
                       variable2$diffmemo.pop,
                       variable2$diffunder.pop,
                       variable2$tummal.pop,
                       variable2$dia.pop,
                       variable2$lesen.pop,
                       variable2$circ.pop,
                       variable2$dige.pop,
                       variable2$diffhear.pop,
                       variable2$cora.pop)


pca2.pc2_data <- scale(pca2.pc2_data)
pca2.pc2_data <- as.data.frame(pca2.pc2_data)
pca2.pc2_score <- rowSums(pca2.pc2_data)/10
pca2.pc2_score <- as.data.frame(pca2.pc2_score)


#checking on missing values in var1
pca2.pc2_score[pca2.pc2_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca2.pc2_score)) {
  pca2.pc2_score[ , i][is.na(pca2.pc2_score[ , i])] <- mean(pca2.pc2_score[ , i], na.rm = TRUE)
}




pca2.pc3_data <- cbind(variable2$suic.pop,
                       variable2$lesen.pop,
                       variable2$tumpan.pop,
                       variable2$diffhear.pop,
                       variable2$ment.pop,
                       variable2$resp.pop,
                       variable2$diffsee.pop,
                       variable2$card.pop,
                       variable2$cere.pop,
                       variable2$diffunder.pop)


pca2.pc3_data <- scale(pca2.pc3_data)
pca2.pc3_data <- as.data.frame(pca2.pc3_data)
pca2.pc3_score <- rowSums(pca2.pc3_data)/10
pca2.pc3_score <- as.data.frame(pca2.pc3_score)


#checking on missing values in var1
pca2.pc3_score[pca2.pc3_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca2.pc3_score)) {
  pca2.pc3_score[ , i][is.na(pca2.pc3_score[ , i])] <- mean(pca2.pc3_score[ , i], na.rm = TRUE)
}




pca2.pc4_data <- cbind(variable2$mengi.pop,
                       variable2$tumpan.pop,
                       variable2$dige.pop,
                       variable2$resp.pop,
                       variable2$dia.pop,
                       variable2$men.pop,
                       variable2$suic.pop,
                       variable2$tummal.pop,
                       variable2$cora.pop,
                       variable2$circ.pop)


pca2.pc4_data <- scale(pca2.pc4_data)
pca2.pc4_data <- as.data.frame(pca2.pc4_data)
pca2.pc4_score <- rowSums(pca2.pc4_data)/10
pca2.pc4_score <- as.data.frame(pca2.pc4_score)


#checking on missing values in var1
pca2.pc4_score[pca2.pc4_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca2.pc4_score)) {
  pca2.pc4_score[ , i][is.na(pca2.pc4_score[ , i])] <- mean(pca2.pc4_score[ , i], na.rm = TRUE)
}




pca2.pc5_data <- cbind(variable2$tumpan.pop,
                       variable2$ment.pop,
                       variable2$nerv.pop,
                       variable2$dia.pop,
                       variable2$mengi.pop,
                       variable2$suic.pop,
                       variable2$resp.pop,
                       variable2$diffhear.pop,
                       variable2$diffmemo.pop,
                       variable2$diffunder.pop)


pca2.pc5_data <- scale(pca2.pc5_data)
pca2.pc5_data <- as.data.frame(pca2.pc5_data)
pca2.pc5_score <- rowSums(pca2.pc5_data)/10
pca2.pc5_score <- as.data.frame(pca2.pc5_score)


#checking on missing values in var1
pca2.pc5_score[pca2.pc5_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca2.pc5_score)) {
  pca2.pc5_score[ , i][is.na(pca2.pc5_score[ , i])] <- mean(pca2.pc5_score[ , i], na.rm = TRUE)
}




pca2.pc6_data <- cbind(variable2$men.pop,
                       variable2$mengi.pop,
                       variable2$cora.pop,
                       variable2$diffsee.pop,
                       variable2$suic.pop,
                       variable2$diffhear.pop,
                       variable2$resp.pop,
                       variable2$cere.pop,
                       variable2$tumpan.pop,
                       variable2$card.pop)


pca2.pc6_data <- scale(pca2.pc6_data)
pca2.pc6_data <- as.data.frame(pca2.pc6_data)
pca2.pc6_score <- rowSums(pca2.pc6_data)/10
pca2.pc6_score <- as.data.frame(pca2.pc6_score)


#checking on missing values in var1
pca2.pc6_score[pca2.pc6_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca2.pc6_score)) {
  pca2.pc6_score[ , i][is.na(pca2.pc6_score[ , i])] <- mean(pca2.pc6_score[ , i], na.rm = TRUE)
}




pca1.pc1_data <- cbind(variable1$X1_apa.bui,
                       variable1$rent.acc,
                       variable1$p_pop_edu2,
                       variable1$mom.fam,
                       variable1$lo.acc,
                       variable1$prop.acc,
                       variable1$p_fam_mono,
                       variable1$X4_apa.bui,
                       variable1$i_pp_pc,
                       variable1$fp.acc)


pca1.pc1_data <- scale(pca1.pc1_data)
pca1.pc1_data <- as.data.frame(pca1.pc1_data)
pca1.pc1_score <- rowSums(pca1.pc1_data)/10
pca1.pc1_score <- as.data.frame(pca1.pc1_score)


#checking on missing values in var1
pca1.pc1_score[pca1.pc1_score == ""] <- NA


#replace NA in all columns 
for(i in 1:ncol(pca1.pc1_score)) {
  pca1.pc1_score[ , i][is.na(pca1.pc1_score[ , i])] <- mean(pca1.pc1_score[ , i], na.rm = TRUE)
}




pca1.pc2_data <- cbind(variable1$park.acc,
                       variable1$nopark.acc,
                       variable1$X50.59m.acc,
                       variable1$p_com_imt,
                       variable1$X60.79m.acc,
                       variable1$C.empl,
                       variable1$X40.49m.acc,
                       variable1$X200.x.acc,
                       variable1$X120.149m.acc,
                       variable1$p_com_ped)


pca1.pc2_data <- scale(pca1.pc2_data)
pca1.pc2_data <- as.data.frame(pca1.pc2_data)
pca1.pc2_score <- rowSums(pca1.pc2_data)/10
pca1.pc2_score <- as.data.frame(pca1.pc2_score)


#checking on missing values in var1
pca1.pc2_score[pca1.pc2_score == ""] <- NA


#replace NA in all columns => https://statisticsglobe.com/replace-missing-values-by-column-mean-in-r
for(i in 1:ncol(pca1.pc2_score)) {
  pca1.pc2_score[ , i][is.na(pca1.pc2_score[ , i])] <- mean(pca1.pc2_score[ , i], na.rm = TRUE)
}




pca1.pc3_data <- cbind(variable1$a_floor.bui,
                       variable1$n_acc.km2,
                       variable1$a_pop.km2,
                       variable1$nomar.fam,
                       variable1$arti.ter,
                       variable1$nonrel.pop,
                       variable1$cat.rel,
                       variable1$rel.pop,
                       variable1$X10.x_apa.bui,
                       variable1$a_acc.bui)


pca1.pc3_data <- scale(pca1.pc3_data)
pca1.pc3_data <- as.data.frame(pca1.pc3_data)
pca1.pc3_score <- rowSums(pca1.pc3_data)/10
pca1.pc3_score <- as.data.frame(pca1.pc3_score)


#checking on missing values in var1
pca1.pc3_score[pca1.pc3_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca1.pc3_score)) {
  pca1.pc3_score[ , i][is.na(pca1.pc3_score[ , i])] <- mean(pca1.pc3_score[ , i], na.rm = TRUE)
}




pca1.pc4_data <- cbind(variable1$X10.15mc.acc,
                       variable1$X80.xmc.acc,
                       variable1$X20.30mc.acc,
                       variable1$X15.20mc.acc,
                       variable1$X60.80mc.acc,
                       variable1$X30.40mc.acc,
                       variable1$X40.49m.acc,
                       variable1$X50.59m.acc,
                       variable1$a_rents,
                       variable1$nomar.fam)


pca1.pc4_data <- scale(pca1.pc4_data)
pca1.pc4_data <- as.data.frame(pca1.pc4_data)
pca1.pc4_score <- rowSums(pca1.pc4_data)/10
pca1.pc4_score <- as.data.frame(pca1.pc4_score)


#checking on missing values in var1
pca1.pc4_score[pca1.pc4_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca1.pc4_score)) {
  pca1.pc4_score[ , i][is.na(pca1.pc4_score[ , i])] <- mean(pca1.pc4_score[ , i], na.rm = TRUE)
}




pca1.pc5_data <- cbind(variable1$L.empl,
                       variable1$R.empl,
                       variable1$mh.acc,
                       variable1$p_com_ped,
                       variable1$ocr.rel,
                       variable1$p_com_imt,
                       variable1$O.empl,
                       variable1$a_rents,
                       variable1$noh.acc,
                       variable1$pro.rel)


pca1.pc5_data <- scale(pca1.pc5_data)
pca1.pc5_data <- as.data.frame(pca1.pc5_data)
pca1.pc5_score <- rowSums(pca1.pc5_data)/10
pca1.pc5_score <- as.data.frame(pca1.pc5_score)


#checking on missing values in var1
pca1.pc5_score[pca1.pc5_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca1.pc5_score)) {
  pca1.pc5_score[ , i][is.na(pca1.pc5_score[ , i])] <- mean(pca1.pc5_score[ , i], na.rm = TRUE)
}




pca1.pc6_data <- cbind(variable1$X2_apa.bui,
                       variable1$prop.acc,
                       variable1$rent.acc,
                       variable1$X3_apa.bui,
                       variable1$mh.acc,
                       variable1$inter.pop,
                       variable1$C.empl,
                       variable1$K.empl,
                       variable1$p_com_imt,
                       variable1$J.empl)


pca1.pc6_data <- scale(pca1.pc6_data)
pca1.pc6_data <- as.data.frame(pca1.pc6_data)
pca1.pc6_score <- rowSums(pca1.pc6_data)/10
pca1.pc6_score <- as.data.frame(pca1.pc6_score)


#checking on missing values in var1
pca1.pc6_score[pca1.pc6_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca1.pc6_score)) {
  pca1.pc6_score[ , i][is.na(pca1.pc6_score[ , i])] <- mean(pca1.pc6_score[ , i], na.rm = TRUE)
}




pca1.pc7_data <- cbind(variable1$p_com_ped,
                       variable1$a_ear_mont,
                       variable1$p_com_imt,
                       variable1$X40.49m.acc,
                       variable1$X50.59m.acc,
                       variable1$X60.79m.acc,
                       variable1$res1.acc,
                       variable1$J.empl,
                       variable1$R.empl,
                       variable1$C.empl)


pca1.pc7_data <- scale(pca1.pc7_data)
pca1.pc7_data <- as.data.frame(pca1.pc7_data)
pca1.pc7_score <- rowSums(pca1.pc7_data)/10
pca1.pc7_score <- as.data.frame(pca1.pc7_score)


#checking on missing values in var1
pca1.pc7_score[pca1.pc7_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca1.pc7_score)) {
  pca1.pc7_score[ , i][is.na(pca1.pc7_score[ , i])] <- mean(pca1.pc7_score[ , i], na.rm = TRUE)
}




pca1.pc8_data <- cbind(variable1$C.empl,
                       variable1$exres.bui,
                       variable1$prres.bui,
                       variable1$p_pop_edu3,
                       variable1$X200.x.acc,
                       variable1$X30.40mc.acc,
                       variable1$n_acc.km2,
                       variable1$J.empl,
                       variable1$p_pop_edu2,
                       variable1$nomar.fam)


pca1.pc8_data <- scale(pca1.pc8_data)
pca1.pc8_data <- as.data.frame(pca1.pc8_data)
pca1.pc8_score <- rowSums(pca1.pc8_data)/10
pca1.pc8_score <- as.data.frame(pca1.pc8_score)


#checking on missing values in var1
pca1.pc8_score[pca1.pc8_score == ""] <- NA


#replace NA in all columns
for(i in 1:ncol(pca1.pc8_score)) {
  pca1.pc8_score[ , i][is.na(pca1.pc8_score[ , i])] <- mean(pca1.pc8_score[ , i], na.rm = TRUE)
}




municipio <- variable1$name

var3.hc1_data <- read.csv("score.clust.csv")
var3.hc1_data <- var3.hc1_data[,-1:-2]
var3.hc1_data <- as.data.frame(cbind(Municipio = municipio, var3.hc1_data,
                                     pca1.pc1_score,
                                     pca1.pc2_score,
                                     pca1.pc3_score,
                                     pca1.pc4_score,
                                     pca1.pc5_score,
                                     pca1.pc6_score,
                                     pca1.pc7_score,
                                     pca1.pc8_score,
                                     pca2.pc1_score,
                                     pca2.pc2_score,
                                     pca2.pc3_score,
                                     pca2.pc4_score,
                                     pca2.pc5_score,
                                     pca2.pc6_score))




var.independent <- as.data.frame(cbind(pca1.pc1_score,
                                       pca1.pc2_score,
                                       pca1.pc3_score,
                                       pca1.pc4_score,
                                       pca1.pc5_score,
                                       pca1.pc6_score,
                                       pca1.pc7_score,
                                       pca1.pc8_score,
                                       pca2.pc1_score,
                                       pca2.pc2_score,
                                       pca2.pc3_score,
                                       pca2.pc4_score,
                                       pca2.pc5_score,
                                       pca2.pc6_score))




#CORRELATION MATRIX 1 

var1_cormat <- round(cor(var.independent),2)
help(round)

#get lower triangle of the correlation matrix
get_lower_tri1 <- function(var1_cormat){
  var1_cormat[upper.tri(var1_cormat)] <- NA
  return(var1_cormat)
}

#get upper triangle of the correlation matrix
get_upper_tri1 <- function(var1_cormat){
  var1_cormat[lower.tri(var1_cormat)]<- NA
  return(var1_cormat)
}

# Use correlation between variables as distance
reorder_cormat <- function(var1_cormat){
  dd <- as.dist((1-var1_cormat)/2)
  hc <- hclust(dd)
  var1_cormat <-var1_cormat[hc$order, hc$order]
}

#reorder the correlation matrix
var1_cormat <- reorder_cormat(var1_cormat)
upper_tri1 <- get_upper_tri1(var1_cormat)

#melt the correlation matrix
melt_cormat1 <- melt(upper_tri1, na.rm = TRUE)

# Create a ggheatmap
var3.cormat <- 
  ggplot(melt_cormat1, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson Correlation Coefficient") +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
  ggtitle("Correlation Matrix") +
  theme(text = element_text(size = 16)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey", size = 0.25, linetype = 1),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0.5, 0),
    legend.position = "bottom",
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 0.5, size = 10, hjust = 1))+
  coord_fixed()


print(var3.cormat)

ggsave(filename = "var3.cormat.jpg", plot = var3.cormat, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 26, units = "cm", dpi = 300)

#write.csv(var3.cormat, "C:/Users/Antonia/Desktop/Documents/R/master/var3.cormat.csv", row.names=TRUE)




#import shapefile
shp_muni <- read_sf("Cont_Mun_CAOP2023.shp")

#delete non useful columns
shp_muni <- shp_muni[,c(-1:-2,-4:-12)]




#Linear Regression Analysis and Geographically Weighted Regression Analysis with all Principal Components per Cluster




# Regression Analysis Cluster 1


#get the scores for cluster 1
var3.data.C1 <- var3.hc1_data %>% group_by(cluster) %>% dplyr::filter(cluster == 1)


#Linear Regression Analysis Cluster 1
linear_model1 <- lm(pca2.pc1_score ~ 
                    pca1.pc1_score +
                    pca1.pc2_score +
                    pca1.pc3_score +
                    pca1.pc4_score +
                    pca1.pc5_score +
                    pca1.pc6_score +
                    pca1.pc7_score +
                    pca1.pc8_score +
                    
                    pca2.pc2_score +
                    pca2.pc3_score +
                    pca2.pc4_score +
                    pca2.pc5_score +
                    pca2.pc6_score,
                    data = var3.data.C1)


summary(linear_model1)


plot(fitted(linear_model1), resid(linear_model1))


#t-test for robust standard errors
coeftest(linear_model1, vcov = vcovHC(linear_model1, type = 'HC0'))


par(mfrow=c(2,2))
plot(linear_model1)


plot(linear_model1, which=1)
grid()
abline(coef = coef(linear_model1))




#name residuals and bind them to the existing data frame
residuals <- as.data.frame(cbind(Municipio = var3.data.C1$Municipio, residuals = residuals(linear_model1)))
residuals$residuals <- as.numeric(residuals$residuals)
residuals.sf  <- merge(shp_muni, residuals,  by="Municipio")
head(residuals.sf)




#map linear_model1 option with gradient fill, plot as A4 portrait pdf
lm.C1 <- ggplot() + 
  geom_sf(data = shp_muni, fill = "grey", color = "lightgrey") +
  geom_sf(data = residuals.sf, aes(fill= residuals)) +
  scale_fill_gradient2(low= "khaki4", mid = "#ffffff", high = "#006666", midpoint=0,
                       breaks = c(-1,0,1)) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("Residual Map Cluster 1") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Value") +
  theme(legend.position="bottom", legend.key.width = unit(2.5, "cm")) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey", "white"),
    height = unit(0.5, "cm"),
    text_cex = 1.5) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(
      fill = c("grey", "grey"),
      line_col = "grey20"),
    height = unit(2, "cm"),
    width = unit(2, "cm")) 

print(lm.C1)

ggsave(filename = "lm.C1.jpg", plot = lm.C1, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




shp_muni.data  <- merge(shp_muni, var3.data.C1,  by="Municipio")


#calculate kernel bandwidth
centro <- st_centroid(residuals.sf) #get the centroid for each polygone
coords.centro <- st_coordinates(centro) #in order to retrieve coordinates in matrix form (lon,lat)




#Geographically Weighted Regression Analysis Cluster 1
gwrbandwidth <- gwr.sel(pca2.pc1_score ~ 
                          pca1.pc1_score +
                          pca1.pc2_score +
                          pca1.pc3_score +
                          pca1.pc4_score +
                          pca1.pc5_score +
                          pca1.pc6_score +
                          pca1.pc7_score +
                          pca1.pc8_score +
                          
                          pca2.pc2_score +
                          pca2.pc3_score +
                          pca2.pc4_score +
                          pca2.pc5_score +
                          pca2.pc6_score,
                        data = shp_muni.data,
                        coords = coords.centro,
                        adapt = TRUE)




gwr.model = gwr(pca2.pc1_score ~ 
                  pca1.pc1_score +
                  pca1.pc2_score +
                  pca1.pc3_score +
                  pca1.pc4_score +
                  pca1.pc5_score +
                  pca1.pc6_score +
                  pca1.pc7_score +
                  pca1.pc8_score +
                  
                  pca2.pc2_score +
                  pca2.pc3_score +
                  pca2.pc4_score +
                  pca2.pc5_score +
                  pca2.pc6_score,
                data = shp_muni.data,
                coords = coords.centro,
                adapt = gwrbandwidth,
                hatmatrix = TRUE,
                se.fit = TRUE)


print(gwr.model)
gwr.model$bandwidth

gwr.results <- as.data.frame(gwr.model$SDF)
gwr.results <- cbind(Municipio = var3.data.C1$Municipio, gwr.results)
summary(gwr.results)


#just for visualisation for interpreting the results in the master thesis (not relevant for code)
gwr.results.localR2 <- cbind(Municipality = gwr.results$Municipio, localR2 = gwr.results$localR2)
################################################################################################


gwr.data <- merge(shp_muni.data, gwr.results, by="Municipio")




lm.C1.R2 <- ggplot() + 
  geom_sf(data = shp_muni, fill = "grey", color = "lightgrey") +
  geom_sf(data = gwr.data, aes(fill= `localR2`)) +
  scale_fill_gradient(low = "#ffff00", high = "#006666") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("R2 Map Cluster 1") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Value") +
  theme(legend.position="bottom", legend.key.width = unit(2.5, "cm")) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey", "white"),
    height = unit(0.5, "cm"),
    text_cex = 1.5) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(
      fill = c("grey", "grey"),
      line_col = "grey20"),
    height = unit(2, "cm"),
    width = unit(2, "cm")) 


print(lm.C1.R2)


ggsave(filename = "lm.C1.R2.jpg", plot = lm.C1.R2, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




# Regression Analysis Cluster 2


#get the scores for cluster 2
var3.data.C2 <- var3.hc1_data %>% group_by(cluster) %>% dplyr::filter(cluster == 2)


#Linear Regression Analysis Cluster 2
linear_model1 <- lm(pca2.pc5_score ~ 
                      pca1.pc1_score +
                      pca1.pc2_score +
                      pca1.pc3_score +
                      pca1.pc4_score +
                      pca1.pc5_score +
                      pca1.pc6_score +
                      pca1.pc7_score +
                      pca1.pc8_score +
                      
                      pca2.pc1_score +
                      pca2.pc2_score +
                      pca2.pc3_score +
                      pca2.pc4_score +
                      
                      pca2.pc6_score,
                    data = var3.data.C2)

summary(linear_model1)

plot(fitted(linear_model1), resid(linear_model1))

coeftest(linear_model1, vcov = vcovHC(linear_model1, type = 'HC0'))


par(mfrow=c(2,2))
plot(linear_model1)


plot(linear_model1, which=1)
grid()
abline(coef = coef(linear_model1))




#name residuals and bind them to the existing data frame
residuals <- as.data.frame(cbind(Municipio = var3.data.C2$Municipio, residuals = residuals(linear_model1)))
residuals$residuals <- as.numeric(residuals$residuals)
residuals.sf  <- merge(shp_muni, residuals,  by="Municipio")
head(residuals.sf)




#map linear_model1 option with gradient fill, plot as A4 portrait pdf
lm.C2 <- ggplot() + 
  geom_sf(data = shp_muni, fill = "grey", color = "lightgrey") +
  geom_sf(data = residuals.sf, aes(fill= residuals)) +
  scale_fill_gradient2(low= "khaki4", mid = "#ffffff", high = "#006666", midpoint=0,
                       breaks = c(-1,0,1)) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("Residual Map Cluster 2") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Value") +
  theme(legend.position="bottom", legend.key.width = unit(2.5, "cm")) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey", "white"),
    height = unit(0.5, "cm"),
    text_cex = 1.5) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(
      fill = c("grey", "grey"),
      line_col = "grey20"),
    height = unit(2, "cm"),
    width = unit(2, "cm")) 


ggsave(filename = "lm.C2.jpg", plot = lm.C2, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




shp_muni.data  <- merge(shp_muni, var3.data.C2,  by="Municipio")


#calculate kernel bandwidth for model2
centro <- st_centroid(residuals.sf) #get the centroid for each polygone
coords.centro <- st_coordinates(centro) #in order to retrieve coordinates in matrix form (lon,lat)




#Geographically Weighted Regression Analysis Cluster 2
gwrbandwidth <- gwr.sel(pca2.pc5_score ~ 
                          pca1.pc1_score +
                          pca1.pc2_score +
                          pca1.pc3_score +
                          pca1.pc4_score +
                          pca1.pc5_score +
                          pca1.pc6_score +
                          pca1.pc7_score +
                          pca1.pc8_score +
                          
                          pca2.pc1_score +
                          pca2.pc2_score +
                          pca2.pc3_score +
                          pca2.pc4_score +
                          
                          pca2.pc6_score,
                        data = shp_muni.data,
                        coords = coords.centro,
                        adapt = T)




gwr.model = gwr(pca2.pc5_score ~ 
                  pca1.pc1_score +
                  pca1.pc2_score +
                  pca1.pc3_score +
                  pca1.pc4_score +
                  pca1.pc5_score +
                  pca1.pc6_score +
                  pca1.pc7_score +
                  pca1.pc8_score +
                  
                  pca2.pc1_score +
                  pca2.pc2_score +
                  pca2.pc3_score +
                  pca2.pc4_score +
                  
                  pca2.pc6_score,
                data = shp_muni.data,
                coords = coords.centro,
                adapt = gwrbandwidth,
                hatmatrix = TRUE,
                se.fit = TRUE)


print(gwr.model)
gwr.model$bandwidth

gwr.results <- as.data.frame(gwr.model$SDF)
gwr.results <- cbind(Municipio = var3.data.C2$Municipio, gwr.results)
summary(gwr.results)


#just for visualisation for interpreting the results in the master thesis (not relevant for code)
gwr.results.localR2 <- cbind(Municipality = gwr.results$Municipio, localR2 = gwr.results$localR2)
################################################################################################


gwr.data <- merge(shp_muni.data, gwr.results, by="Municipio")




lm.C2.R2 <- ggplot() + 
  geom_sf(data = shp_muni, fill = "grey", color = "lightgrey") +
  geom_sf(data = gwr.data, aes(fill= `localR2`)) +
  scale_fill_gradient(low = "#ffff00", high = "#006666") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("R2 Map Cluster 2") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Value") +
  theme(legend.position="bottom", legend.key.width = unit(2.5, "cm")) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey", "white"),
    height = unit(0.5, "cm"),
    text_cex = 1.5) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(
      fill = c("grey", "grey"),
      line_col = "grey20"),
    height = unit(2, "cm"),
    width = unit(2, "cm")) 


ggsave(filename = "lm.C2.R2.jpg", plot = lm.C2.R2, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




# Regression Analysis Cluster 3


#get the scores for cluster 3
var3.data.C3 <- var3.hc1_data %>% group_by(cluster) %>% dplyr::filter(cluster == 3)


#Linear Regression Analysis Cluster 3
linear_model1 <- lm(pca2.pc4_score ~ 
                    pca1.pc1_score +
                    pca1.pc2_score +
                    pca1.pc3_score +
                    pca1.pc4_score +
                    pca1.pc5_score +
                    pca1.pc6_score +
                    pca1.pc7_score +
                    pca1.pc8_score +
                      
                    pca2.pc1_score +
                    pca2.pc2_score +
                    pca2.pc3_score +
                    
                    pca2.pc5_score +
                    pca2.pc6_score, 
                    data = var3.data.C3)

summary(linear_model1)

plot(fitted(linear_model1), resid(linear_model1))

coeftest(linear_model1, vcov = vcovHC(linear_model1, type = 'HC0'))



#NO RESULTS FOR linear_model1 DUE TO SMALL NUMBER OF OBSERVATIONS




# Regression Analysis Cluster 4


#get the scores for cluster 4
var3.data.C4 <- var3.hc1_data %>% group_by(cluster) %>% dplyr::filter(cluster == 4)


#Linear Regression Analysis Cluster 4
linear_model1 <- lm(pca2.pc2_score ~ 
                      pca1.pc1_score +
                      pca1.pc2_score +
                      pca1.pc3_score +
                      pca1.pc4_score +
                      pca1.pc5_score +
                      pca1.pc6_score +
                      pca1.pc7_score +
                      pca1.pc8_score +
                      
                      pca2.pc1_score +
                      
                      pca2.pc3_score +
                      pca2.pc4_score +
                      pca2.pc5_score +
                      pca2.pc6_score,
                    data = var3.data.C4)

summary(linear_model1)

plot(fitted(linear_model1), resid(linear_model1))

coeftest(linear_model1, vcov = vcovHC(linear_model1, type = 'HC0'))


par(mfrow=c(2,2))
plot(linear_model1)


plot(linear_model1, which=1)
grid()
abline(coef = coef(linear_model1))




#name residuals and bind them to the existing data frame
residuals <- as.data.frame(cbind(Municipio = var3.data.C4$Municipio, residuals = residuals(linear_model1)))
residuals$residuals <- as.numeric(residuals$residuals)
residuals.sf  <- merge(shp_muni, residuals,  by="Municipio")
head(residuals.sf)




#map linear_model1 option with gradient fill, plot as A4 portrait pdf
lm.C4 <- ggplot() + 
  geom_sf(data = shp_muni, fill = "grey", color = "lightgrey") +
  geom_sf(data = residuals.sf, aes(fill= residuals)) +
  scale_fill_gradient2(low= "khaki4", mid = "#ffffff", high = "#006666", midpoint=0,
                       breaks = c(-1,0,1)) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("Residual Map Cluster 4") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Value") +
  theme(legend.position="bottom", legend.key.width = unit(2.5, "cm")) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey", "white"),
    height = unit(0.5, "cm"),
    text_cex = 1.5) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(
      fill = c("grey", "grey"),
      line_col = "grey20"),
    height = unit(2, "cm"),
    width = unit(2, "cm")) 

print(lm.C4)

ggsave(filename = "lm.C4.jpg", plot = lm.C4, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




shp_muni.data  <- merge(shp_muni, var3.data.C4,  by="Municipio")


#calculate kernel bandwidth for model2
centro <- st_centroid(residuals.sf) #get the centroid for each polygone
coords.centro <- st_coordinates(centro) #in order to retrieve coordinates in matrix form (lon,lat)




#Geographically Weighted Regression Analysis Cluster 4
gwrbandwidth <- gwr.sel(pca2.pc2_score ~ 
                          pca1.pc1_score +
                          pca1.pc2_score +
                          pca1.pc3_score +
                          pca1.pc4_score +
                          pca1.pc5_score +
                          pca1.pc6_score +
                          pca1.pc7_score +
                          pca1.pc8_score +
                          
                          pca2.pc1_score +
                          
                          pca2.pc3_score +
                          pca2.pc4_score +
                          pca2.pc5_score +
                          pca2.pc6_score ,
                        data = shp_muni.data,
                        coords = coords.centro,
                        adapt = T)




gwr.model = gwr(pca2.pc2_score ~ 
                  pca1.pc1_score +
                  pca1.pc2_score +
                  pca1.pc3_score +
                  pca1.pc4_score +
                  pca1.pc5_score +
                  pca1.pc6_score +
                  pca1.pc7_score +
                  pca1.pc8_score +
                  
                  pca2.pc1_score +
                  
                  pca2.pc3_score +
                  pca2.pc4_score +
                  pca2.pc5_score +
                  pca2.pc6_score ,
                data = shp_muni.data,
                coords = coords.centro,
                adapt = gwrbandwidth,
                hatmatrix = TRUE,
                se.fit = TRUE)


print(gwr.model)
gwr.model$bandwidth

gwr.results <- as.data.frame(gwr.model$SDF)
gwr.results <- cbind(Municipio = var3.data.C4$Municipio, gwr.results)
summary(gwr.results)


#just for visualisation for interpreting the results in the master thesis (not relevant for code)
gwr.results.localR2 <- cbind(Municipality = gwr.results$Municipio, localR2 = gwr.results$localR2)
################################################################################################


gwr.data <- merge(shp_muni.data, gwr.results, by="Municipio")




lm.C4.R2 <- ggplot() + 
  geom_sf(data = shp_muni, fill = "grey", color = "lightgrey") +
  geom_sf(data = gwr.data, aes(fill= `localR2`)) +
  scale_fill_gradient(low = "#ffff00", high = "#006666") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("R2 Map Cluster 4") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Value") +
  theme(legend.position="bottom", legend.key.width = unit(2.5, "cm")) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey", "white"),
    height = unit(0.5, "cm"),
    text_cex = 1.5) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(
      fill = c("grey", "grey"),
      line_col = "grey20"),
    height = unit(2, "cm"),
    width = unit(2, "cm")) 


print(lm.C4.R2)


ggsave(filename = "lm.C4.R2.jpg", plot = lm.C4.R2, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)





# Regression Analysis Cluster 5


#get the scores per cluster
var3.data.C5 <- var3.hc1_data %>% group_by(cluster) %>% dplyr::filter(cluster == 5)


#Linear Regression Analysis Cluster 5
linear_model1 <- lm(pca2.pc2_score ~ 
                      pca1.pc1_score +
                      pca1.pc2_score +
                      pca1.pc3_score +
                      pca1.pc4_score +
                      pca1.pc5_score +
                      pca1.pc6_score +
                      pca1.pc7_score +
                      pca1.pc8_score +
                      
                      pca2.pc1_score +
                      
                      pca2.pc3_score +
                      pca2.pc4_score +
                      pca2.pc5_score +
                      pca2.pc6_score,
                    data = var3.data.C5)

summary(linear_model1)

plot(fitted(linear_model1), resid(linear_model1))

coeftest(linear_model1, vcov = vcovHC(linear_model1, type = 'HC0'))


par(mfrow=c(2,2))
plot(linear_model1)


plot(linear_model1, which=1)
grid()
abline(coef = coef(linear_model1))




#name residuals and bind them to the existing data frame
residuals <- as.data.frame(cbind(Municipio = var3.data.C5$Municipio, residuals = residuals(linear_model1)))
residuals$residuals <- as.numeric(residuals$residuals)
residuals.sf  <- merge(shp_muni, residuals,  by="Municipio")
head(residuals.sf)




#map linear_model1 option with gradient fill, plot as A4 portrait pdf
lm.C5 <- ggplot() + 
  geom_sf(data = shp_muni, fill = "grey", color = "lightgrey") +
  geom_sf(data = residuals.sf, aes(fill= residuals)) +
  scale_fill_gradient2(low= "khaki4", mid = "#ffffff", high = "#006666", midpoint=0,
                       breaks = c(-1,0,1)) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("Residual Map Cluster 5") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Value") +
  theme(legend.position="bottom", legend.key.width = unit(2.5, "cm")) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey", "white"),
    height = unit(0.5, "cm"),
    text_cex = 1.5) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(
      fill = c("grey", "grey"),
      line_col = "grey20"),
    height = unit(2, "cm"),
    width = unit(2, "cm")) 

print(lm.C5)

ggsave(filename = "lm.C5.jpg", plot = lm.C5, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




shp_muni.data  <- merge(shp_muni, var3.data.C5,  by="Municipio")


#calculate kernel bandwidth for model2
centro <- st_centroid(residuals.sf) #get the centroid for each polygone
coords.centro <- st_coordinates(centro) #in order to retrieve coordinates in matrix form (lon,lat)




#Geographically Weighted Reggression Analysis Cluster 5
gwrbandwidth <- gwr.sel(pca2.pc2_score ~ 
                           pca1.pc1_score +
                           pca1.pc2_score +
                           pca1.pc3_score +
                           pca1.pc4_score +
                           pca1.pc5_score +
                           pca1.pc6_score +
                           pca1.pc7_score +
                           pca1.pc8_score +
                           
                           pca2.pc1_score +
                           
                           pca2.pc3_score +
                           pca2.pc4_score +
                           pca2.pc5_score +
                           pca2.pc6_score,
                        data = shp_muni.data,
                        coords = coords.centro,
                        adapt = T)




gwr.model = gwr(pca2.pc2_score ~ 
                   pca1.pc1_score +
                   pca1.pc2_score +
                   pca1.pc3_score +
                   pca1.pc4_score +
                   pca1.pc5_score +
                   pca1.pc6_score +
                   pca1.pc7_score +
                   pca1.pc8_score +
                   
                   pca2.pc1_score +
                   
                   pca2.pc3_score +
                   pca2.pc4_score +
                   pca2.pc5_score +
                   pca2.pc6_score,
                 data = shp_muni.data,
                coords = coords.centro,
                adapt = gwrbandwidth,
                hatmatrix = TRUE,
                se.fit = TRUE)


print(gwr.model)
gwr.model$bandwidth

gwr.results <- as.data.frame(gwr.model$SDF)
gwr.results <- cbind(Municipio = var3.data.C5$Municipio, gwr.results)
summary(gwr.results)



#just for visualisation for interpreting the results in the master thesis (not relevant for code)
gwr.results.localR2 <- cbind(Municipality = gwr.results$Municipio, localR2 = gwr.results$localR2)
################################################################################################


#gwr.data <- cbind(residuals.sf, as.matrix(gwr.results))
gwr.data <- merge(shp_muni.data, gwr.results, by="Municipio")




lm.C5.R2 <- ggplot() + 
  geom_sf(data = shp_muni, fill = "grey", color = "lightgrey") +
  geom_sf(data = gwr.data, aes(fill= `localR2`)) +
  scale_fill_gradient(low = "#ffff00", high = "#006666") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("R2 Map Cluster 5") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Value") +
  theme(legend.position="bottom", legend.key.width = unit(2.5, "cm")) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  ggspatial::annotation_scale(
    location = "br",
    bar_cols = c("grey", "white"),
    height = unit(0.5, "cm"),
    text_cex = 1.5) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering(
      fill = c("grey", "grey"),
      line_col = "grey20"),
    height = unit(2, "cm"),
    width = unit(2, "cm")) 


print(lm.C5.R2)


ggsave(filename = "lm.C5.R2.jpg", plot = lm.C5.R2, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




# Regression Analysis Cluster 5


#get the scores per cluster
var3.data.C6 <- var3.hc1_data %>% group_by(cluster) %>% dplyr::filter(cluster == 6)

#no results, cluster is only one observation



