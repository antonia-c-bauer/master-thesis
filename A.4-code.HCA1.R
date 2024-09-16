# Universidade de Coimbra - master thesis - Hierarchical Cluster Analysis 1


library(base)
library(tmap)
library(GWmodel)      
library(sp)           
library(spdep)       
library(gstat)        
library(RColorBrewer) 
library(classInt)     
library(raster)       
library(gridExtra)   
library(ggplot2)      
devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
library(ggradar)
library(ggdendro)
library(dendextend)
library(dplyr)
library(sp)
library(ggspatial)
library(scales)
library(optmatch)
library(zoo)
library(reshape2)
library(factoextra)
library(FactoMineR)
library(ggrepel)
library(grDevices)
library(psych)
library(purrr)
library(cluster)
library(NbClust)
library(heatmaply)


#avoid scientific notation (promt=Eingabeaufforderung =>??)
options(prompt="R> ", scipen=999)


#set working directory (=> wd) and load data
setwd("C:/Users/Antonia/Desktop/Documents/R/master/")


#load data
variable1 <- read.csv("health determinants data.csv")
var1_scores <- read.csv("var1_scores.csv")
shp_muni <- read_sf("Cont_Mun_CAOP2023.shp")
municipio <- variable1$name


#delete the first column 
var1_scores <- var1_scores[,-1]


#delete non useful columns
shp_muni <- shp_muni[,c(-1:-2,-4:-12)]


#custom color palette
MAcol_biglight = colorRampPalette(c("#C29DE3", "#85C5FF", "#95BFC5", "#C1FFE1", "#D7D0BB", "#FFFBEF","#FFE497", "#FBC27D", "#F9A1A5") , bias = 1, space = "rgb")




#Hierarchical Clustering


#calculate distance bandwidth
scores1_dist <- dist(var1_scores) #method = default is euclidean proximity matrix, power of the Minkowski distance p=2 by default
scores1_dist


#selecting the optimal clustering algorithm with the agnes function
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward") #ward in agnes() is the same as ward.D2 in hclust()

ac <- function(x) {
  agnes(var1_scores, method = x)$ac
}

map_dbl(m, ac)
#values closer to 1 suggest strong clustering structure => ward.D2 suggests strongest clustering


#selecting the optimal number of clusters
gap_stat <- clusGap(var1_scores, FUN = hcut, nstart = 25, K.max = 10, B = 500)

#print results
print(gap_stat, method = "firstmax")


jpeg(filename = "C:/Users/Antonia/Desktop/Documents/R/master/plots/hc1.var1.gap.jpeg", width = 25, height = 10.5, units = "cm", res = 300)
fviz_gap_stat(gap_stat, linecolor = "#00447F")
dev.off()


#hierarchical cluster analysis with hclust function
var1.hc1 <- hclust(scores1_dist, method="ward.D2") #simple dendrogram #choose wards methods #ward.D2 implements Ward's (1963) clustering criterion
var1.hc1$labels <- municipio


hc1.dend <- as.dendrogram(var1.hc1)

jpeg(filename = "C:/Users/Antonia/Desktop/Documents/R/master/plots/hc1.var1.dend.jpeg", width = 100, height = 42, units = "cm", res = 300)
par(mar=c(12,8,4,2))
hc1.dend %>%
  set("labels_col", value = c("#BE0C15", "#FFC425",  "#4D868E", "#A99A6F", "#53247F"), k=5) %>%
  set("branches_k_color", value = c("#BE0C15", "#FFC425",  "#4D868E", "#A99A6F", "#53247F"), k = 5) %>%
  set("branches_lwd", 4) %>%
  plot(axes=TRUE, cex.main = 4, cex.axis = 2.5, lwd = 2, main = "Cluster Dendrogram", xlab = "", sub = "")
title(ylab = "Height", cex.lab = 3, line = 5)
rect.dendrogram(tree = hc1.dend, k=5, border ="grey", lty = 1, lwd = 2)
grid(nx = NA, ny = NULL, lty = 1, col = "gray", lwd = 1)
legend("topright",legend = paste("Cluster",1:5), fill= c("#C29DE3", "#95BFC5", "#D7D0BB","#FFE497", "#F9A1A5"), bty="n", border="white", cex = 3 )
dev.off()


#cut the cluster tree into groups
help("cutree")
var1.hc1.groups <- as.factor(cutree(tree = var1.hc1, k = 5))


#bind with municipality names
var1_scores <- cbind(Municipio = municipio, var1_scores)


#append cluster labels to original data
var1.hc1.data <- as.data.frame(cbind(var1_scores, cluster = as.numeric(var1.hc1.groups), deparse.level = 0))


#display first six rows of final data
head(var1.hc1.data)
write.csv(var1.hc1.data, file = "C:/Users/Antonia/Desktop/Documents/R/master/var1.hc1_data.csv")


#find mean values for each cluster
var1.hc1_mean <- aggregate(var1.hc1.data, by=list(cluster= var1.hc1.data$cluster), mean)
var1.hc1_mean <- var1.hc1_mean[,-2]
var1.hc1_mean <- var1.hc1_mean[,-10]
var1.hc1_mean


write.csv(var1.hc1_mean, file = "C:/Users/Antonia/Desktop/Documents/R/master/var1.hc1_mean.csv")


#merge with shapefile
hc1.sf <- merge(shp_muni, var1.hc1.data, by="Municipio")


#cluster map
hc1.var1.C5 <- ggplot(data = hc1.sf) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = hc1.sf, aes(fill = as.factor(cluster))) +
  scale_fill_manual(values = MAcol_biglight(5)) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("Cluster Map") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Cluster") +
  theme(legend.position="right", legend.key.width = unit(1, "cm")) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
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

ggsave(filename = "hc1.var1.C5.jpg", plot = hc1.var1.C5, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




var1.hc1_mean <- round(x = var1.hc1_mean, digits = 2)


lcols <- c("#53247F", "#4D868E", "#A99A6F", "#FFC425",  "#BE0C15")


#plot radarplot gesamt 
hc1.rad.ges <- ggradar(var1.hc1_mean,
                       plot.title = "PC Score Mean in Clusters 1 to 5",
                       grid.min = -120,
                       grid.mid = 0,
                       grid.max = 120,
                       values.radar = c(-120, 0, 120),
                       axis.line.colour = "gray60",
                       gridline.min.colour = "gray60",
                       gridline.mid.colour = "gray60",
                       gridline.max.colour = "gray60",
                       group.colours = lcols,
                       legend.title = "Cluster",
                       legend.position = "bottom")

ggsave(filename = "hc1.rad.ges.jpg", plot = hc1.rad.ges, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var1.hc1_meanC1 <- as.data.frame(var1.hc1_mean[1,1:9])

hc1.rad.C1 <- ggradar(var1.hc1_meanC1,
                      plot.title = "PC Score Mean in Cluster 1",
                      grid.min = -18,
                      grid.mid = 0,
                      grid.max = 18,
                      values.radar = c(-18, 0, 18),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#53247F",
                      legend.title = "Cluster",
                      legend.position = "bottom")

ggsave(filename = "hc1.rad.C1.jpg", plot = hc1.rad.C1, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var1.hc1_meanC2 <- as.data.frame(var1.hc1_mean[2,1:9])

hc1.rad.C2 <- ggradar(var1.hc1_meanC2,
                      plot.title = "PC Score Mean in Cluster 2",
                      grid.min = -12,
                      grid.mid = 0,
                      grid.max = 12,
                      values.radar = c(-12, 0, 12),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#4D868E")

ggsave(filename = "hc1.rad.C2.jpg", plot = hc1.rad.C2, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var1.hc1_meanC3 <- as.data.frame(var1.hc1_mean[3,1:9])

hc1.rad.C3 <- ggradar(var1.hc1_meanC3,
                      plot.title = "PC Score Mean in Cluster 3",
                      grid.min = -45,
                      grid.mid = 0,
                      grid.max = 45,
                      values.radar = c(-45, 0, 45),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#A99A6F")

ggsave(filename = "hc1.rad.C3.jpg", plot = hc1.rad.C3, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var1.hc1_meanC4 <- as.data.frame(var1.hc1_mean[4,1:9])

hc1.rad.C4 <- ggradar(var1.hc1_meanC4,
                      plot.title = "PC Score Mean in Cluster 4",
                      grid.min = -50,
                      grid.mid = 0,
                      grid.max = 50,
                      values.radar = c(-50, 0, 50),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#FFC425")

ggsave(filename = "hc1.rad.C4.jpg", plot = hc1.rad.C4, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var1.hc1_meanC5 <- as.data.frame(var1.hc1_mean[5,1:9])

hc1.rad.C5 <- ggradar(var1.hc1_meanC5,
                      plot.title = "PC Score Mean in Cluster 5",
                      grid.min = -116,
                      grid.mid = 0,
                      grid.max = 116,
                      values.radar = c(-116, 0, 116),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours =  "#BE0C15")

ggsave(filename = "hc1.rad.C5.jpg", plot = hc1.rad.C5, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)
