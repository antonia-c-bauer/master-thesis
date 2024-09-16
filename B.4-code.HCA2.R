#Universidade de Coimbra - master thesis - Hierarchical Cluster Analysis 2


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
library(dendextend)
library(purrr)
library(cluster)
library(NbClust)
library(heatmaply)


#avoid scientific notation (promt=Eingabeaufforderung =>??)
options(prompt="R> ", scipen=999)


#set working directory (=> wd) and load data
setwd("C:/Users/Antonia/Desktop/Documents/R/master/")


#load data
variable2 <- read.csv("health outcomes data.csv")
var2_scores <- read.csv("var2_scores_CID.csv")
shp_muni <- read_sf("Cont_Mun_CAOP2023.shp")
municipio <- variable2$name


#delete the first column 
var2_scores <- var2_scores[,-1]


#delete non useful columns
shp_muni <- shp_muni[,c(-1:-2,-4:-12)]



#custom color palette
MAcol_biglight = colorRampPalette(c("#C29DE3", "#85C5FF", "#95BFC5", "#C1FFE1", "#D7D0BB", "#FFFBEF","#FFE497", "#FBC27D", "#F9A1A5") , bias = 1, space = "rgb")




#Hierarchical Clustering


#calculate distance bandwidth
scores2_dist <- dist(var2_scores) #method = default is euclidean proximity matrix, power of the Minkowski distance p=2 by default
scores2_dist


#optional: selecting the optimal clustering algorithm with the agnes function
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward") #ward in agnes() is the same as ward.D2 in hclust()

ac <- function(x) {
  agnes(var2_scores, method = x)$ac
}

map_dbl(m, ac)
#values closer to 1 suggest strong clustering structure => ward.D2 suggests strongest clustering


#selecting the optimal number of clusters
gap_stat <- clusGap(var2_scores, FUN = hcut, nstart = 25, K.max = 10, B = 500)
# Print the result
print(gap_stat, method = "firstmax")


jpeg(filename = "C:/Users/Antonia/Desktop/Documents/R/master/plots/hc2.var2.gap.jpeg", width = 25, height = 10.5, units = "cm", res = 300)
fviz_gap_stat(gap_stat, linecolor = "#00447F")
dev.off()


#hierarchical cluster analysis with hclust function
var2.hc1 <- hclust(scores2_dist, method="ward.D2") # simple dendrogram
var2.hc1$labels <- municipio


hc2.dend <- as.dendrogram(var2.hc1)

jpeg(filename = "C:/Users/Antonia/Desktop/Documents/R/master/plots/hc2.var2.dend.jpeg", width = 100, height = 42, units = "cm", res = 300)
par(mar=c(12,8,4,2))
hc2.dend %>%
  set("labels_col", value = c("#BE0C15", "#A99A6F", "#FFC425",  "#4D868E", "#00447F", "#F99D31", "#53247F",  "#00703C" ), k = 8) %>%
  set("branches_k_color", value = c("#BE0C15", "#A99A6F", "#FFC425",  "#4D868E", "#00447F", "#F99D31", "#53247F",  "#00703C"), k = 8) %>%
  set("branches_lwd", 4) %>%
  plot(axes=TRUE, cex.main = 4, cex.axis = 2.5, lwd = 2, main = "Cluster Dendrogram", xlab = "", sub = "")
title(ylab = "Height", cex.lab = 3, line = 5)
rect.dendrogram(tree = hc2.dend, k=8, border ="grey", lty = 1, lwd = 2)
grid(nx = NA, ny = NULL, lty = 1, col = "gray", lwd = 1)
legend("topright",legend = paste("Cluster",1:8), fill= c("#C29DE3", "#85C5FF", "#95BFC5", "#99C5B1FF", "#D7D0BB","#FFE497", "#FBC27D", "#F9A1A5"), bty="n", border="white", cex = 3 )
dev.off()


#cut the cluster tree into number of groups (cluster)
help("cutree")
var2.hc1.groups <- as.factor(cutree(tree = var2.hc1, k = 8))


#bind with municipality names
var2_scores <- cbind(Municipio = municipio, var2_scores)


#append cluster labels to original data
var2.hc1.data <- as.data.frame(cbind(var2_scores, cluster = as.numeric(var2.hc1.groups), deparse.level = 0))


#display first six rows of final data
head(var2.hc1.data)
write.csv(var2.hc1.data, file = "C:/Users/Antonia/Desktop/Documents/R/master/var2.hc1_data.csv")


#find mean values for each cluster
var2.hc1_mean <- aggregate(var2.hc1.data, by=list(cluster= var2.hc1.data$cluster), mean)
var2.hc1_mean <- var2.hc1_mean[,-2]
var2.hc1_mean <- var2.hc1_mean[,-8]
var2.hc1_mean


write.csv(var2.hc1_mean, file = "C:/Users/Antonia/Desktop/Documents/R/master/var2.hc1_mean.csv")


#merge with shapefile
hc1.sf <- merge(shp_muni, var2.hc1.data, by="Municipio")


#cluster map
hc2.var2.C8 <- ggplot(data = hc1.sf) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = hc1.sf, aes(fill = as.factor(cluster))) +
  scale_fill_manual(values = MAcol_biglight(8)) +
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

ggsave(filename = "hc2.var2.C8.jpg", plot = hc2.var2.C8, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




var2.hc1_mean <- round(x = var2.hc1_mean, digits = 2)


lcols <- c("#53247F", "#00447F", "#4D868E", "#00703C", "#A99A6F", "#FFC425", "#F99D31", "#BE0C15")


#plot radarplot gesamt 
hc2.rad.ges <- ggradar(var2.hc1_mean,
                       plot.title = "PC Score Mean in Clusters 1 to 8",
                       grid.min = -16,
                       grid.mid = 0,
                       grid.max = 16,
                       values.radar = c(-16, 0, 16),
                       axis.line.colour = "gray60",
                       gridline.min.colour = "gray60",
                       gridline.mid.colour = "gray60",
                       gridline.max.colour = "gray60",
                       group.colours = lcols,
                       legend.title = "Cluster",
                       legend.position = "bottom")

ggsave(filename = "hc2.rad.ges.jpg", plot = hc2.rad.ges, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var2.hc1_meanC1 <- as.data.frame(var2.hc1_mean[1,1:7])

hc2.rad.C1 <- ggradar(var2.hc1_meanC1,
                      plot.title = "PC Score Mean in Cluster 1",
                      grid.min = -2,
                      grid.mid = 0,
                      grid.max = 2,
                      values.radar = c(-2, 0, 2),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#53247F",
                      legend.title = "Cluster",
                      legend.position = "bottom")

ggsave(filename = "hc2.rad.C1.jpg", plot = hc2.rad.C1, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var2.hc1_meanC2 <- as.data.frame(var2.hc1_mean[2,1:7])

hc2.rad.C2 <- ggradar(var2.hc1_meanC2,
                      plot.title = "PC Score Mean in Cluster 2",
                      grid.min = -1,
                      grid.mid = 0,
                      grid.max = 1,
                      values.radar = c(-1, 0, 1),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#00447F")

ggsave(filename = "hc2.rad.C2.jpg", plot = hc2.rad.C2, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var2.hc1_meanC3 <- as.data.frame(var2.hc1_mean[3,1:7])

hc2.rad.C3 <- ggradar(var2.hc1_meanC3,
                      plot.title = "PC Score Mean in Cluster 3",
                      grid.min = -1,
                      grid.mid = 0,
                      grid.max = 1,
                      values.radar = c(-1, 0, 1),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#4D868E")

ggsave(filename = "hc2.rad.C3.jpg", plot = hc2.rad.C3, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var2.hc1_meanC4 <- as.data.frame(var2.hc1_mean[4,1:7])

hc2.rad.C4 <- ggradar(var2.hc1_meanC4,
                      plot.title = "PC Score Mean in Cluster 4",
                      grid.min = -3,
                      grid.mid = 0,
                      grid.max = 3,
                      values.radar = c(-3, 0, 3),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#00703C")

ggsave(filename = "hc2.rad.C4.jpg", plot = hc2.rad.C4, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var2.hc1_meanC5 <- as.data.frame(var2.hc1_mean[5,1:7])

hc2.rad.C5 <- ggradar(var2.hc1_meanC5,
                      plot.title = "PC Score Mean in Cluster 5",
                      grid.min = -7,
                      grid.mid = 0,
                      grid.max = 7,
                      values.radar = c(-7, 0, 7),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#A99A6F")

ggsave(filename = "hc2.rad.C5.jpg", plot = hc2.rad.C5, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var2.hc1_meanC6 <- as.data.frame(var2.hc1_mean[6,1:7])

hc2.rad.C6 <- ggradar(var2.hc1_meanC6,
                      plot.title = "PC Score Mean in Cluster 6",
                      grid.min = -4,
                      grid.mid = 0,
                      grid.max = 4,
                      values.radar = c(-4, 0, 4),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#FFC425")

ggsave(filename = "hc2.rad.C6.jpg", plot = hc2.rad.C6, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var2.hc1_meanC7 <- as.data.frame(var2.hc1_mean[7,1:7])

hc2.rad.C7 <- ggradar(var2.hc1_meanC7,
                      plot.title = "PC Score Mean in Cluster 7",
                      grid.min = -3,
                      grid.mid = 0,
                      grid.max = 3,
                      values.radar = c(-3, 0, 3),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#F99D31")

ggsave(filename = "hc2.rad.C7.jpg", plot = hc2.rad.C7, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var2.hc1_meanC8 <- as.data.frame(var2.hc1_mean[8,1:7])

hc2.rad.C8 <- ggradar(var2.hc1_meanC8,
                      plot.title = "PC Score Mean in Cluster 8",
                      grid.min = -16,
                      grid.mid = 0,
                      grid.max = 16,
                      values.radar = c(-16, 0, 16),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours =  "#BE0C15")

ggsave(filename = "hc2.rad.C8.jpg", plot = hc2.rad.C8, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)
