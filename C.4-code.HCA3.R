#Universidade de Coimbra - master thesis - Hierarchical Cluster Analysis 3 (with previous scaling of the scores)


library(base)
library(fastmap)
library(usethis)
library(devtools)
library(tmap)

library(robustbase)
library(sp)           
library(Rcpp)
library(GWmodel)      

library(spData)
library(sf)
library(spdep)        

library(gstat)        
library(RColorBrewer) 
library(classInt)     
library(raster)       
library(gridExtra)    
library(ggplot2)     

library(ggradar)
library(dplyr)

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

library(plotly)
library(viridis)
library(viridisLite)
library(heatmaply)


#avoid scientific notation (promt=Eingabeaufforderung =>??)
options(prompt="R> ", scipen=999)


#set working directory (=> wd) and load data
setwd("C:/Users/Antonia/Desktop/Documents/R/master/")


#load data
variable1 <- read.csv("health determinants data.csv")
var1_scores <- read.csv("var1_scores.csv")
var2_scores <- read.csv("var2_scores_CID.csv")
shp_muni <- read_sf("Cont_Mun_CAOP2023.shp")
municipio <- variable1$name


#delete the first column
var1_scores <- var1_scores[,-1]


#data preparation var1
colnames(var1_scores) <- paste(c("GPCA1"),c(colnames(var1_scores)), sep = ".")


#delete first column
var2_scores <- var2_scores[,-1]


#data preparation var2
colnames(var2_scores) <- paste(c("GPCA2"),c(colnames(var2_scores)), sep = ".")


#bind the two score data frames
var3_scores.original <- cbind(var1_scores, var2_scores)
var3_scores <- scale(var3_scores.original)
var3_scores <- as.data.frame(var3_scores)


#delete non useful columns
shp_muni <- shp_muni[,c(-1:-2,-4:-12)]


#continuous colors for discrete scales
MAcol_biglight = colorRampPalette(c("#C29DE3", "#85C5FF", "#95BFC5", "#C1FFE1", "#D7D0BB", "#FFFBEF","#FFE497", "#FBC27D", "#F9A1A5") , bias = 1, space = "rgb")




#Hierarchical Clustering


#calculate distance bandwidth
scores3_dist <- dist(var3_scores) #method = default is euclidean proximity matrix, power of the Minkowski distance p=2 by default
scores3_dist


#optional: selecting the optimal clustering algorithm with the agnes function
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward") #ward in agnes() is the same as ward.D2 in hclust()

ac <- function(x) {
  agnes(var3_scores, method = x)$ac
}

map_dbl(m, ac)
#values closer to 1 suggest strong clustering structure => ward.D2 suggests strongest clustering


#selecting the optimal number of clusters
gap_stat <- clusGap(var3_scores, FUN = hcut, nstart = 25, K.max = 10, B = 500)
# Print the result
print(gap_stat, method = "firstmax")


jpeg(filename = "C:/Users/Antonia/Desktop/Documents/R/master/plots/hc3.var3.gap.jpeg", width = 25, height = 10.5, units = "cm", res = 300)
fviz_gap_stat(gap_stat, linecolor = "#00447F")
dev.off()




#hierarchical cluster analysis with hclust function
var3.hc1 <- hclust(scores3_dist, method="ward.D2") # simple dendrogram
var3.hc1$labels <- municipio




hc3.dend <- as.dendrogram(var3.hc1)

jpeg(filename = "C:/Users/Antonia/Desktop/Documents/R/master/plots/hc3.var3.dend.jpeg", width = 100, height = 42, units = "cm", res = 300)
par(mar=c(12,8,4,2))
hc3.dend %>%
  set("labels_col", value = c("#BE0C15", "#00447F", "#FFC425",  "#A99A6F",  "#00703C", "#53247F" ), k = 6) %>%
  set("branches_k_color", value = c("#BE0C15", "#00447F", "#FFC425",  "#A99A6F",  "#00703C", "#53247F"), k = 6) %>%
  set("branches_lwd", 4) %>%
  plot(axes=TRUE, cex.main = 4, cex.axis = 2.5, lwd = 2, main = "Cluster Dendrogram", xlab = "", sub = "")
title(ylab = "Height", cex.lab = 3, line = 5)
rect.dendrogram(tree = hc3.dend, k=6, border ="grey", lty = 1, lwd = 2)
grid(nx = NA, ny = NULL, lty = 1, col = "gray", lwd = 1)
legend("topright",legend = paste("Cluster",1:6), fill= c("#C29DE3", "#85C5FF", "#99C5B1FF", "#D7D0BB","#FFE497", "#F9A1A5"), bty="n", border="white", cex = 3 )
dev.off()




#cut the cluster tree into groups
var3.hc1.groups <- as.factor(cutree(tree = var3.hc1, k = 6))


#bind with municipality names
var3_scores <- cbind(Municipio = municipio, var3_scores)

#write.csv(var3_scores, file = "C:/Users/Antonia/Desktop/Documents/R/master/var3_scores.csv")


#append cluster labels to original data
var3.hc1.data <- as.data.frame(cbind(var3_scores, cluster = as.numeric(var3.hc1.groups), deparse.level = 0))


#display first six rows of final data
head(var3.hc1.data)
write.csv(var3.hc1.data, file = "C:/Users/Antonia/Desktop/Documents/R/master/var3.hc1_data.csv")


#data preparation for GWR (see extra rsheet) ######################################################################
var3_scores.original <- cbind(Municipio = municipio, var3_scores.original)
score.clust <- as.data.frame(cbind(var3_scores.original, cluster = as.numeric(var3.hc1.groups), deparse.level = 0))
write.csv(score.clust, file = "C:/Users/Antonia/Desktop/Documents/R/master/score.clust.csv")
###################################################################################################################


#find mean values for each cluster
var3.hc1_mean <- aggregate(var3.hc1.data, by=list(cluster= var3.hc1.data$cluster), mean)
var3.hc1_mean <- var3.hc1_mean[,-2]
var3.hc1_mean <- var3.hc1_mean[,-16]


var3.hc1_mean


write.csv(var3.hc1_mean, file = "C:/Users/Antonia/Desktop/Documents/R/master/var3.hc1_mean_CID.csv")


#merge with shapefile
hc1.sf <- merge(shp_muni, var3.hc1.data, by="Municipio")




#map clusters, save as A4 portrait pdf
hc3.var3.C6 <- ggplot(data = hc1.sf) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = hc1.sf, aes(fill = as.factor(cluster))) +
  scale_fill_manual(values = MAcol_biglight(6)) +
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

ggsave(filename = "hc3.var3.C6.jpg", plot = hc3.var3.C6, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




var3.hc1_mean <- round(x = var3.hc1_mean, digits = 2)


lcols <- c("#53247F", "#00447F", "#00703C", "#A99A6F", "#FFC425", "#F99D31")


#plot radarplot gesamt 
hc3.rad.ges <- ggradar(var3.hc1_mean,
        plot.title = "PC Score Mean in Clusters 1 to 6",
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

ggsave(filename = "hc3.rad.ges.jpg", plot = hc3.rad.ges, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var3.hc1_meanC1 <- as.data.frame(var3.hc1_mean[1,1:15])

hc3.rad.C1 <- ggradar(var3.hc1_meanC1,
        plot.title = "PC Score Mean in Cluster 1",
        grid.min = -1,
        grid.mid = 0,
        grid.max = 1,
        values.radar = c(-1, 0, 1),
        axis.line.colour = "gray60",
        gridline.min.colour = "gray60",
        gridline.mid.colour = "gray60",
        gridline.max.colour = "gray60",
        group.colours = "#53247F",
        legend.title = "Cluster",
        legend.position = "bottom")

ggsave(filename = "hc3.rad.C1.jpg", plot = hc3.rad.C1, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var3.hc1_meanC2 <- as.data.frame(var3.hc1_mean[2,1:15])

hc3.rad.C2 <- ggradar(var3.hc1_meanC2,
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

ggsave(filename = "hc3.rad.C2.jpg", plot = hc3.rad.C2, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var3.hc1_meanC3 <- as.data.frame(var3.hc1_mean[3,1:15])

hc3.rad.C3 <- ggradar(var3.hc1_meanC3,
                      plot.title = "PC Score Mean in Cluster 3",
                      grid.min = -7,
                      grid.mid = 0,
                      grid.max = 7,
                      values.radar = c(-7, 0, 7),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#00703C")

ggsave(filename = "hc3.rad.C3.jpg", plot = hc3.rad.C3, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var3.hc1_meanC4 <- as.data.frame(var3.hc1_mean[4,1:15])

hc3.rad.C4 <- ggradar(var3.hc1_meanC4,
                      plot.title = "PC Score Mean in Cluster 4",
                      grid.min = -2,
                      grid.mid = 0,
                      grid.max = 2,
                      values.radar = c(-2, 0, 2),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#A99A6F")

ggsave(filename = "hc3.rad.C4.jpg", plot = hc3.rad.C4, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var3.hc1_meanC5 <- as.data.frame(var3.hc1_mean[5,1:15])

hc3.rad.C5 <- ggradar(var3.hc1_meanC5,
                      plot.title = "PC Score Mean in Cluster 5",
                      grid.min = -2,
                      grid.mid = 0,
                      grid.max = 2,
                      values.radar = c(-2, 0, 2),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#FFC425")

ggsave(filename = "hc3.rad.C5.jpg", plot = hc3.rad.C5, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)




var3.hc1_meanC6 <- as.data.frame(var3.hc1_mean[6,1:15])

hc3.rad.C6 <- ggradar(var3.hc1_meanC6,
                      plot.title = "PC Score Mean in Cluster 6",
                      grid.min = -16,
                      grid.mid = 0,
                      grid.max = 16,
                      values.radar = c(-16, 0, 16),
                      axis.line.colour = "gray60",
                      gridline.min.colour = "gray60",
                      gridline.mid.colour = "gray60",
                      gridline.max.colour = "gray60",
                      group.colours = "#F99D31")

ggsave(filename = "hc3.rad.C6.jpg", plot = hc3.rad.C6, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 21, units = "cm", dpi = 300)
