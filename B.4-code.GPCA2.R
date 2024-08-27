#Universidade de Coimbra - master thesis - Principal Component Analysis 2


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




#avoid
options(prompt="R> ", scipen=999)


#set working directory
setwd("C:/Users/Antonia/Desktop/Documents/R/master/")


#load data
variable2 <- read.csv("health outcomes data CID R.csv")
municipio <- variable2$name




#MISSING VALUES

#delete the first two columns ID and name
variable2 <- variable2[,-1:-2]

#substitute blanks with NA
variable2[variable2 == ""] <- NA

#see how many NA are in each variable / municipality in absolute numbers
colSums(is.na(variable2)) #variable
rowSums(is.na(variable2)) #municipality

#calculate the percentage of missing values per per column (2) and row (1)
pmiss <- function(x){sum(is.na(x))/length(x)*100}
apply(variable2,2,pmiss) #variable
apply(variable2,1,pmiss) #municipality

#delete variables with more than 5 % missing values from variable 2 (n = 9)
variable2 <- subset(variable2, select = -c(tumcol.pop,
tumpro.pop,
tumtec.pop,
pneu.pop,
cron.pop,
geni.pop,
rimur.pop,
peri.pop,
acid.pop))

  
#for all other variables with missing values a single imputation method (here mean substitution) is applied
  
#replace NA in all columns
for(i in 1:ncol(variable2)) {
variable2[ , i][is.na(variable2[ , i])] <- mean(variable2[ , i], na.rm = TRUE)
}




variable2 <- as.data.frame(scale(variable2))




#MULTIDIMENSIONAL SCALING

distance <- dist(variable2) 
mds_var2<-cmdscale(distance, k=2)

mds_var2 <- as.data.frame(mds_var2)
mds_var2 <- bind_cols(Municipio=municipio, mds_var2)


pca2.var2.mds <- 
  ggplot(mds_var2) +
  geom_point(aes(x = V1, y = V2)) + 
  geom_text_repel(data = mds_var2, aes(V1, V2, label = Municipio)) + #, position = position_nudge_repel(x = 0.1, y = 0.05)) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "grey", linewidth = 0.25, linetype = 1),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
  ) +
  ggtitle("Multidimensional Scaling") +
  theme(plot.title = element_text(size = 18)) +
  xlab("V1") +
  theme(axis.title.x = element_text(size = 14)) +
  ylab("V2") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(vjust = 1, size = 14, hjust = 1))+ 
  theme(axis.text.y = element_text(vjust = 1, size = 14, hjust = 1))

print(pca2.var2.mds)

ggsave(filename = "pca2.var2.mds.jpg", plot = pca2.var2.mds, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 29.7, height = 21, units = "cm", dpi = 300)




#CORRELATION MATRIX 1

var2_cormat1 <- round(cor(variable2),2)

#get lower triangle of the correlation matrix
get_lower_tri1 <- function(var2_cormat1){
  var2_cormat1[upper.tri(var2_cormat1)] <- NA
  return(var2_cormat1)
}

#get upper triangle of the correlation matrix
get_upper_tri1 <- function(var2_cormat1){
  var2_cormat1[lower.tri(var2_cormat1)]<- NA
  return(var2_cormat1)
}

# Use correlation between variables as distance
reorder_cormat1 <- function(var2_cormat1){
  dd <- as.dist((1-var2_cormat1)/2)
  hc <- hclust(dd)
  var2_cormat1 <-var2_cormat1[hc$order, hc$order]
}

#reorder the correlation matrix
var2_cormat1 <- reorder_cormat1(var2_cormat1)
upper_tri1 <- get_upper_tri1(var2_cormat1)

#melt the correlation matrix
melt_cormat1 <- melt(upper_tri1, na.rm = TRUE)

# Create a ggheatmap
pca2.var2.cormat1 <- 
  ggplot(melt_cormat1, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson Correlation Coefficient") +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 0.8) +
  ggtitle("Correlation Matrix") +
  theme(text = element_text(size = 16)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.25, linetype = 1),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0.5, 0),
    legend.position = "bottom",
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 0.5, size = 6, hjust = 1))+
  coord_fixed()


print(pca2.var2.cormat1)

ggsave(filename = "pca2.var2.cormat1.jpg", plot = pca2.var2.cormat1, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 26, units = "cm", dpi = 300)

write.csv(var2_cormat1, "C:/Users/Antonia/Desktop/Documents/R/master/var2_cormat1.csv", row.names=TRUE)





#GLOBAL PRINCIPAL COMPONENT ANALYSIS 2 - PCA 1

#var2 n = 62

var2_pca1 <- pca(r = variable2, nfactors = 62, residuals = FALSE, rotate = "none", n.obs = 278, covar = FALSE, scores = TRUE, oblique.scores = FALSE, cor = "cor")
summary(var2_pca1)

var2_pca1_load1 <- var2_pca1$loadings
print(var2_pca1_load1)

var2_pca1_eig1 <- var2_pca1$values
print(var2_pca1_eig1)


write.csv(var2_pca1_load1, file = "C:/Users/Antonia/Desktop/Documents/R/master/var2_pca1_load1_CID.csv")


#variables with eigenvalue >=1 are n=24
#delete variables that in all first 24 components provide loadings between x[-0.5;0.5] (n = 39)

var2_1 <- subset(variable2, select = -c(
  par.pop,
  dip.pop,
  hiv.pop,
  tum.pop,
  tumlab.pop,
  tumeso.pop,
  tumrec.pop,
  tumfig.pop,
  tumlar.pop,
  melpel.pop,
  tummam.pop,
  tumcolut.pop,
  tumute.pop,
  tumova.pop,
  tumrim.pop,
  tumbex.pop,
  sang.pop,
  endo.pop,
  mental.pop,
  drog.pop,
  grip.pop,
  asma.pop,
  ilc.pop,
  figa.pop,
  pele.pop,
  oste.pop,
  art.pop,
  grav.pop,
  mfcrom.pop,
  mfnerv.pop,
  mfcirc.pop,
  sint.pop,
  mlact.pop,
  mdes.pop,
  acidtr.pop,
  queda.pop,
  enven.pop,
  hom.pop,
  leso.pop
))




#GLOBAL PRINCIPAL COMPONENT ANALYSIS 2 - PCA 2

#var2 n = 23

var2_pca2 <- pca(r = var2_1, nfactors = 23, residuals = FALSE, rotate = "none", n.obs = 278, covar = FALSE, scores = TRUE, oblique.scores = FALSE, cor = "cor")
summary(var2_pca2)

var2_pca2_load2 <- var2_pca2$loadings
print(var2_pca2_load2)

var2_pca2_eig2 <- var2_pca2$values
print(var2_pca2_eig2)


write.csv(var2_pca2_load2, file = "C:/Users/Antonia/Desktop/Documents/R/master/var2_pca2_load2_CID.csv")


#variables with eigenvalue >=1 are n=6
#delete variables that in all first 6 components provide loadings between x[-0.5;0.5] (n = 1)

var2_2 <- subset(var2_1, select = -c(
  tumest.pop
))




#GLOBAL PRINCIPAL COMPONENT ANALYSIS 2 - PCA 3

#var2 n = 22

var2_pca3 <- pca(r = var2_2, nfactors = 22, residuals = FALSE, rotate = "none", n.obs = 278, covar = FALSE, scores = TRUE, oblique.scores = FALSE, cor = "cor")
summary(var2_pca3)

var2_pca3_load3 <- var2_pca3$loadings
print(var2_pca3_load3)

var2_pca3_eig3 <- var2_pca3$values
print(var2_pca3_eig3)


write.csv(var2_pca3_load3, file = "C:/Users/Antonia/Desktop/Documents/R/master/var2_pca3_load3_CID.csv")


#variables with eigenvalue >=1 are n=6
#delete variables that in all first 6 components provide loadings between x[-0.5;0.5] (n = 1)

var2_3 <- subset(var2_2, select = -c(
  hv.pop
))




#GLOBAL PRINCIPAL COMPONENT ANALYSIS 2 - PCA 4

#var2 n = 21

var2_pca4 <- pca(r = var2_3, nfactors = 21, residuals = FALSE, rotate = "none", n.obs = 278, covar = FALSE, scores = TRUE, oblique.scores = FALSE, cor = "cor")
print.psych(var2_pca4)

#get the loadings
var2_pca4_load4 <- var2_pca4$loadings
print(var2_pca4_load4)

var2_pca4_load4 <- as.data.frame(var2_pca4_load4[1:21,1:21])

#get the quality of representation (cos2)
cos2 <- as.data.frame(var2_pca4_load4^2)

#get the communalities
communalities <- as.data.frame(rowSums(cos2[1:21,1:6]))

#see proportion of total variance explained per PC
prop.table(var2_pca4$values)

#get eigenvalues
var2_pca4_eig4 <- var2_pca4$values

#interchange columns and rows and create seperate dataframe
eig4.df <- data.frame(t(var2_pca4_eig4))

#copy column values 21 times
eig4.copy <- data.frame(eig4.df[,1:21], ntimes=c(21))
eig4.copy <- as.data.frame(lapply(eig4.copy, rep, eig4.copy$ntimes))
eig4.copy <- eig4.copy[,-22]

#get the contribution of the variables
contribution <- as.data.frame(((var2_pca4_load4^2)/eig4.copy)*100)

var2_pca4_eig4 <- as.data.frame(var2_pca4_eig4[1:10])
print(var2_pca4_eig4)


#data preparation screeplot
var2_pca4_eig4 <- data.frame(comp10 = paste(c("PC"),c(1:10), sep = ""), value = round(var2_pca4_eig4[,1],1))


var2_scores <- var2_pca4$scores[,1:6]
var2_scores <- as.data.frame(var2_scores)


write.csv(var2_pca4_load4, file = "C:/Users/Antonia/Desktop/Documents/R/master/var2_pca4_load4_CID.csv")
write.csv(cos2, file = "C:/Users/Antonia/Desktop/Documents/R/master/var2_pca4_cos2.csv")
write.csv(communalities, file = "C:/Users/Antonia/Desktop/Documents/R/master/var2_pca4_communalities.csv")
write.csv(contribution, file = "C:/Users/Antonia/Desktop/Documents/R/master/var2_pca4_contribution.csv")
write.csv(var2_scores, file = "C:/Users/Antonia/Desktop/Documents/R/master/var2_scores_CID.csv")




#BIPLOT PCA 4

pca2.var2.biplot <- ggplot(var2_pca4_load4, aes(x = PC1, y = PC2, group = 1)) +
  geom_point(colour = "#000000") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), linewidth = 0.75, arrow = arrow(length = unit(0.3, "cm"))) +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Biplot") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  geom_text_repel(aes(label = rownames(var2_pca4_load4), vjust= 1, hjust = 0, size = 3)) +
  xlab("PC1 (33.6 %)") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("PC2 (10.1 %)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(expand = FALSE, xlim = c(-0.1,1.25), ylim = c(-0.8,0.8))

print(pca2.var2.biplot)

ggsave(filename = "pca2.var2.biplot.jpg", plot = pca2.var2.biplot, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 20, height = 20, units = "cm", dpi = 300)




#SCREEPLOT (EIGENVALUES) PCA 4
pca2.var2.screeplot <- ggplot(var2_pca4_eig4, aes(x = reorder(comp10,-value), y = value, fill = value, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  geom_line(colour = "#000000") +
  geom_point(colour = "#000000") +
  geom_hline(yintercept = 1, linewidth = 1, colour =  "#BE0C15") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Scree Plot") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  geom_text(aes(label = value), vjust= -0.5, hjust = 0, size = 3) +
  xlab("Princiapl Components") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Eigenvalues (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.5, 10))

print(pca2.var2.screeplot)

ggsave(filename = "pca2.var2.screeplot.jpg", plot = pca2.var2.screeplot, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




variance <- data.frame(
  comp5 = c("PC1","PC2","PC3","PC4","PC5","PC6"),
  PTV = round(c(0.336, 0.101, 0.064, 0.053, 0.052, 0.048)*100,1),
  CPTV = round(c(0.336, 0.436, 0.501, 0.553, 0.605, 0.653)*100,1)
)




#BARPLOT PVT PCA 4

pca2.var2.PTV <- ggplot(variance, aes(x = comp5, y = PTV, fill = PTV)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Proportion of Total Variance") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  geom_text(aes(label = PTV), vjust= -0.5, size = 3) +
  xlab("Princiapl Components") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("PTV (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.5, 42))

print(pca2.var2.PTV)

ggsave(filename = "pca2.var2.PTV.jpg", plot = pca2.var2.PTV, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#BARPLOT CPVT PCA 4
pca2.var2.CPTV <- ggplot(variance, aes(x = "", y = PTV, fill = comp5, group = PTV)) +
  geom_col() +
  scale_fill_manual(values = c("#C29DE3", "#85C5FF", "#95BFC5", "#C1FFE1", "#D7D0BB","#FFE497", "#FBC27D", "#F9A1A5") ) +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Cumulative Proportion of Total Variance") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(title="Principal \n Components")) +
  theme(panel.border = element_blank()) +
  geom_text(aes(label = PTV), size = 5, colour = "#000000", position = position_stack(vjust = 0.5)) +
  xlab(label = "Stacked Princiapl Components") +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("CPTV (%)") +
  theme(axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)))

print(pca2.var2.CPTV)

ggsave(filename = "pca2.var2.CPTV.jpg", plot = pca2.var2.CPTV, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#BARPLOT COMMUNALITIES PC1-PC6

pca2.var2.communalities <- ggplot(communalities, aes(x = reorder(rownames(communalities), -rowSums(cos2[1:21,1:6])), y = rowSums(cos2[1:21,1:6]), fill = rowSums(cos2[1:21,1:6]), group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Communalities of Variables in PC1 to PC6") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  #  geom_text(aes(label = var1names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Communality") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 1))

print(pca2.var2.communalities)

ggsave(filename = "pca2.var2.communalities.jpg", plot = pca2.var2.communalities, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




PC1.contribution <- as.data.frame(cbind(var2names = row.names(contribution), PC1cont = contribution$PC1))
PC1.contribution$PC1cont <- as.numeric(as.character(PC1.contribution$PC1cont))
str(PC1.contribution)


PC2.contribution <- as.data.frame(cbind(var2names = row.names(contribution), PC2cont = contribution$PC2))
PC2.contribution$PC2cont <- as.numeric(as.character(PC2.contribution$PC2cont))
str(PC2.contribution)


PC3.contribution <- as.data.frame(cbind(var2names = row.names(contribution), PC3cont = contribution$PC3))
PC3.contribution$PC3cont <- as.numeric(as.character(PC3.contribution$PC3cont))
str(PC3.contribution)


PC4.contribution <- as.data.frame(cbind(var2names = row.names(contribution), PC4cont = contribution$PC4))
PC4.contribution$PC4cont <- as.numeric(as.character(PC4.contribution$PC4cont))
str(PC4.contribution)


PC5.contribution <- as.data.frame(cbind(var2names = row.names(contribution), PC5cont = contribution$PC5))
PC5.contribution$PC5cont <- as.numeric(as.character(PC5.contribution$PC5cont))
str(PC5.contribution)


PC6.contribution <- as.data.frame(cbind(var2names = row.names(contribution), PC6cont = contribution$PC6))
PC6.contribution$PC6cont <- as.numeric(as.character(PC6.contribution$PC6cont))
str(PC6.contribution)




#PLOT VARIABLE CONTRIBUTION PC1

pca2.var2.PC1cont <- ggplot(PC1.contribution, aes(x = reorder(var2names, -PC1cont), y = PC1cont, fill = PC1cont, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Proportion of Variable Contribution to PC1") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  #  geom_text(aes(label = var2names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Contribution (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 15))

print(pca2.var2.PC1cont)

ggsave(filename = "pca2.var2.PC1cont.jpg", plot = pca2.var2.PC1cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC2

pca2.var2.PC2cont <- ggplot(PC2.contribution, aes(x = reorder(var2names, -PC2cont), y = PC2cont, fill = PC2cont, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Proportion of Variable Contribution to PC2") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  #  geom_text(aes(label = var1names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Contribution (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 15))

print(pca2.var2.PC2cont)

ggsave(filename = "pca2.var2.PC2cont.jpg", plot = pca2.var2.PC2cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC3

pca2.var2.PC3cont <- ggplot(PC3.contribution, aes(x = reorder(var2names, -PC3cont), y = PC3cont, fill = PC3cont, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Proportion of Variable Contribution to PC3") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  #  geom_text(aes(label = var1names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Contribution (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 30))

print(pca2.var2.PC3cont)

ggsave(filename = "pca2.var2.PC3cont.jpg", plot = pca2.var2.PC3cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC4

pca2.var2.PC4cont <- ggplot(PC4.contribution, aes(x = reorder(var2names, -PC4cont), y = PC4cont, fill = PC4cont, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Proportion of Variable Contribution to PC4") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  #  geom_text(aes(label = var1names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Contribution (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 65))

print(pca2.var2.PC4cont)

ggsave(filename = "pca2.var2.PC4cont.jpg", plot = pca2.var2.PC4cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC5

pca2.var2.PC5cont <- ggplot(PC5.contribution, aes(x = reorder(var2names, -PC5cont), y = PC5cont, fill = PC5cont, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Proportion of Variable Contribution to PC5") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  #  geom_text(aes(label = var1names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Contribution (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 40))

print(pca2.var2.PC5cont)

ggsave(filename = "pca2.var2.PC5cont.jpg", plot = pca2.var2.PC5cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC6

pca2.var2.PC6cont <- ggplot(PC6.contribution, aes(x = reorder(var2names, -PC6cont), y = PC6cont, fill = PC6cont, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Proportion of Variable Contribution to PC6") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  #  geom_text(aes(label = var1names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Contribution (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 90))

print(pca2.var2.PC6cont)

ggsave(filename = "pca2.var2.PC6cont.jpg", plot = pca2.var2.PC6cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#bind with municipality names
var2_scores <- bind_cols(Municipio=municipio, var2_scores)

write.csv(var2_scores, file = "C:/Users/Antonia/Desktop/Documents/R/master/var2_scores_muninames.csv")


#read shapefiles
shp_muni <- read_sf("Cont_Mun_CAOP2023.shp")

#delete non useful columns
shp_muni <- shp_muni[,c(-1:-2,-4:-12)]

#data2 = standardized variables of health outcomes merged with shape file

#merge scaled variable data frame with shapefile
pca4.sfdata <- merge(shp_muni,var2_scores, by="Municipio")

#delete first row, which is of type "character"
#data2 <- data2[,-1]
#str(data2)




#PLOT SCORE MAP PC1

pca2.var2.PC1score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca4.sfdata, aes(fill = PC1)) +
  scale_fill_gradient2(low = "#F99D31", high = "#53247F", mid = "#ffffff", 
                       midpoint = 0) + #limit = c(-1,10), space = "Lab"
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("PC1") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Score Value") +
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

ggsave(filename = "pca2.var2.PC1score.jpg", plot = pca2.var2.PC1score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC2

pca2.var2.PC2score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca4.sfdata, aes(fill = PC2)) +
  scale_fill_gradient2(low = "#F99D31", high = "#53247F", mid = "#ffffff", 
                       midpoint = 0) + #limit = c(-1,10), space = "Lab"
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("PC2") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Score Value") +
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

ggsave(filename = "pca2.var2.PC2score.jpg", plot = pca2.var2.PC2score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC3

pca2.var2.PC3score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca4.sfdata, aes(fill = PC3)) +
  scale_fill_gradient2(low = "#F99D31", high = "#53247F", mid = "#ffffff", 
                       midpoint = 0) + #limit = c(-1,10), space = "Lab"
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("PC3") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Score Value") +
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

ggsave(filename = "pca2.var2.PC3score.jpg", plot = pca2.var2.PC3score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC4

pca2.var2.PC4score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca4.sfdata, aes(fill = PC4)) +
  scale_fill_gradient2(low = "#F99D31", high = "#53247F", mid = "#ffffff", 
                       midpoint = 0) + #limit = c(-1,10), space = "Lab"
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("PC4") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Score Value") +
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

ggsave(filename = "pca2.var2.PC4score.jpg", plot = pca2.var2.PC4score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC5

pca2.var2.PC5score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca4.sfdata, aes(fill = PC5)) +
  scale_fill_gradient2(low = "#F99D31", high = "#53247F", mid = "#ffffff", 
                       midpoint = 0) + #limit = c(-1,10), space = "Lab"
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("PC5") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Score Value") +
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

ggsave(filename = "pca2.var2.PC5score.jpg", plot = pca2.var2.PC5score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC6

pca2.var2.PC6score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca4.sfdata, aes(fill = PC6)) +
  scale_fill_gradient2(low = "#F99D31", high = "#53247F", mid = "#ffffff", 
                       midpoint = 0) + #limit = c(-1,10), space = "Lab"
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("PC6") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(fill="Score Value") +
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

ggsave(filename = "pca2.var2.PC6score.jpg", plot = pca2.var2.PC6score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT LOADING AND SCORE POSITIONS

plot.scores1 <- as.data.frame(cbind(municipio = municipio, scores = var2_scores$PC1))
plot.scores1$scores <- as.numeric(as.character(plot.scores1$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-5, 5)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca2.var2.PC1pos <- ggplot() +
  geom_point(data = var2_pca4_load4, aes(x = 1, y = PC1)) + 
  geom_text_repel(data = var2_pca4_load4, aes(1, PC1, label = row.names(var2_pca4_load4)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
  geom_point(data = plot.scores1, aes(x = 2, y = a+scores*b)) +
  geom_text_repel(data = plot.scores1, aes(2, a+scores*b, label = municipio), position = position_nudge_repel(x = 0.1, y = 0), size = 3) +
  scale_y_continuous("Loadings", sec.axis = sec_axis(~ (. - a)/b, name = "Scores")) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "grey", size = 0.25, linetype = 1),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
  ) +
  ggtitle("Position of variables and observations on PC1") +
  theme(plot.title = element_text(size = 18, hjust=0)) +
  xlab("PC1") +
  theme(axis.title.x = element_text(size = 14)) +
  ylab("Loadings") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(vjust = 1, size = 14, hjust = 1)) + 
  theme(axis.text.y = element_text(vjust = 1, size = 14, hjust = 1)) +
  coord_cartesian(expand = FALSE, xlim = c(0,3), ylim = c(-1.2,1.2))

print(pca2.var2.PC1pos)

ggsave(filename = "pca2.var2.PC1pos.jpg", plot = pca2.var2.PC1pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores2 <- as.data.frame(cbind(municipio = municipio, scores = var2_scores$PC2))
plot.scores2$scores <- as.numeric(as.character(plot.scores2$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-6, 6)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca2.var2.PC2pos <- ggplot() +
  geom_point(data = var2_pca4_load4, aes(x = 1, y = PC2)) + 
  geom_text_repel(data = var2_pca4_load4, aes(1, PC2, label = row.names(var2_pca4_load4)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
  geom_point(data = plot.scores2, aes(x = 2, y = a+scores*b)) +
  geom_text_repel(data = plot.scores2, aes(2, a+scores*b, label = municipio), position = position_nudge_repel(x = 0.1, y = 0), size = 3) +
  scale_y_continuous("Loadings", sec.axis = sec_axis(~ (. - a)/b, name = "Scores")) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "grey", size = 0.25, linetype = 1),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
  ) +
  ggtitle("Position of variables and observations on PC2") +
  theme(plot.title = element_text(size = 18, hjust=0)) +
  xlab("PC1") +
  theme(axis.title.x = element_text(size = 14)) +
  ylab("Loadings") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(vjust = 1, size = 14, hjust = 1)) + 
  theme(axis.text.y = element_text(vjust = 1, size = 14, hjust = 1)) +
  coord_cartesian(expand = FALSE, xlim = c(0,3), ylim = c(-1.2,1.2))

print(pca2.var2.PC2pos)

ggsave(filename = "pca2.var2.PC2pos.jpg", plot = pca2.var2.PC2pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores3 <- as.data.frame(cbind(municipio = municipio, scores = var2_scores$PC3))
plot.scores3$scores <- as.numeric(as.character(plot.scores3$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-6, 6)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca2.var2.PC3pos <- ggplot() +
  geom_point(data = var2_pca4_load4, aes(x = 1, y = PC3)) + 
  geom_text_repel(data = var2_pca4_load4, aes(1, PC3, label = row.names(var2_pca4_load4)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
  geom_point(data = plot.scores3, aes(x = 2, y = a+scores*b)) +
  geom_text_repel(data = plot.scores3, aes(2, a+scores*b, label = municipio), position = position_nudge_repel(x = 0.1, y = 0), size = 3) +
  scale_y_continuous("Loadings", sec.axis = sec_axis(~ (. - a)/b, name = "Scores")) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "grey", size = 0.25, linetype = 1),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
  ) +
  ggtitle("Position of variables and observations on PC3") +
  theme(plot.title = element_text(size = 18, hjust=0)) +
  xlab("PC1") +
  theme(axis.title.x = element_text(size = 14)) +
  ylab("Loadings") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(vjust = 1, size = 14, hjust = 1)) + 
  theme(axis.text.y = element_text(vjust = 1, size = 14, hjust = 1)) +
  coord_cartesian(expand = FALSE, xlim = c(0,3), ylim = c(-1.2,1.2))

print(pca2.var2.PC3pos)

ggsave(filename = "pca2.var2.PC3pos.jpg", plot = pca2.var2.PC3pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores4 <- as.data.frame(cbind(municipio = municipio, scores = var2_scores$PC4))
plot.scores4$scores <- as.numeric(as.character(plot.scores4$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-10, 10)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca2.var2.PC4pos <- ggplot() +
  geom_point(data = var2_pca4_load4, aes(x = 1, y = PC4)) + 
  geom_text_repel(data = var2_pca4_load4, aes(1, PC4, label = row.names(var2_pca4_load4)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
  geom_point(data = plot.scores4, aes(x = 2, y = a+scores*b)) +
  geom_text_repel(data = plot.scores4, aes(2, a+scores*b, label = municipio), position = position_nudge_repel(x = 0.1, y = 0), size = 3) +
  scale_y_continuous("Loadings", sec.axis = sec_axis(~ (. - a)/b, name = "Scores")) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "grey", size = 0.25, linetype = 1),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
  ) +
  ggtitle("Position of variables and observations on PC4") +
  theme(plot.title = element_text(size = 18, hjust=0)) +
  xlab("PC1") +
  theme(axis.title.x = element_text(size = 14)) +
  ylab("Loadings") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(vjust = 1, size = 14, hjust = 1)) + 
  theme(axis.text.y = element_text(vjust = 1, size = 14, hjust = 1)) +
  coord_cartesian(expand = FALSE, xlim = c(0,3), ylim = c(-1.2,1.2))

print(pca2.var2.PC4pos)

ggsave(filename = "pca2.var2.PC4pos.jpg", plot = pca2.var2.PC4pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores5 <- as.data.frame(cbind(municipio = municipio, scores = var2_scores$PC5))
plot.scores5$scores <- as.numeric(as.character(plot.scores5$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-6, 6)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca2.var2.PC5pos <- ggplot() +
  geom_point(data = var2_pca4_load4, aes(x = 1, y = PC5)) + 
  geom_text_repel(data = var2_pca4_load4, aes(1, PC5, label = row.names(var2_pca4_load4)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
  geom_point(data = plot.scores5, aes(x = 2, y = a+scores*b)) +
  geom_text_repel(data = plot.scores5, aes(2, a+scores*b, label = municipio), position = position_nudge_repel(x = 0.1, y = 0), size = 3) +
  scale_y_continuous("Loadings", sec.axis = sec_axis(~ (. - a)/b, name = "Scores")) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "grey", size = 0.25, linetype = 1),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
  ) +
  ggtitle("Position of variables and observations on PC5") +
  theme(plot.title = element_text(size = 18, hjust=0)) +
  xlab("PC1") +
  theme(axis.title.x = element_text(size = 14)) +
  ylab("Loadings") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(vjust = 1, size = 14, hjust = 1)) + 
  theme(axis.text.y = element_text(vjust = 1, size = 14, hjust = 1)) +
  coord_cartesian(expand = FALSE, xlim = c(0,3), ylim = c(-1.2,1.2))

print(pca2.var2.PC5pos)

ggsave(filename = "pca2.var2.PC5pos.jpg", plot = pca2.var2.PC5pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores6 <- as.data.frame(cbind(municipio = municipio, scores = var2_scores$PC6))
plot.scores6$scores <- as.numeric(as.character(plot.scores6$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-18, 18)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca2.var2.PC6pos <- ggplot() +
  geom_point(data = var2_pca4_load4, aes(x = 1, y = PC6)) + 
  geom_text_repel(data = var2_pca4_load4, aes(1, PC6, label = row.names(var2_pca4_load4)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
  geom_point(data = plot.scores6, aes(x = 2, y = a+scores*b)) +
  geom_text_repel(data = plot.scores6, aes(2, a+scores*b, label = municipio), position = position_nudge_repel(x = 0.1, y = 0), size = 3) +
  scale_y_continuous("Loadings", sec.axis = sec_axis(~ (. - a)/b, name = "Scores")) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "grey", size = 0.25, linetype = 1),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
  ) +
  ggtitle("Position of variables and observations on PC6") +
  theme(plot.title = element_text(size = 18, hjust=0)) +
  xlab("PC1") +
  theme(axis.title.x = element_text(size = 14)) +
  ylab("Loadings") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(vjust = 1, size = 14, hjust = 1)) + 
  theme(axis.text.y = element_text(vjust = 1, size = 14, hjust = 1)) +
  coord_cartesian(expand = FALSE, xlim = c(0,3), ylim = c(-1.2,1.2))

print(pca2.var2.PC6pos)

ggsave(filename = "pca2.var2.PC6pos.jpg", plot = pca2.var2.PC6pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#CORRELATION MATRIX 2

var2_cormat2 <- round(cor(var2_3),2)

#get lower triangle of the correlation matrix
get_lower_tri2 <- function(var2_cormat2){
  var2_cormat2[upper.tri(var2_cormat2)] <- NA
  return(var2_cormat2)
}

#get upper triangle of the correlation matrix
get_upper_tri2 <- function(var2_cormat2){
  var2_cormat2[lower.tri(var2_cormat2)]<- NA
  return(var2_cormat2)
}

# Use correlation between variables as distance
reorder_cormat2 <- function(var2_cormat2){
  dd <- as.dist((1-var2_cormat2)/2)
  hc <- hclust(dd)
  var2_cormat2 <-var2_cormat2[hc$order, hc$order]
}

#reorder the correlation matrix
var2_cormat2 <- reorder_cormat2(var2_cormat2)
upper_tri2 <- get_upper_tri2(var2_cormat2)

#melt the correlation matrix
melt_cormat2 <- melt(upper_tri2, na.rm = TRUE)

# Create a ggheatmap
pca2.var2.cormat2 <- 
  ggplot(melt_cormat2, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson Correlation Coefficient") +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2.5) +
  ggtitle("Correlation Matrix") +
  theme(text = element_text(size = 16)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.25, linetype = 1),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0.5, 0),
    legend.position = "bottom",
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 0.5, size = 12, hjust = 1))+
  coord_fixed()


print(pca2.var2.cormat2)

ggsave(filename = "pca2.var2.cormat2.jpg", plot = pca2.var2.cormat2, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 26, units = "cm", dpi = 300)

write.csv(var2_cormat2, "C:/Users/Antonia/Desktop/Documents/R/master/var2_cormat2.csv", row.names=TRUE)
