#Universidade de Coimbra - master thesis - Global Principal Component Analysis 1


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




#avoid scientific notation
options(prompt="R> ", scipen=999)


#set working directory
setwd("C:/Users/Antonia/Desktop/Documents/R/master/")


#load data
variable1 <- read.csv("health determinants data.csv")
municipio <- variable1$name




#MISSING VALUES

#delete the first two columns ID and name
variable1 <- variable1[,-1:-2]

#substitute blanks with NA
variable1[variable1 == ""] <- NA

#see how many NA are in each variable / municipality in absolute numbers
colSums(is.na(variable1)) #variable
rowSums(is.na(variable1)) #municipality

#calculate the percentage of missing values per per column (2) and row (1)
pmiss <- function(x){sum(is.na(x))/length(x)*100}
apply(variable1,2,pmiss) #variable
apply(variable1,1,pmiss) #municipality

#delete variables with more than 5 % missing values from variable 1 (n = 9))
variable1 <- subset(variable1, select = -c(ncr.house,
                                           r_edu_sec,
                                           robb.cri,
                                           thef.cri,
                                           alco.cri,
                                           stat.cri,
                                           pets.cri,
                                           legis.cri,
                                           lice.cri))


#for all other variables with missing values a single imputation method (here mean substitution) is applied

#replace NA in all columns
for(i in 1:ncol(variable1)) {
  variable1[ , i][is.na(variable1[ , i])] <- mean(variable1[ , i], na.rm = TRUE)
}




variable1 <- as.data.frame(scale(variable1))




#MULTIDIMENSIONAL SCALING

distance <- dist(variable1)
mds_var1<-cmdscale(distance, k=2)

mds_var1 <- as.data.frame(mds_var1)
mds_var1 <- bind_cols(Municipio=municipio, mds_var1)


pca1.var1.mds <- 
  ggplot(mds_var1) +
  geom_point(aes(x = V1, y = V2)) + 
  geom_text_repel(data = mds_var1, aes(V1, V2, label = Municipio), position = position_nudge_repel(x = 0.1, y = 0.05)) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "grey", size = 0.25, linetype = 1),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
  ) +
  ggtitle("Multidimensional Scaling") +
  theme(plot.title = element_text(size = 18, hjust=0)) +
  xlab("V1") +
  theme(axis.title.x = element_text(size = 14)) +
  ylab("V2") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(vjust = 1, size = 14, hjust = 1))+ 
  theme(axis.text.y = element_text(vjust = 1, size = 14, hjust = 1))

print(pca1.var1.mds)

ggsave(filename = "pca1.var1.mds.jpg", plot = pca1.var1.mds, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 29.7, height = 21, units = "cm", dpi = 300)




#CORRELATION MATRIX 1

var1_cormat <- round(cor(variable1),2)

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
pca1.var1.cormat1 <- 
  ggplot(melt_cormat1, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson Correlation Coefficient") +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 0.5) +
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 0.5, size = 5, hjust = 1))+
  coord_fixed()


print(pca1.var1.cormat1)

ggsave(filename = "pca1.var1.cormat1.jpg", plot = pca1.var1.cormat1, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 26, units = "cm", dpi = 300)

write.csv(var1_cormat, "C:/Users/Antonia/Desktop/Documents/R/master/var1_cormat.csv", row.names=TRUE)




#GLOBAL PRINCIPAL COMPONENT ANALYSIS 1 - PCA 1

#var1 n = 144

var1_pca1 <- principal(r = variable1, nfactors = 144, residuals = FALSE, rotate = "none", n.obs = 278, covar = FALSE, scores = TRUE, oblique.scores = FALSE, cor = "cor")
summary(var1_pca1)
print(var1_pca1)

var1_pca1_load1 <- var1_pca1$loadings
print(var1_pca1_load1)

var1_pca1_eig1 <- var1_pca1$values
print(var1_pca1_eig1)

write.csv(var1_pca1_load1, file = "C:/Users/Antonia/Desktop/Documents/R/master/var1_pca1_load1.csv")


#variables with eigenvalue >=1 are n=23
#delete variables that in all first 23 components provide loadings between x[-0.6;0.6] (n = 63)

var1_1 <- subset(variable1, select = -c(
  nores.bui,
  to1945.bui,
  X1946.1980.bui,
  X1981.2000.bui,
  X2001.2010.bui,
  X2011.2021.bui,
  p_bui_rep,
  ch.acc,
  fh.acc,
  x.10mc.acc,
  X40.60mc.acc,
  x.30m.acc,
  X30.39m.acc,
  X80.99m.acc,
  X100.119m.acc,
  X150.199m.acc,
  res2.acc,
  p_acc_empty,
  p_house_unem,
  dad.fam,
  p_pop_desem,
  A.empl,
  B.empl,
  D.empl,
  E.empl,
  F.empl,
  G.empl,
  H.empl,
  P.empl,
  Q.empl,
  S.empl,
  T.empl,
  U.empl,
  a_com_tim,
  marf.pop,
  marm.pop,
  mardiv.pop,
  p_exp_cul,
  i_gini_inco,
  r_edu_bas,
  n_fire_rural,
  cri.pop,
  pers.cri,
  inju.cri,
  vio.cri,
  prop.cri,
  soci.cri,
  dis_ear_sex,
  agri.ter,
  past.ter,
  agro.ter,
  fore.ter,
  scru.ter,
  open.ter,
  wetl.ter,
  wate.ter,
  a_temp,
  ort.rel,
  wit.rel,
  bud.rel,
  hin.rel,
  jew.rel,
  onc.rel
))




#GLOBAL PRINCIPAL COMPONENT ANALYSIS 1 - PCA 2

#var1 = 81 because 63 were deleted

var1_pca2 <- pca(r = var1_1, nfactors = 81, residuals = FALSE, rotate = "none", n.obs = 278, covar = FALSE, scores = TRUE, oblique.scores = FALSE, cor = "cor")
summary(var1_pca2)

var1_pca2_load2 <- var1_pca2$loadings
print(var1_pca2_load2)

var1_pca2_eig2 <- var1_pca2$values
print(var1_pca2_eig2)

write.csv(var1_pca2_load2, file = "C:/Users/Antonia/Desktop/Documents/R/master/var1_pca2_load2.csv")


#variables with eigenvalue >=1 are n=12
#delete variables that in all first 12 components provide loadings between x[-0.6;0.6] (n = 8)

var1_2 <- subset(var1_1, select = -c(
  ac.acc,
  noac.acc,
  ho.acc,
  p_pop_for,
  I.empl,
  phys.cri,
  dis_ear_sec,
  mus.rel
))




#GLOBAL PRINCIPAL COMPONENT ANALYSIS 1 - PCA 3

#var1 = 73 because 8 were deleted

var1_pca3 <- principal(r = var1_2, nfactors = 73, residuals = FALSE, rotate = "none", n.obs = 278, covar = FALSE, scores = TRUE, oblique.scores = FALSE, cor = "cor")
print.psych(var1_pca3)

#get the loadings
var1_pca3_load3 <- var1_pca3$loadings
print(var1_pca3_load3)

var1_pca3_load3 <- as.data.frame(var1_pca3_load3[1:73,1:73])

#get the quality of representation (cos2)
cos2 <- as.data.frame(var1_pca3_load3^2)

#get the communalities
communalities <- as.data.frame(rowSums(cos2[1:73,1:8]))

#see proportion of total variance explained per PC
prop.table(var1_pca3$values)

#get the eigenvalues
var1_pca3_eig3 <- as.data.frame(var1_pca3$values)

#interchange columns and rows 
eig3.df <- data.frame(t(var1_pca3_eig3))

#copy column values 73 times
eig3.copy <- data.frame(eig3.df[,1:73], ntimes=c(73))
eig3.copy <- as.data.frame(lapply(eig3.copy, rep, eig3.copy$ntimes))
eig3.copy <- eig3.copy[,-74]

#get the contribution of the variables
contribution <- as.data.frame(((var1_pca3_load3^2)/eig3.copy)*100)

var1_pca3_eig3 <- as.data.frame(var1_pca3_eig3[1:10,1])


#data preparation screeplot
var1_pca3_eig3 <- data.frame(comp10 = paste(c("PC"),c(1:10), sep = ""), value = round(var1_pca3_eig3[,1],1))
print(var1_pca3_eig3)


var1_scores <- var1_pca3$scores[,1:8]
var1_scores <- as.data.frame(var1_scores)


write.csv(var1_pca3_load3, file = "C:/Users/Antonia/Desktop/Documents/R/master/var1_pca3_load3.csv")
write.csv(cos2, file = "C:/Users/Antonia/Desktop/Documents/R/master/var1_pca3_cos2.csv")
write.csv(communalities, file = "C:/Users/Antonia/Desktop/Documents/R/master/var1_pca3_communalities.csv")
write.csv(contribution, file = "C:/Users/Antonia/Desktop/Documents/R/master/var1_pca3_contribution.csv")
write.csv(var1_scores, file = "C:/Users/Antonia/Desktop/Documents/R/master/var1_scores.csv")




#BIPLOT PCA 3

pca1.var1.biplot <- ggplot(var1_pca3_load3, aes(x = PC1, y = PC2, group = 1)) +
  geom_point(colour = "#000000") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), linewidth = 0.75, arrow = arrow(length = unit(0.3, "cm"))) +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Biplot") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  geom_text_repel(aes(label = rownames(var1_pca3_load3), vjust= 1, hjust = 0, size = 3)) +
  xlab("PC1 (47.4 %)") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("PC2 (16.0 %)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(expand = FALSE, xlim = c(-1.25,1.25), ylim = c(-1.25,1.25))

print(pca1.var1.biplot)

ggsave(filename = "pca1.var1.biplot.jpg", plot = pca1.var1.biplot, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 20, height = 20, units = "cm", dpi = 300)




#SCREEPLOT (EIGENVALUES) PCA 3

pca1.var1.screeplot <- ggplot(var1_pca3_eig3, aes(x = reorder(comp10,-value), y = value, fill = value, group = 1)) +
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
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.5, 41))

print(pca1.var1.screeplot)

ggsave(filename = "pca1.var1.screeplot.jpg", plot = pca1.var1.screeplot, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




variance <- data.frame(
  comp8 = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
  PTV = round(c(0.474,0.160,0.079,0.072,0.030,0.021,0.018,0.014)*100,1),
  CPTV = round(c(0.474,0.634, 0.713, 0.785, 0.815, 0.836, 0.854,0.867)*100,1),
  x = (c(1)))
str(variance)




#BARPLOT PVT PCA 3

pca1.var1.PTV <- ggplot(variance, aes(x = comp8, y = PTV, fill = PTV)) +
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
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.5, 51))

print(pca1.var1.PTV)

ggsave(filename = "pca1.var1.PTV.jpg", plot = pca1.var1.PTV, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#BARPLOT CPVT PCA 3

pca1.var1.CPTV <- ggplot(variance, aes(x = "", y = PTV, fill = comp8, group = PTV)) + #sort(PTV, decreasing = TRUE)/ sort(comp8, decreasing = TRUE)
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c( "#C29DE3",  "#85C5FF", "#95BFC5",  "#C1FFE1",  "#D7D0BB" ,"#FFE497",  "#FBC27D", "#F9A1A5") ) +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Cumulative Proportion of Total Variance") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(title="Principal \n Components")) +
  theme(panel.border = element_blank()) +
  geom_text(aes(label = sort(PTV, decreasing = TRUE)), size = 5, colour = "#000000", position = position_stack(vjust = 0.5)) +
  xlab(label = "Stacked Princiapl Components") +
  theme(axis.title.x = element_text(size = 16, margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("CPTV (%)") +
  theme(axis.title.y = element_text(size = 16, margin = margin(t = 0, r = 10, b = 0, l = 0)))

print(pca1.var1.CPTV)

ggsave(filename = "pca1.var1.CPTV.jpg", plot = pca1.var1.CPTV, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#BARPLOT COMMUNALITIES PC1-PC8

pca1.var1.communalities <- ggplot(communalities, aes(x = reorder(rownames(communalities), -`rowSums(cos2[1:73, 1:8])`), y = `rowSums(cos2[1:73, 1:8])`, fill = `rowSums(cos2[1:73, 1:8])`, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Communalities of Variables in PC1 to PC8") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  #  geom_text(aes(label = var1names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Communality") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 1))

print(pca1.var1.communalities)

ggsave(filename = "pca1.var1.communalities.jpg", plot = pca1.var1.communalities, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




PC1.contribution <- as.data.frame(cbind(var1names = row.names(contribution), PC1cont = contribution$PC1))
PC1.contribution$PC1cont <- as.numeric(as.character(PC1.contribution$PC1cont))
str(PC1.contribution)


PC2.contribution <- as.data.frame(cbind(var1names = row.names(contribution), PC2cont = contribution$PC2))
PC2.contribution$PC2cont <- as.numeric(as.character(PC2.contribution$PC2cont))
str(PC2.contribution)


PC3.contribution <- as.data.frame(cbind(var1names = row.names(contribution), PC3cont = contribution$PC3))
PC3.contribution$PC3cont <- as.numeric(as.character(PC3.contribution$PC3cont))
str(PC3.contribution)


PC4.contribution <- as.data.frame(cbind(var1names = row.names(contribution), PC4cont = contribution$PC4))
PC4.contribution$PC4cont <- as.numeric(as.character(PC4.contribution$PC4cont))
str(PC4.contribution)


PC5.contribution <- as.data.frame(cbind(var1names = row.names(contribution), PC5cont = contribution$PC5))
PC5.contribution$PC5cont <- as.numeric(as.character(PC5.contribution$PC5cont))
str(PC5.contribution)


PC6.contribution <- as.data.frame(cbind(var1names = row.names(contribution), PC6cont = contribution$PC6))
PC6.contribution$PC6cont <- as.numeric(as.character(PC6.contribution$PC6cont))
str(PC6.contribution)


PC7.contribution <- as.data.frame(cbind(var1names = row.names(contribution), PC7cont = contribution$PC7))
PC7.contribution$PC7cont <- as.numeric(as.character(PC7.contribution$PC7cont))
str(PC7.contribution)


PC8.contribution <- as.data.frame(cbind(var1names = row.names(contribution), PC8cont = contribution$PC8))
PC8.contribution$PC8cont <- as.numeric(as.character(PC8.contribution$PC8cont))
str(PC8.contribution)




#PLOT VARIABLE CONTRIBUTION PC1

pca1.var1.PC1cont <- ggplot(PC1.contribution, aes(x = reorder(var1names, -PC1cont), y = PC1cont, fill = PC1cont, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Proportion of Variable Contribution to PC1") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
#  geom_text(aes(label = var1names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Contribution (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 3))

print(pca1.var1.PC1cont)

ggsave(filename = "pca1.var1.PC1cont.jpg", plot = pca1.var1.PC1cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC2

pca1.var1.PC2cont <- ggplot(PC2.contribution, aes(x = reorder(var1names, -PC2cont), y = PC2cont, fill = PC2cont, group = 1)) +
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
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 10))

print(pca1.var1.PC2cont)

ggsave(filename = "pca1.var1.PC2cont.jpg", plot = pca1.var1.PC2cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC3

pca1.var1.PC3cont <- ggplot(PC3.contribution, aes(x = reorder(var1names, -PC3cont), y = PC3cont, fill = PC3cont, group = 1)) +
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
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 15))

print(pca1.var1.PC3cont)

ggsave(filename = "pca1.var1.PC3cont.jpg", plot = pca1.var1.PC3cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC4

pca1.var1.PC4cont <- ggplot(PC4.contribution, aes(x = reorder(var1names, -PC4cont), y = PC4cont, fill = PC4cont, group = 1)) +
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
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 10))

print(pca1.var1.PC4cont)

ggsave(filename = "pca1.var1.PC4cont.jpg", plot = pca1.var1.PC4cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC5

pca1.var1.PC5cont <- ggplot(PC5.contribution, aes(x = reorder(var1names, -PC5cont), y = PC5cont, fill = PC5cont, group = 1)) +
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
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 15))

print(pca1.var1.PC5cont)

ggsave(filename = "pca1.var1.PC5cont.jpg", plot = pca1.var1.PC5cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC6

pca1.var1.PC6cont <- ggplot(PC6.contribution, aes(x = reorder(var1names, -PC6cont), y = PC6cont, fill = PC6cont, group = 1)) +
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
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 20))

print(pca1.var1.PC6cont)

ggsave(filename = "pca1.var1.PC6cont.jpg", plot = pca1.var1.PC6cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC7

pca1.var1.PC7cont <- ggplot(PC7.contribution, aes(x = reorder(var1names, -PC7cont), y = PC7cont, fill = PC7cont, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Proportion of Variable Contribution to PC7") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  #  geom_text(aes(label = var1names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Contribution (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 15))

print(pca1.var1.PC7cont)

ggsave(filename = "pca1.var1.PC7cont.jpg", plot = pca1.var1.PC7cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC7

pca1.var1.PC7cont <- ggplot(PC7.contribution, aes(x = reorder(var1names, -PC7cont), y = PC7cont, fill = PC7cont, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Proportion of Variable Contribution to PC7") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  #  geom_text(aes(label = var1names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Contribution (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 15))

print(pca1.var1.PC7cont)

ggsave(filename = "pca1.var1.PC7cont.jpg", plot = pca1.var1.PC7cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#PLOT VARIABLE CONTRIBUTION PC8

pca1.var1.PC8cont <- ggplot(PC8.contribution, aes(x = reorder(var1names, -PC8cont), y = PC8cont, fill = PC8cont, group = 1)) +
  geom_bar(stat = "identity", width = 0.9) +
  scale_fill_steps(high = "#4D868E", low = "#4D868E") +
  theme_bw() + 
  theme(text = element_text(size = 16)) +
  ggtitle("Proportion of Variable Contribution to PC8") +
  theme(legend.position = "none") +
  theme(panel.border = element_blank()) +
  #  geom_text(aes(label = var1names), vjust= -0.5, size = 3) +
  xlab("Variables") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  ylab("Contribution (%)") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.1, size = 10, angle = 90)) +
  coord_cartesian(expand = FALSE, xlim = c(0.5,NA), ylim = c(-0.1, 10))

print(pca1.var1.PC8cont)

ggsave(filename = "pca1.var1.PC8cont.jpg", plot = pca1.var1.PC8cont, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 30, height = 15, units = "cm", dpi = 300)




#bind with municipality names
var1_scores <- bind_cols(Municipio=municipio, var1_scores)

write.csv(var1_scores, file = "C:/Users/Antonia/Desktop/Documents/R/master/var1_scores_muninames.csv")


#read shapefiles
shp_muni <- read_sf("Cont_Mun_CAOP2023.shp")

#delete non useful columns
shp_muni <- shp_muni[,c(-1:-2,-4:-12)]

#data2 = standardized variables of health outcomes merged with shape file

#merge scaled variable data frame with shapefile
pca3.sfdata <- merge(shp_muni,var1_scores, by="Municipio")

#delete first row, which is of type "character"
#data2 <- data2[,-1]
#str(data2)




#PLOT SCORE MAP PC1

pca1.var1.PC1score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca3.sfdata, aes(fill = PC1)) +
  scale_fill_gradient2(low = "#F99D31", high = "#00447F", mid = "#ffffff", 
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

ggsave(filename = "pca1.var1.PC1score.jpg", plot = pca1.var1.PC1score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC2

pca1.var1.PC2score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca3.sfdata, aes(fill = PC2)) +
  scale_fill_gradient2(low = "#F99D31", high = "#00447F", mid = "#ffffff", 
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

ggsave(filename = "pca1.var1.PC2score.jpg", plot = pca1.var1.PC2score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC3

pca1.var1.PC3score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca3.sfdata, aes(fill = PC3)) +
  scale_fill_gradient2(low = "#F99D31", high = "#00447F", mid = "#ffffff", 
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

ggsave(filename = "pca1.var1.PC3score.jpg", plot = pca1.var1.PC3score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC4

pca1.var1.PC4score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca3.sfdata, aes(fill = PC4)) +
  scale_fill_gradient2(low = "#F99D31", high = "#00447F", mid = "#ffffff", 
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

ggsave(filename = "pca1.var1.PC4score.jpg", plot = pca1.var1.PC4score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC5

pca1.var1.PC5score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca3.sfdata, aes(fill = PC5)) +
  scale_fill_gradient2(low = "#F99D31", high = "#00447F", mid = "#ffffff", 
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

ggsave(filename = "pca1.var1.PC5score.jpg", plot = pca1.var1.PC5score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC6

pca1.var1.PC6score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca3.sfdata, aes(fill = PC6)) +
  scale_fill_gradient2(low = "#F99D31", high = "#00447F", mid = "#ffffff", 
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

ggsave(filename = "pca1.var1.PC6score.jpg", plot = pca1.var1.PC6score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC7

pca1.var1.PC7score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca3.sfdata, aes(fill = PC7)) +
  scale_fill_gradient2(low = "#F99D31", high = "#00447F", mid = "#ffffff", 
                       midpoint = 0) + #limit = c(-1,10), space = "Lab"
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("PC7") +
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

ggsave(filename = "pca1.var1.PC7score.jpg", plot = pca1.var1.PC7score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT SCORE MAP PC8

pca1.var1.PC8score <- ggplot(data = shp_muni) +
  geom_sf(color = "black", fill = "grey") +
  geom_sf(data = pca3.sfdata, aes(fill = PC8)) +
  scale_fill_gradient2(low = "#F99D31", high = "#00447F", mid = "#ffffff", 
                       midpoint = 0) + #limit = c(-1,10), space = "Lab"
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(text = element_text(size = 20)) +
  ggtitle("PC8") +
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

ggsave(filename = "pca1.var1.PC8score.jpg", plot = pca1.var1.PC8score, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#PLOT LOADING AND SCORE POSITIONS

plot.scores1 <- as.data.frame(cbind(municipio = municipio, scores = var1_scores$PC1))
plot.scores1$scores <- as.numeric(as.character(plot.scores1$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-150, 150)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca1.var1.PC1pos <- ggplot() +
  geom_point(data = var1_pca3_load3, aes(x = 1, y = PC1)) + 
  geom_text_repel(data = var1_pca3_load3, aes(1, PC1, label = row.names(var1_pca3_load3)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
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

print(pca1.var1.PC1pos)

ggsave(filename = "pca1.var1.PC1pos.jpg", plot = pca1.var1.PC1pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores2 <- as.data.frame(cbind(municipio = municipio, scores = var1_scores$PC2))
plot.scores2$scores <- as.numeric(as.character(plot.scores2$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-50, 50)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca1.var1.PC2pos <- ggplot() +
  geom_point(data = var1_pca3_load3, aes(x = 1, y = PC2)) + 
  geom_text_repel(data = var1_pca3_load3, aes(1, PC2, label = row.names(var1_pca3_load3)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
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

print(pca1.var1.PC2pos)

ggsave(filename = "pca1.var1.PC2pos.jpg", plot = pca1.var1.PC2pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores3 <- as.data.frame(cbind(municipio = municipio, scores = var1_scores$PC3))
plot.scores3$scores <- as.numeric(as.character(plot.scores3$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-40, 40)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca1.var1.PC3pos <- ggplot() +
  geom_point(data = var1_pca3_load3, aes(x = 1, y = PC3)) + 
  geom_text_repel(data = var1_pca3_load3, aes(1, PC3, label = row.names(var1_pca3_load3)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
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

print(pca1.var1.PC3pos)

ggsave(filename = "pca1.var1.PC3pos.jpg", plot = pca1.var1.PC3pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores4 <- as.data.frame(cbind(municipio = municipio, scores = var1_scores$PC4))
plot.scores4$scores <- as.numeric(as.character(plot.scores4$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-20, 20)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca1.var1.PC4pos <- ggplot() +
  geom_point(data = var1_pca3_load3, aes(x = 1, y = PC4)) + 
  geom_text_repel(data = var1_pca3_load3, aes(1, PC4, label = row.names(var1_pca3_load3)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
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

print(pca1.var1.PC4pos)

ggsave(filename = "pca1.var1.PC4pos.jpg", plot = pca1.var1.PC4pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores5 <- as.data.frame(cbind(municipio = municipio, scores = var1_scores$PC5))
plot.scores5$scores <- as.numeric(as.character(plot.scores5$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-10, 10)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca1.var1.PC5pos <- ggplot() +
  geom_point(data = var1_pca3_load3, aes(x = 1, y = PC5)) + 
  geom_text_repel(data = var1_pca3_load3, aes(1, PC5, label = row.names(var1_pca3_load3)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
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

print(pca1.var1.PC5pos)

ggsave(filename = "pca1.var1.PC5pos.jpg", plot = pca1.var1.PC5pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores6 <- as.data.frame(cbind(municipio = municipio, scores = var1_scores$PC6))
plot.scores6$scores <- as.numeric(as.character(plot.scores6$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-6, 6)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca1.var1.PC6pos <- ggplot() +
  geom_point(data = var1_pca3_load3, aes(x = 1, y = PC6)) + 
  geom_text_repel(data = var1_pca3_load3, aes(1, PC6, label = row.names(var1_pca3_load3)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
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

print(pca1.var1.PC6pos)

ggsave(filename = "pca1.var1.PC6pos.jpg", plot = pca1.var1.PC6pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores7 <- as.data.frame(cbind(municipio = municipio, scores = var1_scores$PC7))
plot.scores7$scores <- as.numeric(as.character(plot.scores7$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-6, 6)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca1.var1.PC7pos <- ggplot() +
  geom_point(data = var1_pca3_load3, aes(x = 1, y = PC7)) + 
  geom_text_repel(data = var1_pca3_load3, aes(1, PC7, label = row.names(var1_pca3_load3)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
  geom_point(data = plot.scores7, aes(x = 2, y = a+scores*b)) +
  geom_text_repel(data = plot.scores7, aes(2, a+scores*b, label = municipio), position = position_nudge_repel(x = 0.1, y = 0), size = 3) +
  scale_y_continuous("Loadings", sec.axis = sec_axis(~ (. - a)/b, name = "Scores")) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "grey", size = 0.25, linetype = 1),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
  ) +
  ggtitle("Position of variables and observations on PC7") +
  theme(plot.title = element_text(size = 18, hjust=0)) +
  xlab("PC1") +
  theme(axis.title.x = element_text(size = 14)) +
  ylab("Loadings") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(vjust = 1, size = 14, hjust = 1)) + 
  theme(axis.text.y = element_text(vjust = 1, size = 14, hjust = 1)) +
  coord_cartesian(expand = FALSE, xlim = c(0,3), ylim = c(-1.2,1.2))

print(pca1.var1.PC7pos)

ggsave(filename = "pca1.var1.PC7pos.jpg", plot = pca1.var1.PC7pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




plot.scores8 <- as.data.frame(cbind(municipio = municipio, scores = var1_scores$PC8))
plot.scores8$scores <- as.numeric(as.character(plot.scores8$scores))

ylim.prim <- c(-1.2, 1.2)   # in this example, Loadings
ylim.sec <- c(-6, 6)    # in this example, Scores

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]


pca1.var1.PC8pos <- ggplot() +
  geom_point(data = var1_pca3_load3, aes(x = 1, y = PC8)) + 
  geom_text_repel(data = var1_pca3_load3, aes(1, PC8, label = row.names(var1_pca3_load3)), position = position_nudge_repel(x = -0.1, y = 0), size = 3) +
  geom_point(data = plot.scores8, aes(x = 2, y = a+scores*b)) +
  geom_text_repel(data = plot.scores8, aes(2, a+scores*b, label = municipio), position = position_nudge_repel(x = 0.1, y = 0), size = 3) +
  scale_y_continuous("Loadings", sec.axis = sec_axis(~ (. - a)/b, name = "Scores")) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "grey", size = 0.25, linetype = 1),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
  ) +
  ggtitle("Position of variables and observations on PC8") +
  theme(plot.title = element_text(size = 18, hjust=0)) +
  xlab("PC1") +
  theme(axis.title.x = element_text(size = 14)) +
  ylab("Loadings") +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(vjust = 1, size = 14, hjust = 1)) + 
  theme(axis.text.y = element_text(vjust = 1, size = 14, hjust = 1)) +
  coord_cartesian(expand = FALSE, xlim = c(0,3), ylim = c(-1.2,1.2))

print(pca1.var1.PC8pos)

ggsave(filename = "pca1.var1.PC8pos.jpg", plot = pca1.var1.PC8pos, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 29.7, units = "cm", dpi = 300)




#CORRELATION MATRIX 2

var1_cormat2 <- round(cor(var1_2),2)

#get lower triangle of the correlation matrix
get_lower_tri2 <- function(var1_cormat2){
  var1_cormat2[upper.tri(var1_cormat2)] <- NA
  return(var1_cormat2)
}

#get upper triangle of the correlation matrix
get_upper_tri2 <- function(var1_cormat2){
  var1_cormat2[lower.tri(var1_cormat2)]<- NA
  return(var1_cormat2)
}

# Use correlation between variables as distance
reorder_cormat2 <- function(var1_cormat2){
  dd <- as.dist((1-var1_cormat2)/2)
  hc <- hclust(dd)
  var1_cormat2 <-var1_cormat2[hc$order, hc$order]
}

#reorder the correlation matrix
var1_cormat2 <- reorder_cormat2(var1_cormat2)
upper_tri2 <- get_upper_tri2(var1_cormat2)

#melt the correlation matrix
melt_cormat2 <- melt(upper_tri2, na.rm = TRUE)

# Create a ggheatmap
pca1.var1.cormat2 <- 
  ggplot(melt_cormat2, aes(Var2, Var1, fill = value))+
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
    panel.grid.major.y = element_line(color = "grey", size = 0.25, linetype = 1),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(0.5, 0),
    legend.position = "bottom",
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7, hjust = 1))+
  theme(axis.text.y = element_text(vjust = 0.5, size = 7, hjust = 1))+
  coord_fixed()


print(pca1.var1.cormat2)

ggsave(filename = "pca1.var1.cormat2.jpg", plot = pca1.var1.cormat2, device = "jpg", path = "C:/Users/Antonia/Desktop/Documents/R/master/plots/", width = 21, height = 26, units = "cm", dpi = 300)

write.csv(var1_cormat2, "C:/Users/Antonia/Desktop/Documents/R/master/var1_cormat2.csv", row.names=TRUE)
