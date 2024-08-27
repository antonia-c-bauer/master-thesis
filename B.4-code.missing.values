#Universidade de Coimbra - master thesis - missing values


library(base)




#avoid scientific notation
options(prompt="R> ", scipen=999)

#set working directory
setwd("C:/Users/Antonia/Desktop/Documents/R/master/")

#load data
variable2 <- read.csv("health outcomes data CID R.csv")
municipio <- variable2$name

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
