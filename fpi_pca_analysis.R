setwd("C:/Users/Owner/Documents/2025 spring/multi_analysis/exam1")
#chp.5 Load packages, functions and data ===============================
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(missMDA)
library(FactoMineR)
# Source additional functions that will be used later in this
# Chapter. Our scripts assume that files to be read are in
# the working directory.
setwd("C:/Users/Owner/Documents/2025 spring/multi_analysis/NEwR-2ed_code_data/NEwR2-Functions")
source("cleanplot.pca.R")
source("PCA.newr.R")
source("CA.newr.R")
setwd("C:/Users/Owner/Documents/2025 spring/multi_analysis/exam1")
fpi <- read.csv("FPI.csv", row.names = 1,
                header=TRUE)
#PCA analysis
# A reminder of the content of the env dataset
summary(fpi)

# PCA based on a correlation matrix
# Argument scale=TRUE calls for a standardization of the variables
fpi.pca <- rda(fpi, scale = TRUE)
fpi.pca
summary(fpi.pca) # Default scaling 2
summary(fpi.pca, scaling = 1)

# Eigenvalues
(ev_fpi <- fpi.pca$CA$eig)

# Scree plot and broken stick model
dev.new(title = "Scree plot of PCA eigenvalues", noRStudioGD = TRUE)
screeplot(fpi.pca, bstick = TRUE, npcs = length(fpi.pca$CA$eig))

# Plots using cleanplot.pca
# A rectangular graphic window is needed to draw the plots together
dev.new(width = 12,
        height = 6,
        title = "PCA biplots - environmental variables - cleanplot.pca", 
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(fpi.pca, scaling = 1, mar.percent = 0.08)
cleanplot.pca(fpi.pca, scaling = 2, mar.percent = 0.04)



2.724891+0.2751086
#PCA1&2=2.724891
#total= 3
(2.724891/3)*100
#90.8297%

#2 questions
#a. the proportion of the variance accounted for by the first two principal components
#(PC1 and PC2) is 90.8297%

#b. According to the broken stick model only 1 principal axis is substantial
# and that is PC1

#c. The estimated loading scores for each variable are listed below 
#  Com: PC1:-2.2  PC2: 2.6
#  Econ: PC1:-2.5 PC2: 0
#  Ecol: PC1: -2  PC2: -2.5
#Ecol variable has the strongest corelation with PC1 and the Com variable has
#the stongest correlation with PC2

#d. According to the circle of equilibrium the Com indicator can be
#interpreted with higher confidence

#e. The direction of each indicator
#will rank the highest to lowest performance for individual fisheries
#the results showed that
#Econ: best 3: A1,A2,B14   worst 3: B19,B37,B35
#Com:  best 3: B14,A8,A9   Worst 3: B19,B37,B35
#Ecol: best 3: A1,A2,A4    Worst 3: B39,B40,B33
# the results do coincide with the results of the TBL indicator-fishery table.

#f. Based on the small 45 degree angles between the indicators there seems to be strong
#correlation between the TBL indicators

#g. Based on the results of the PCA analysis we can look at fisheries
# inside the angles that are shared between each indicator, This area would share
#fisheries that have relatively high scores in all three of the TBL indicators, meaning
#this will be our "sustainability quadrant". Using this sustainability
#quadrant the current state of a fishery can be assessed 
#and areas where improvements are needed can be identified.
  