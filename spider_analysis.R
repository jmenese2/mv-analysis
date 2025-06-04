
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(missMDA)
library(FactoMineR)
setwd("C:/Users/Owner/Documents/2025 spring/multi_analysis/NEwR-2ed_code_data/NEwR2-Functions")
source("cleanplot.pca.R")
source("PCA.newr.R")
source("CA.newr.R")
setwd("C:/Users/Owner/Documents/2025 spring/multi_analysis/exam_3")

spe <- read.csv("Spe.csv", row.names = 1,header=TRUE)
env <- read.csv("Environment.csv", row.names = 1, header=TRUE)
spe<- log1p(spe)

#3 PCA
#row 359
(spe.pca <- rda(spe))

# Scree plot and broken stick model
dev.new(
  title = "Scree plot of PCA eigenvalues - species data", 
  noRStudioGD = TRUE
)
screeplot(spe.pca,
          bstick = TRUE, 
          npcs = length(spe.pca$CA$eig)
)

# PCA biplots
spe.pca.sc1 <- scores(spe.pca, display = "species", scaling = 1)
spe.pca.sc2 <- scores(spe.pca, display = "species", scaling = 2)

dev.new(title = "PCA on spider species",
        width = 12,
        height = 6,
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
cleanplot.pca(spe.pca, scaling = 1, mar.percent = 0.06)
cleanplot.pca(spe.pca, scaling = 2, mar.percent = 0.06)

#5
#389
# A posteriori projection of environmental variables in a PCA
dev.new(
  title = "PCA biplot scaling 2 with environment",
  noRStudioGD = TRUE
)
biplot(spe.pca, main = "PCA spider abundances - scaling 2")
(spe.pca.env <-
    envfit(spe.pca, env, scaling = 2)) # Scaling 2 is default
# Plot significant variables with a user-selected colour
plot(spe.pca.env, p.max = 0.05, col = 3)

# My major conclusions by looking at the results of this posteriori
#projection of the environmental variables are listed below

# Illuminance was positively related with the species A. peri,A.fabr,A.acce

#herb cover and calamagrostis cover was positvely related with the species
#p. pull, p. nigr

#Water content was positively related with the species Z. spin & T. terr

#Lastly, tree cover, leaves and twigs were positively related with the species P. lugu




#3CA
#515

spe.ca <- cca(spe)
summary(spe.ca)

# Scree plot and broken stick model using vegan's screeplot.cca()
dev.new(title = "Scree plot of CA eigenvalues", noRStudioGD = TRUE)
screeplot(spe.ca, bstick = TRUE, npcs = length(spe.ca$CA$eig))

# CA biplots
dev.new(title = "CA biplots",
        width = 14,
        height = 7,
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
# Scaling 1: sites are centroids of species
plot(spe.ca, 
     scaling = 1, 
     main = "CA spider abundances - biplot scaling 1"
)
# Scaling 2 (default): species are centroids of sites
plot(spe.ca, main = "CA spider abundances - biplot scaling 2")

#3NMDS
#752


## NMDS applied to the Doubs spider species - percentage difference
## dissimilarity matrix

spe.nmds <- metaMDS(spe, distance = "bray")
spe.nmds
spe.nmds$stress
dev.new(title = "NMDS on spider species - Percentage difference",
        noRStudioGD = TRUE
)
plot(
  spe.nmds,
  type = "t",
  main = paste(
    "NMDS/Percentage difference - Stress =",
    round(spe.nmds$stress, 3)
  )
)

# Shepard plot and goodness of fit
dev.new(title = "NMDS - Shepard plot",
        width = 12,
        height = 6,
        noRStudioGD = TRUE
)
par(mfrow = c(1, 2))
stressplot(spe.nmds, main = "Shepard plot")
gof <- goodness(spe.nmds)
plot(spe.nmds, type = "t", main = "Goodness of fit")
points(spe.nmds, display = "sites", cex = gof * 300)

#3 answers

# I did not observe any major differences in terms of the results across
# different analyses. Most of the species were accurately shown to be related to
#certain sites across the analyses. for example the species pard.lugu
#was found primarily in sites 15,17,18, & 20 and this was represented across
#the different analyses. 

#4
#316
# Combining clustering and ordination results

# Clustering the objects using the environmental data: Euclidean
# distance after standardizing the variables, followed by Ward
# clustering
spe.w <- hclust(dist(scale(spe)), "ward.D")
# Cut the dendrogram to yield 4 groups
gr <- cutree(spe.w, k = 4)
grl <- levels(factor(gr))

# Extract the site scores, scaling 1
sit.sc1 <- scores(spe.pca, display = "wa", scaling = 1)

# Plot the sites with cluster symbols and colours (scaling 1)
dev.new(title = "Ordination and clustering", noRStudioGD = TRUE)
p <- plot(
  spe.pca,
  display = "wa",
  scaling = 1,
  type = "n",
  main = "PCA correlation + clusters"
)
abline(v = 0, lty = "dotted")
abline(h = 0, lty = "dotted")
for (i in 1:length(grl)) {
  points(sit.sc1[gr == i, ],
         pch = (14 + i),
         cex = 2,
         col = i + 1)
}
text(sit.sc1, row.names(spe), cex = 0.7, pos = 3)
# Add the dendrogram
ordicluster(p, spe.w, col = "dark grey")
# Add legend interactively
legend(
  locator(1),
  paste("Cluster", c(1:length(grl))),
  pch = 14 + c(1:length(grl)),
  col = 1 + c(1:length(grl)),
  pt.cex = 2
)



