############################
### Multivariate lecture ###
###   phytoplankton gp   ###
############################

# Lucie Kuczynski
# lucie.kuczynski@hotmail.com
# last edit: Jan 17, 2022

# DONT RUN IF YOU HAVE THINGS IN YOUR ENVIRONMENT YOU WANT TO KEEP!
# cleaning the environment and the memory
for(i in 1:10) gc(reset = T) ; rm(list = ls())

# Data --------------------------------------------------------------------

# importing data from the library ade4
library(ade4)
data(aravo)

# This dataset describe the distribution of 82 species of Alpine plants in 75 sites. 
# Choler, P. (2005) Consistent shifts in Alpine plant traits along a mesotopographical gradient. 
#     Arctic, Antarctic, and Alpine Research, 37,444–453.

help(aravo)
summary(aravo)

# community matrix
str(aravo$spe)
# environmental data
str(aravo$env)
# trait data
str(aravo$traits)

# PCA ---------------------------------------------------------------------
library(FactoMineR)
res.pca <- PCA(aravo$env, scale.unit = TRUE, graph = T)

res.pca <- PCA(aravo$env[, -c(3, 5)], scale.unit = TRUE, graph = T)
print.PCA(res.pca)

# eigenvalues to know how much of the variance is captures
# by the newly created dimensions
res.pca$eig

# interesting to look at the correlation, 
# especially if a very small arrow on the plot axis1:axis2
res.pca$var$cor

# the values for each statistical individual
res.pca$ind$coord
plot.PCA(res.pca)
plot(res.pca, choix = 'ind', habillage = 1)
plot(res.pca, choix = 'ind', habillage = 2)
plot(res.pca, choix = 'ind', habillage = 3)

plot(res.pca, axes = c(1, 3), choix = 'ind', habillage = 3)
plot(res.pca, axes = c(2, 3), choix = 'ind', habillage = 3)

dev.off() ; rm(res.pca)

# RDA & CCA ---------------------------------------------------------------
# If you expected linear responses, use RDA. 
# If you expect unimodal responses you should chose CCA. 

library(vegan)

res.rda <- rda(aravo$spe ~ as.matrix(aravo$env[, -c(3, 5)]))

# sp scores are red
# sites scores are black
# env variable contributions are the blue arrows
plot(res.rda)

# note that because it's constrained, 
# the env variables are not 'grouped' the same way
# as in the PCA

# the code for the cca follows the exact same syntax :)
# res.cca <- cca(aravo$spe ~ as.matrix(aravo$env[, -c(3, 5)]))
# plot(res.cca)

dev.off() ; rm(res.rda)

# Co-inertia --------------------------------------------------------------
# reminder: https://pbil.univ-lyon1.fr/R/pdf/course6.pdf 
# for a nice tutorial

library(ade4)

# first step: run the two pcas
# # here it's done on sp abundances and their traits
# dudi1 <- dudi.pca(t(aravo$spe), scale = TRUE, scan = TRUE)
# dudi2 <- dudi.pca(aravo$traits, scale = TRUE, scan = FALSE, nf = 3)

# here it's done on sp abundances and environmental conditions at sites
dudi1 <- dudi.pca(aravo$spe, scale = TRUE, scan = TRUE)
dudi2 <- dudi.pca(aravo$env[, -c(3, 5)], scale = TRUE, scan = FALSE, nf = 3)

# then, run the coinertie
coin1 <- coinertia(dudi1, dudi2, scan = TRUE)

# look at the results
coin1
summary(coin1)

# the RV coefficent varies from -1 to 1
coin1$RV

# plots of the two PCAs
# clab = XXX ; tag size
# com matrix based PCA
s.arrow(coin1$c1, clab = 0.7)
# env based PCA
s.arrow(coin1$l1, clab = 0.7)


par(mfrow = c(1, 2))
# sp abundance based PCA in the coinertia space
s.corcircle(coin1$aX)
# env based PCA in the coinertia space
s.corcircle(coin1$aY)

par(mfrow = c(1, 1))
plot(coin1)

dev.off() ; rm(coin1, dudi1, dudi2)

# 4th corner --------------------------------------------------------------
library(ade4)

res.4th <- fourthcorner(aravo$env, 
             aravo$spe,
             aravo$traits, 
             modeltype = 6, p.adjust.method.G = 'fdr',
             p.adjust.method.D = 'fdr', nrepet = 999)

# 'fdr' is a way to correct for multiple comparisons
# 'nrepet' needs to be fairly high to have enough power in corrected tests


plot(res.4th, alpha = .05, stat = "D2")

# Blue cells correspond to negative significant relationships and
# red cells correspond to positive significant relationships. 

# 'alpha' is a the significance level

# In this example, there are some associations between categorical traits and 
# quantitative environmental variables which can be measured in three different ways. 
# These methods correspond to three possible values of the 'stat' argument in the plot and print functions:
  
# stat = 'D2': the association is measured between the quantitative variable and 
#               each category seperately. A correlation coefficient is used to indicate 
#               the strength of the association between the given category and 
#               the small or large values of the quantitative variables.
# stat = 'G': the association between the quantitative variable and the whole categorical variable 
#             is measured by a global statistic (F).
# stat = 'D': the association is estimated between the quantitative variable and 
#               each category separately by a measure of the within-group homogeneity. 
#               The strength of the association is indicated by the dispersion of the values 
#               of the quantitative variable for a given category.

dev.off() ; rm(res.4th)





