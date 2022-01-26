####################################
### Functional Diversity lecture ###
###       phytoplankton gp       ###
####################################

# Lucie Kuczynski
# lucie.kuczynski@hotmail.com
# last edit: Jan 14, 2022

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

# Univariate indices ------------------------------------------------------

# Computing relative abundances
# structure of the apply function
# apply(X, MARGIN, FUN )
#   X: what you will apply the function FUN on
#   MARGIN: the way you will: 
#           1 on each row individually, 
#           2 on each col individually.
#   FUN: the function you apply on X the way said with MARGIN
# Here: apply(aravo$spe, 1, sum) means:
# you apply the function sum on each row of aravo$spe
# basically, it's the total abundance of each site
Relative.abund <- aravo$spe / apply(aravo$spe, 1, sum)

# Here the use of apply is a little bit more complicated
# but basically, you apply the weighted mean function 
# on each row of Relative.abund
CWM.height <- apply(Relative.abund, 1, function(z) weighted.mean(x = aravo$traits$Height, w = z))
CWM.area <- apply(Relative.abund, 1, function(z) weighted.mean(x = aravo$traits$Area, w = z))

par(pty = 's', mfrow = c(2, 2))
hist(aravo$traits$Height, col = 'forestgreen', border = 'white', xlab = '', ylab = '', las = 1)
hist(CWM.height, col = 'forestgreen', border = 'white', xlab = '', ylab = '', las = 1)
hist(aravo$traits$Area, col = 'darkmagenta', border = 'white', xlab = '', ylab = '', las = 1)
hist(CWM.area, col = 'darkmagenta', border = 'white', xlab = '', ylab = '', las = 1)

rm(CWM.area, CWM.height)

# Multivariate indices ----------------------------------------------------
# We will need three packages
library(fundiversity)
library(cluster)
library(ape)

# useful vignette for fundiversity:
# https://cran.r-project.org/web/packages/fundiversity/vignettes/fundiversity.html

# there's also a very recent package: mFD
# https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.05904
# and the associated vignette/tutorial:
# https://cmlmagneville.github.io/mFD/articles/mFD_general_workflow.html

# running the distance matrix using a Gower distance
traits.Gower <- cluster::daisy(x = aravo$traits, metric = 'gower')

# running the PCoA itself
traits.PCoA <- ape::pcoa(traits.Gower)

# to know the amount of variance captured by the PCoA axes
round(100 * (traits.PCoA$values$Eigenvalues)/sum(traits.PCoA$values$Eigenvalues))

# based on the eigenvalues, four axes makes sense
traits.PCoA <- traits.PCoA$vectors[, 1:4]

# The second step is to simply run the functions! :)
functional_indices <- data.frame(
  FRic = fundiversity::fd_fric(traits = traits.PCoA, sp_com = aravo$spe)[, 2], 
  FEve = fundiversity::fd_feve(traits = traits.PCoA, sp_com = aravo$spe)[, 2],
  FDiv = fundiversity::fd_fdiv(traits = traits.PCoA, sp_com = aravo$spe)[, 2],
  Rao = fundiversity::fd_raoq(dist_matrix = traits.Gower, 
                              sp_com = as.matrix(aravo$spe))[, 2])

# to quickly look at the relationship between indices
library(GGally)
GGally::ggpairs(functional_indices, title = 'Functional indices') 







