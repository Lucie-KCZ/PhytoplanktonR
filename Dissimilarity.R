###########################
### Dissimilarity class ###
###  phytoplankton gp.  ###
###########################

# Lucie Kuczynski
# lucie.kuczynski@hotmail.com
# last edit: Nov 25, 2021

# corresponding recording can be found here:
# https://osf.io/8uqcv/?view_only=03935ee3d6ac4567aa3180c6cdeb71ef

# DONT RUN IF YOU HAVE THINGS IN YOUR ENVIRONMENT YOU WANT TO KEEP!
# cleaning the environment and the memory
for(i in 1:10) gc(reset = T) ; rm(list = ls())

# Data --------------------------------------------------------------------
# For Section Beta div., we just need a community matrix
# A community matrix is a table with sites as rows (i.e., statistical individuals/observations)
# and species as columns; it can be presence/absence or abundance

n.sites <- 5

# we are creating a random community matrix 
com.matrix <- matrix(data = rbinom(n.sites*10, 60, 0.05), 
                     nrow = n.sites, 
                     dimnames = list(paste('site', 1:n.sites, sep = '_'), 
                                     paste('sp', 1:10, sep = '_')))

# let's also have it as a presence-absence matrix
pa.matrix <- ifelse(com.matrix > 0, 1, 0)

# For Section Environmental gradient, we need environmental variables
# Again, sites are the rows, and we will have three environmental variables
# again, we are creating random variables 
# for the different environmental variables
env <- data.frame(site = paste('site', 1:n.sites, sep = '_'), 
                  elevation = round(runif(n.sites, min = 10, max = 50), 1),
                  temp.anomalies = rnorm(n = n.sites, mean = 5, sd = 2), 
                  water = rbinom(n = n.sites, size = 1, prob = .3), 
                  land.use = as.factor(sample(LETTERS[1:3], n.sites, replace = T)))
# let's check it out
summary(env)

# Finally, for Section Functional dissimilarity,
# we need a species trait matrix
# in which species are rows and traits are columns
traits <- data.frame(species = paste('sp', 1:10, sep = '_'), 
                     size = rnorm(n = 10, mean = 25, sd = 2), 
                     longevity = rpois(n = 10, lambda = 10), 
                     optimal.temp = rnorm(n = 10, mean = 7, sd = 2),
                     trophic = c('herb', 'herb', 'omn', 'pred', 'omn', 'herb', 'pred', 'omn', 'herb', 'pred'))
# let's have a look at it
summary(traits)

# Beta diversity ----------------------------------------------------------
# 1. vegan
# install the package
install.packages('vegan')
# load the package
library(vegan)

# quick off topic - to compute alpha diversity measures:
help(diversity)
# calculate Simpson's 1-D Index of Diversity for each site. 
# closer to 1 = greater diversity
diversity(com.matrix, index = 'simpson') 

# for beta diversity:
help(vegdist)
# Bray-Curtis is good for abundances
# and not considering joint absences
# puts the emphasis on compostion and 
# relative abundances differences
(bray <- vegdist(com.matrix, 'bray'))

# 2. betapart
# install the package
install.packages('betapart')
# load the package
library(betapart)

# Let's compute the Sorensen metric first which is
# based on presence absence: we will thus use the pa.matrix
sorensen.metrics <- beta.pair(pa.matrix, index.family = 'sorensen')

# Here is the Bray Curtis metric
bray.metrics <- bray.part(com.matrix)

# comparing metrics
# get the plot square, you don't have to run it if you don't like it! 
# here it is nice because the two axes will both go from 0 to 1, 
# so a square helps somehow with the visualisation
par(pty = 's') 

# here we create the biplot
# time to get (not too) crazy with colors and shapes and size to make a nice plot
# NOTE: the '$' only works if you are using it with a dataframe or a list
# Check class(sorensen.metrics): it's a list
plot(sorensen.metrics$beta.sor, bray.metrics$bray, pch = 20, 
     xlab = 'Sorensen (P/A)', ylab = 'Bray Curtis (Abundances)', 
     main = 'Comparison of turnover metrics',
     xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, col = 'grey')

# looking at the different beta diversity components
# Here is the global beta diversity (i.e., Sorensen distance)
hist(sorensen.metrics$beta.sor, col = 'grey', border = 'white',
     las = 1, xlab = 'Beta diversity', main = '', ylim = c(0, length(sorensen.metrics$beta.sor)))
# here we add (because of the 'add = T' argument) in blue the turnover component (i.e., Simpson index)
hist(sorensen.metrics$beta.sim, col=rgb(0,0,1,1/4), add = T) 
# and here we add the nestedness (in red)
hist(sorensen.metrics$beta.sne, col=rgb(1,0,0,1/4), add = T) 
# clean the plot area
dev.off()

# Environmental gradient --------------------------------------------------
# we don't need to install any packages for 
# common distances such as the euclidean distance
dist(env$elevation, method = 'euclidean')
dist(env$water, method = 'binary')
dist(table(env$site, env$land.use), method = 'binary')

# Gower's distance
install.packages('cluster')
library(cluster)

# Here we compute the distance, 
# but we don't want the site ID to be used 
# that's why we remove the first col.
# the notation package::function is sometimes used when 
# packages have functions with the same name or
# when you want to remember which package you need for some functions 
# (useful when you have to write the methods and cite the packages)
gower.dist <- cluster::daisy(env[, -1], metric = 'gower')

plot(sorensen.metrics$beta.sim ~ gower.dist, pch = 21, las = 1, cex = 1.5,
     xlab = 'Environmental distance', ylab = 'Turnover', bg = 'tomato3')
dev.off()

# Mantel test
install.packages('ape')
library(ape)
# the mantel test tests whether there is a significant correlation
# between two distance matrices (m1 and m2)
# using randomisation (nperm = 999) to estimate the significance

# H0: the correlation between m1 and m2 is equal to zero
# H1: the correlation between m1 and m2 is different from zero

# a good overview: https://mb3is.megx.net/gustame/hypothesis-tests/the-mantel-test

# here we do it with the turnover based on occurrences (the Simpon index)
mantel.test(m1 = as.matrix(sorensen.metrics$beta.sim), 
            m2 = as.matrix(gower.dist), 
            nperm = 999, graph = F)

# and here the one with the Bray Curtis dissimilarity
mantel.test(m1 = as.matrix(bray.metrics$bray), 
            m2 = as.matrix(gower.dist), 
            nperm = 999, graph = F)

# Variance Partitioning
library(vegan)

# the variance partitioning look how much variance in turnover (Y = bray.metrics$bray)
# can be explained by other variables (the three distance matrices)
# a good resource:
# https://mb3is.megx.net/gustame/constrained-analyses/variation-partitioning 
VP <- varpart(Y = bray.metrics$bray, 
              dist(env$elevation, method = 'euclidean'), 
              dist(env$water, method = 'binary'),
              dist(table(env$site, env$land.use), method = 'binary'))

# Functional dissimilarity ------------------------------------------------
# Here it is more an coarse first step really rather than  going in the rabbit hole of functional diversity

# For the daisy function, you need factors rather than character
# you thus can transform all character-variable into
# factor-variable by doing as done the next line
# (you'll have to do it for each individual variable though!)
traits$trophic <- as.factor(traits$trophic)

# We don't want the species names included in the distance computation
# hense the [, -1] after traits
fal.dist <- cluster::daisy(traits[, -1], metric = 'gower')

# For the visualisation
install.packages('ape')
library(ape)

# Here we are running a PCoA
# We'll talk about it when we'll talk about functional diversity
# let's have it as a black box for now
help(pcoa)
fal.pcoa <- pcoa(fal.dist)
rownames(fal.pcoa$vectors) <- traits$species
# The biplot shows you how similar species are
# the closer, the more similar
biplot(fal.pcoa)
