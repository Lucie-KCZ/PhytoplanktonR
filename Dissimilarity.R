###########################
### Dissimilarity class ###
###  phytoplankton gp.  ###
###########################

# Lucie Kuczynski
# lucie.kuczynski@hotmail.com
# last edit: Nov 25, 2021

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

sorensen.metrics <- beta.pair(pa.matrix, index.family = 'sorensen')
bray.metrics <- bray.part(com.matrix)

# comparing metrics
par(pty = 's')
plot(sorensen.metrics$beta.sor, bray.metrics$bray, pch = 20, 
     xlab = 'Sorensen (P/A)', ylab = 'Bray Curtis (Abundances)', 
     main = 'Comparison of turnover metrics',
     xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, col = 'grey')

# looking at the different beta diversity components
hist(sorensen.metrics$beta.sor, col = 'grey', border = 'white',
     las = 1, xlab = 'Beta diversity', main = '', ylim = c(0, length(sorensen.metrics$beta.sor)))
hist(sorensen.metrics$beta.sim, col=rgb(0,0,1,1/4), add = T) 
hist(sorensen.metrics$beta.sne, col=rgb(1,0,0,1/4), add = T) 

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

gower.dist <- cluster::daisy(env[, -1], metric = 'gower')

plot(sorensen.metrics$beta.sim ~ gower.dist, pch = 20, las = 1, 
     xlab = 'Environmental distance', ylab = 'Turnover')

# Mantel test
install.packages('ape')
library(ape)
mantel.test(m1 = as.matrix(sorensen.metrics$beta.sim), 
            m2 = as.matrix(gower.dist), 
            nperm = 999, graph = F)

mantel.test(m1 = as.matrix(bray.metrics$bray), 
            m2 = as.matrix(gower.dist), 
            nperm = 999, graph = F)

# Variance Partitioning
library(vegan)
VP <- varpart(Y = bray.metrics$bray, 
              dist(env$elevation, method = 'euclidean'), 
              dist(env$water, method = 'binary'),
              dist(table(env$site, env$land.use), method = 'binary'))
showvarparts(2, bg = c("tomato3","deepskyblue2"))

# Functional dissimilarity ------------------------------------------------

traits$trophic <- as.factor(traits$trophic)
fal.dist <- cluster::daisy(traits[, -1], metric = 'gower')

# For the visualisation
install.packages('ape')
library(ape)

help(pcoa)
fal.pcoa <- pcoa(fal.dist)
rownames(fal.pcoa$vectors) <- traits$species
biplot(fal.pcoa)
