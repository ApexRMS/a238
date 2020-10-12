library(prioritizr)
library(raster)
library(sp)
library(rasterVis)
library(viridis)
library(ggplot2)
library(vegan)
library(gstat)
#library("lpsymphony")
library(Rsymphony)
set.seed(10)

# create structure
xy <- expand.grid(1:100, 1:100)
names(xy) <- c("x","y")
 
# define the gstat object (spatial model)
g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=vgm(psill=0.025, model="Exp", range=2), nmax=5)
yy <- predict(g.dummy, newdata=xy, nsim=1)
yy[ , 3] <- (1-decostand(yy[, 3], "range"))^2
#yy[, 3] <- ifelse(yy[, 3] < 0.3, 0, yy[, 3]) 
gridded(yy) = ~x+y
yyr <- raster(yy)

g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=vgm(psill=0.025, model="Exp", range=5), nmax=5)
yy1 <- predict(g.dummy1, newdata=xy, nsim=1)
yy1[, 3] <- decostand(yy1[, 3], "range")^1.5 
#yy1[, 3] <- ifelse(yy1[, 3] < 0.4, 0, yy1[, 3]) 
gridded(yy1) = ~x+y
yyr1 <- raster(yy1)

g.dummy2 <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=vgm(psill=0.025, model="Exp", range=20), nmax=5)
yy2 <- predict(g.dummy2, newdata=xy, nsim=1)
yy2[, 3] <- decostand(yy2[, 3], "range")
#yy2[, 3] <- ifelse(yy2[, 3] <0, 0, yy2[, 3])
gridded(yy2) = ~x+y
yyr2 <- raster(yy2)

#
costlayer <- yyr
values(costlayer) <- 1
speciesHabitatSuit <- stack(yyr, yyr1, yyr2)
names(speciesHabitatSuit) <- c("Sp1", "Sp2", "Sp3")
spplot(speciesHabitatSuit, main="Habitat suitability", ncol=3)

##Scenarios where you minimze the cost of the relative feature targets
minprob30 <- problem(costlayer, speciesHabitatSuit) %>% #input is the cost surface + features 
      add_min_set_objective() %>% #minimize cost surface
      add_relative_targets(0.35) %>% 
      add_binary_decisions() #inclusion vs no-inclusion	
minSol30 <- solve(minprob30)
plot(minSol30)
feature_representation(minprob30, minSol30)

fr30 <- freq(minSol30)[2, "count"]/c(freq(minSol30)[1, "count"] + freq(minSol30)[2, "count"])
fr30

## Scenarios which maximize feature representation conditioned on budget


maxfeat30 <- problem(costlayer, speciesHabitatSuit) %>% #input is the cost surface + features 
      add_max_features_objective(0.225 * 100 * 100) %>% #minimize cost surface
	  add_relative_targets(0.35) %>% 
      add_binary_decisions() #inclusion vs no-inclusion	
maxfeatsol30 <- solve(maxfeat30)
plot(maxfeatsol30)
feature_representation(maxfeat30, maxfeatsol30)
frmaxFeat30 <- freq(maxfeatsol30)[2, "count"]/c(freq(maxfeatsol30)[1, "count"] + freq(maxfeatsol30)[2, "count"])

## Scenarios which maximizes coverage across all features given the budget
  #  as much of the features as possible without exceeding a budget.
maxUtility30 <- problem(costlayer, speciesHabitatSuit) %>% #input is the cost surface + features 
      add_max_utility_objective(0.225 * 100 * 100) %>% #minimize cost surface
      add_binary_decisions() #inclusion vs no-inclusion	
maxUtilitysol30 <- solve(maxUtility30)
plot(maxUtilitysol30)
feature_representation(maxUtility30, maxUtilitysol30)
frmaxUtil30 <- freq(maxUtilitysol30)[2, "count"]/c(freq(maxUtilitysol30)[1, "count"] + freq(maxUtilitysol30)[2, "count"])



## Scenarios which minimize shortfalls
minShortfall30 <- problem(costlayer, speciesHabitatSuit) %>% #input is the cost surface + features 
      add_min_shortfall_objective(0.225 * 100 * 100) %>% #minimize cost surface
      add_relative_targets(0.35) %>% 
      add_binary_decisions() #inclusion vs no-inclusion	
minShortfallSol30 <- solve(minShortfall30)
plot(minShortfallSol30)
feature_representation(minShortfall30, minShortfallSol30)
frmaxShort30 <- freq(minShortfallSol30)[2, "count"]/c(freq(minShortfallSol30)[1, "count"] + freq(minShortfallSol30)[2, "count"])



spplot(stack(minSol30, maxfeatsol30, maxUtilitysol30, minShortfallSol30))
plot( minSol30 - minShortfallSol30)
plot(maxUtilitysol30-minShortfallSol30)
plot(maxfeatsol30-minShortfallSol30)