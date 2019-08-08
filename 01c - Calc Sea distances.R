# R code from the manuscript

# Title: "A multi-locus demo-genetic simulator to model local adaptation in heterogeneous landscapes"
# Authors: Marco Andrello, Christelle Noirot, Florence Débarre and Stéphanie Manel
# Submitted to: Molecular Ecology Resources
# Date: 05 Aug 2019

# "01c - Calc distancesR"
# Code to calculate least-cost marine distances and salinity differences
# between demes in the example on Mullus surmuletus
# Marco Andrello
# 05/08/2019


rm(list=ls())

library(raster)
library(rgdal)
library(rgeos)
library(gdistance)



#########################################################################
# calculating least-cost marine distances
#########################################################################

# Load the raster of marine areas
# This should be already a conductance layer: Sea = 1, others (land) = 0
SeaRaster <- raster("Medi.tif")

# Warning: these calculation can take some minutes to complete!
# Create transition matrices based on conductance rasters
t <- transition(SeaRaster, function(x) mean(x), directions=8) 
# Correct it for lat-long
t <- geoCorrection(t, type="c")
save(t,file="Transition_Medi.RData")

# Extract coordinates of the centroids of the demes
deme <- readOGR("Demes.shp")
# p <- gCentroid(deme,byid=T)
# write.csv(p,file="Centroids.csv")
# writeOGR(p,"centroids.shp",driver="ESRI Shapefile")
# # Shift centroids falling in land into water manually in QGis
# rm(p)
p <- readOGR("Centroids.shp")
# points(p@coords[,1],p@coords[,2])

# Calculating least-cost distances
# Warning: this calculation can take some minutes to complete!
sea_dist <- costDistance(t, p@coords) 
sea_dist <- sea_dist / 1000
sea_dist <- as.matrix(sea_dist)


#########################################################################
# calculating salinity differences
#########################################################################

load("Salinity_data.RData")
salinity_diff <- dist(Salinity_per_deme[3,])
salinity_diff <- as.matrix(salinity_diff)


#########################################################################
# Save
#########################################################################
save(sea_dist,
     salinity_diff,
     file="Distances.RData")

