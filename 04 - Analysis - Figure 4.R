# R code from the manuscript

# Title: "A multi-locus demo-genetic simulator to model local adaptation in heterogeneous landscapes"
# Authors: Marco Andrello, Christelle Noirot, Florence Débarre and Stéphanie Manel
# Submitted to: Molecular Ecology Resources
# Date: 05 Aug 2019

# "04 - Analysis - Figure 4.R"
# Analyse the patterns of adaptation in the Eastern Gulf of Lion and Dardanelles
# of the demo-genetic simulations on Mullus surmuletus
# Patterns in Local adaptation
# Marco Andrello
# 05/08/2019

rm(list=ls())

# library(maps)
# library(mapdata)
library(RColorBrewer)
library(rgdal)

# # Load only useful data
load("Results_base.RData")
sigma_ij <- t(surv.trans) # NB : sigma_ij is the survival of individuals originating in i (row) transplanted in j (column)
rm(fst_A, fst_B, fst_C, fst_D, Globfst, Nt, Nt_A, Nt_B, Nt_C, Nt_D, ro, surv.trans)

deme <- readOGR(dsn="Demes.shp") 

load("Distances.RData")


# Calculating mean survival in each deme
sigma_i <- diag(sigma_ij)                               

# Calculating local adaptation LAij for each pair of demes
LA_ij <- matrix(NA, nrow=100, ncol=100)                 
for (i in 1 : 100) {
  for (j in 1 : 100) {
    LA_ij[i,j] <- sigma_i[i] - sigma_ij[i,j]
  }
}
rm(i,j,AL)


# Create data frame with local adaptation, sea distance and salinity distance
# ("as.vector" transforms a matrix into a vector column by column)
# Local adaptation of individuals originating in i and translocated into j
LA_vec <- as.vector(LA_ij)
# Spatial marine distance in km between i and j
sea_dist_vec <- as.vector(as.matrix(sea_dist)) 
# Salinity difference in km between i and j
salinity_diff_vec <- as.vector(as.matrix(salinity_diff))
site_or <- rep(c(1:100),100)
site_dest <- rep(c(1:100),each=100)
data <- data.frame(LA_vec,sea_dist_vec,salinity_diff_vec,site_or,site_dest)


# Fit 4th order polynomial curves to each population
# As a function of sea distances
modp <- list()
for (i in 1 : 100) {
  data.pop <- data[which(data$site_or==i),]
  modp[[i]] <- lm(LA_vec~poly(sea_dist_vec,4,raw=TRUE),data.pop)
}

# Plot polynomial models
xnew <- data.frame(sea_dist_vec=seq(0,4000,100))
iplot <- 1
for (i in 1 : 100) {
  if (i %in% c(1,31,61,91)) {
    pdf(paste0("plot_seadist_",iplot,".pdf"),8,12)
    par(mar=c(2,2,1,1),mfrow=c(10,3))
    iplot <- iplot +1
  }
  data.pop <- data[which(data$site_or==i),]
  plot(LA_vec~sea_dist_vec,data.pop,pch=16)
  if (is.na(deme$Name[i])) {
    title(paste0(i,", R^2 = ",format(summary(modp[[i]])$adj.r.squared,digits=2)))
  } else {
    title(paste0(i,", ",deme$Name[i],", R^2 = ",format(summary(modp[[i]])$adj.r.squared,digits=2)))
  }
  abline(h=0,lty=2)
  lines(xnew$sea_dist_vec, predict(modp[[i]],xnew), col="red" )
  if (i %in% c(30,60,90,100)) dev.off()
}

r.squared.sea_dist <- vector()
for (i in 1 : 100){
  r.squared.sea_dist[i] <- summary(modp[[i]])$adj.r.squared
}


### As a function of salinity difference
# Fit 4th order polynomial curves to each population
mods <- list()
for (i in 1 : 100) {
  data.pop <- data[which(data$site_or==i),]
  mods[[i]] <- lm(LA_vec~poly(salinity_diff_vec,4,raw=TRUE),data.pop)
}

# Plot polynomial models
xnew <- data.frame(salinity_diff_vec=seq(0,7,0.1))
iplot <- 1
for (i in 1 : 100) {
  if (i %in% c(1,31,61,91)) {
    # png(paste0("plot_salinitydiff_",iplot,".png"),20,30,units="cm",res=600)
    pdf(paste0("plot_salinitydiff_",iplot,".pdf"),8,12)
    par(mar=c(2,2,1,1),mfrow=c(10,3))
    iplot <- iplot +1
  }
  data.pop <- data[which(data$site_or==i),]
  plot(LA_vec~salinity_diff_vec,data.pop,pch=16)
  if (is.na(deme$Name[i])) {
    title(paste0(i,", R^2 = ",format(summary(mods[[i]])$adj.r.squared,digits=2)))
  } else {
    title(paste0(i,", ",deme$Name[i],", R^2 = ",format(summary(mods[[i]])$adj.r.squared,digits=2)))
  }
  abline(h=0,lty=2)
  lines(xnew$salinity_diff_vec, predict(mods[[i]],xnew), col="red" )
  if (i %in% c(30,60,90,100)) dev.off()
}

r.squared.salinity_diff <- vector()
for (i in 1 : 100){
  r.squared.salinity_diff[i] <- summary(mods[[i]])$adj.r.squared
}


# Statistics
summary(r.squared.sea_dist)
sd(r.squared.sea_dist)
summary(r.squared.salinity_diff)
sd(r.squared.salinity_diff)


# Biplot with color gradient
# Figure 4
# Define colors
data$col.level <- cut(data$LA_vec,seq(-0.8,0.8,0.2))
colors.graph <- brewer.pal(8,"RdBu")

data.1 <- data[order(data$LA_vec,decreasing=T),] # To highlight red points

# Plot figure
png("Figure_4.png",width=8,height=5,res=300,units="cm")
par(mar=c(4,4,1,1))
plot(data.1$sea_dist_vec, data.1$salinity_diff_vec,
     xlab = "Spatial distance (sea)", ylab = "Salinity difference",
     col = colors.graph[data.1$col.level], pch = 16, cex = 0.5,las=1,cex.axis=0.75)
dev.off()
# Plot legend
png("Figure_4_legend.png",width=8,height=5,res=300,units="cm")
par(mar=c(0,0,0,0))
plot.new()
legend(0,1,levels(data.1$col.level),pch=16,col=colors.graph,ncol=2,
       title="Local adaptation")
dev.off()

# Find one negative LAij point to show as an example
tail(data.1,20)
deme$Name[44]
deme$Name[20]
plot(data.1$sea_dist_vec, data.1$salinity_diff_vec,
     xlab = "Spatial distance (sea)", ylab = "Salinity difference",
     col = colors.graph[data.1$col.level], pch = 16, cex = 0.5,las=1,cex.axis=0.75)
id <- which(data.1$site_or==44 & data.1$site_dest==20)
points(data.1$sea_dist_vec[id], data.1$salinity_diff_vec[id],col="black",cex=1)

# Correlation between spatial distance and environmental distance (Mantel test)
library(ape)
mantel.test(sea_dist,salinity_diff,graph=T)


