# R code from the manuscript

# Title: "A multi-locus demo-genetic simulator to model local adaptation in heterogeneous landscapes"
# Authors: Marco Andrello, Christelle Noirot, Florence Débarre and Stéphanie Manel
# Submitted to: Molecular Ecology Resources
# Date: 05 Aug 2019

# "01a - Simulations.R"
# Code to parametrize and execute the example simulation on Mullus surmuletus
# Marco Andrello
# 05/08/2019


rm(list=ls())

library(MetaPopGen)


###########################################################
# Defining simulation paramters
###########################################################

# Number of alleles per locus
List_gene<-c(2,2,2,2)


# Mutation rate per locus
mu<-c(1e-06,1e-06,1e-06,1e-06)


# Recombination rate (used only if there are only two loci)
r<-0.5


# Matrix defining the index of each genotype in the MetaPopGen objects
index <- genotype.index.multilocus(List_gene)


# Matrix giving the probability of each gamete type for each parental genotype
Proba <- create.probability.matrix(index,List_gene,r,mu)


# General parameters
n <- 100 # Number of demes
l <- dim(Proba)[1] # Number of multilocus gametes
m <- dim(Proba)[2] # Number of genotypes
z <- 5 # Number of age-classes
T_max <- 100 # Number of years of simulation


# Number of phenotypes
K = 2* (length(List_gene) - 1) + 1 # length(List_gene) - 1 because the fourth locus is neutral


# Initialize survival and fecundity objects
sigma_M <- array(0,dim=c(m,n,z,T_max))
sigma_F <- array(0,dim=c(m,n,z,T_max))
phi_M <- array(0,dim=c(m,n,z,T_max))
phi_F <- array(0,dim=c(m,n,z,T_max))


# Load salinity
# Salinity_per_deme[3,] is salinity in PSU
load("Salinity_data.RData")


# Defining survival as a function of phenotype and salinity using an exponential function
# See eq. 1 in the manuscript
# omega = 1 and xi ranging from 36 to 39 with 0.5 steps
xi <- seq(36,39,0.5) # Optima i.e. phenotype
x <- Salinity_per_deme[3,] # Salinity
omega <- 1
Surv <- array(NA,dim=c(n,K))
for(i in 1 : K) {
  Surv[,i] <- 0.8 * exp( -(xi[i] - x)^2 / (2*omega^2) )
}


# Computing the number of "+" alleles  per multilocus genotype
# See also Supplementary Table 1. Note that Suppl Table 1 contains 36 multilocus genotypes while here we have 136 multilocus genotypes
# because here the double, triple and quadruple heterozygotes are distinct (i.e A1B1/A2B2 is distinct from A1B2/A2B1 and so on)
al<-array(0,dim=m)
for (k in 1:m){
  if(substr(colnames(Proba)[k],2,2)=="1"){
    al[k]<-al[k]+1
  }
  if(substr(colnames(Proba)[k],4,4)=="1"){
    al[k]<-al[k]+1
  }
  if(substr(colnames(Proba)[k],6,6)=="1"){
    al[k]<-al[k]+1
  }
  if(substr(colnames(Proba)[k],11,11)=="1"){
    al[k]<-al[k]+1
  }
  if(substr(colnames(Proba)[k],13,13)=="1"){
    al[k]<-al[k]+1
  }
  if(substr(colnames(Proba)[k],15,15)=="1"){
    al[k]<-al[k]+1
  }
}


# Assigning survival to multilocus genotypes as a function of the number of "+" alleles
for (x in 1 : z) {
  for (t in 1 : T_max) {
    for(k in 1:m) {
      sigma_M[k,,x,t] <- Surv[,al[k]+1]
      sigma_F[k,,x,t] <- Surv[,al[k]+1]
    }
  }
}


# Assigning fecundities
# Age-class 1 has zero fecundity
for (geno in 1 : m) {
  for (t in 1 : T_max) {
    for(deme in 1 : n) {
      phi_M[geno,deme,2:5,t] <- 1000000 # Male fecundity
      phi_F[geno,deme,2:5,t] <- c(6852,11089,16052,20974) # Female fecundity
    }
  }
}


# Carrying capacity of each deme
kappa0<-array(5000,dim=c(n,T_max))


# Loading larval dispersal matrix
load("Dispersal_SC.RData")
delta <- array(Dispersal_SC,dim=c(n,n,T_max)) 


# Save the simulation paramters 
save(m,n,z,T_max,
     sigma_M,sigma_F,phi_M,phi_F,
     List_gene,
     mu,r,Proba,
     delta,kappa0,
     al,
     file="Simulation_parameters.RData")



###################################################################################
# Run the simulations
###################################################################################

rm(list=ls())

library(MetaPopGen)

# Load dataset (created in Part 1)
load("Simulation_parameters.RData")

# Simulation
nrepl <- 10

# Run the simulations
for (rr in 1 : nrepl) {
  
  #  Initial population, at random
  print("Initializing populations")
  N1_M <- array(NA,dim=c(m,n,z)) # Initial number of male individuals
  N1_F <- array(NA,dim=c(m,n,z)) # Initial number of female individuals
  for (deme in 1 : n){
    for (age in 1 : z) {
      for (genotype in 1 : m) {
        N1_M[genotype,deme,age] <- round(runif(1)*150)
        N1_F[genotype,deme,age] <- round(runif(1)*150)
      }
    }
  }
  
  
  cat("Replicate",rr,"\n") # Printing on screen the progress of the simulation
  
  sim.metapopgen.dioecious.multilocus(
    N1_M = N1_M, N1_F = N1_F,
    sigma_M = sigma_M, sigma_F = sigma_F,
    phi_M = phi_M, phi_F = phi_F,
    List_gene=List_gene, mu=mu,r=r,
    delta = delta, recr.dd = "adults", kappa0 = kappa0,
    T_max = T_max,
    save.res = T, verbose = F)
}




