# R code from the manuscript

# Title: "A multi-locus demo-genetic simulator to model local adaptation in heterogeneous landscapes"
# Authors: Marco Andrello, Christelle Noirot, Florence Débarre and Stéphanie Manel
# Submitted to: Molecular Ecology Resources
# Date: 05 Aug 2019

# "01b - Reading_results.R"
# Code to read the results of the demo-genetic simulations on Mullus surmuletus
# Marco Andrello
# 05/08/2019

rm(list=ls())

library(MetaPopGen)


# Load the simulation parameters
# load("Analysis_Marco.RData")
load("Simulation_parameters.RData")


#definition de la direction des outputs des simulations
# Read the names of the folders containing the simulation results
list.dir <- list.dirs()
list.dir <- list.dir[-1]
name<-getwd() 


# Matrix defining the index of each genotype in the MetaPopGen objects
index<-genotype.index.multilocus(List_gene)


# Number of replicates
nrepl <- length(list.dir)


# Create the array for fst
Globfst<-array(NA,dim=c(nrepl,T_max))
fst_A <- array(NA,dim=c(nrepl,T_max))
fst_B <- array(NA,dim=c(nrepl,T_max))
fst_C <- array(NA,dim=c(nrepl,T_max))
fst_D <- array(NA,dim=c(nrepl,T_max))


# Nt groups the results of all replicates at each time step disregarding sex and age-class
Nt<-array(NA,dim=c(m,n,T_max,nrepl))
Nt_A<-array(NA,dim=c(4,n,T_max,nrepl))
Nt_B<-array(NA,dim=c(4,n,T_max,nrepl))
Nt_C<-array(NA,dim=c(4,n,T_max,nrepl))
Nt_D<-array(NA,dim=c(4,n,T_max,nrepl))


# Calculating Global FST and locus-specific Fst for each replicate at each time-step
for (rr in 1 : nrepl) {
  
  dir.name <- list.dir[rr]
  
  # Show the progress of the calculation
  print(dir.name)
  flush.console()
  
  for (t in 1 : T_max){

    load(paste0(name,"/",dir.name,"/N",t,".RData"))
    
    N<-N_F+N_M
    N_z<-apply(N,c(1,2),sum)
    dimnames(N)<-list(colnames(Proba),deme=c(1:n))
    Nt[,,t,rr]<-N_z
    
    mA<-1+List_gene[1]*(List_gene[1]+1)/2
    mB<-1+List_gene[2]*(List_gene[2]+1)/2
    mC<-1+List_gene[3]*(List_gene[3]+1)/2
    mD<-1+List_gene[4]*(List_gene[4]+1)/2
    
    N_A<-array(NA,dim=c(mA,n))
    N_B<-array(NA,dim=c(mB,n))
    N_C<-array(NA,dim=c(mC,n))
    N_D<-array(NA,dim=c(mD,n))
    
    name_A<-c("A1A1","A1A2","A2A1","A2A2")
    name_B<-c("B1B1","B1B2","B2B1","B2B2")
    name_C<-c("C1C1","C1C2","C2C1","C2C2")
    name_D<-c("D1D1","D1D2","D2D1","D2D2")
    
    dimnames(N_A)<-list(name_A,c(1:n))
    dimnames(N_B)<-list(name_B,c(1:n))
    dimnames(N_C)<-list(name_C,c(1:n))
    dimnames(N_D)<-list(name_D,c(1:n))
    
    for (k in 1:mA){
      A<-which(rownames(N_A)[k]== paste0(substr(colnames(Proba),1,2),substr(colnames(Proba),10,11)))
      N_A[k,]<-apply(N_z[A,],2,sum)
    }
    for (k in 1:mB){
      B<-which(rownames(N_B)[k]== paste0(substr(colnames(Proba),3,4),substr(colnames(Proba),12,13)))
      N_B[k,]<-apply(N_z[B,],2,sum)
    }
    for (k in 1:mC){
      C<-which(rownames(N_C)[k]== paste0(substr(colnames(Proba),5,6),substr(colnames(Proba),14,15)))
      N_C[k,]<-apply(N_z[C,],2,sum)
    }
    for (k in 1:mC){
      D<-which(rownames(N_D)[k]== paste0(substr(colnames(Proba),7,8),substr(colnames(Proba),16,17)))
      N_D[k,]<-apply(N_z[D,],2,sum)
    }
    Nt_A[,,t,rr]<-N_A
    Nt_B[,,t,rr]<-N_B
    Nt_C[,,t,rr]<-N_C
    Nt_D[,,t,rr]<-N_D
    fst_A[rr,t]<- fst(N_A)
    fst_B[rr,t]<-fst(N_B)
    fst_C[rr,t]<- fst(N_C)
    fst_D[rr,t]<-fst(N_D)
    Globfst[rr,t]<-fst(Nt[,,t,rr])
    
  }
}


# FST at the end of the simulation
# Mean and standard deviation over replicates
mean(fst_A[,T_max])
mean(fst_B[,T_max])
mean(fst_C[,T_max])
mean(fst_D[,T_max])
sd(fst_A[,T_max])
sd(fst_B[,T_max])
sd(fst_C[,T_max])
sd(fst_D[,T_max])


# Calculating average survival in each deme for individuals originating each deme
# This is needed to calculate local adaptation as per eq. 3 in the manuscript
surv_trans <- array(NA,dim=c(n,n,nrepl))
t <- T_max
for(rr in 1: nrepl){
  for (i in 1:n){
    sigmaik <- sigma_F[,i,1,1]
    for (j in 1:n){
      vjk <- Nt[,j,t,rr] / sum(Nt[,j,t,rr])
      surv_trans[i,j,rr] <- sum(sigma_F[,i,1,1] * vjk)
    }
  }
}
surv.trans <- apply(surv_trans,c(1,2),mean)


# AL: number of individual of a given phenotype in each deme
AL<-array(0,dim=c(100,7))
for (i in 1:n){
  for(k in 1:m){
    AL[i,al[k]+1]<-AL[i,al[k]+1]+N_z[k,i]
  }
}


# Calculating reproductive output from the age structure and total fecundity in each deme
# This is needed to calculate larval flow (Figure 3 in the manuscript)
# (Using only one replicate because replicates are very similar)
name.dir <- "/2019-Aug-05-10.07.27" # Replace this with one of the folder created by the simulation 
fectot <- array(NA,c(n,50)) # Inspecting the last 50 years, when the population genetics are stable
for (i in 1 : 50) {
  name.file <- paste0(getwd(),name.dir,"/N",(i+50),".RData")
  load(name.file)
  nftot <- apply(N_F,c(2,3),sum)
  fec <- phi_F[1,1,,1]
  fectot[,i] <- as.matrix(nftot) %*% as.matrix(fec)
}
ro <- rowMeans(fectot)

save(AL,
     Globfst,
     Nt, Nt_A, Nt_B, Nt_C, Nt_D,
     fst_A, fst_B, fst_C, fst_D,
     surv.trans,
     ro,
     file="Results_Base.RData" )


