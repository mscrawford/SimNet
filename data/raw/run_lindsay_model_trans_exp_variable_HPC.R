#!/usr/bin/env Rscript
#catcherror
rm(list=ls())
#setwd("~/Dropbox/Projects/047_Barry_BEF/SIMNET_meeting/src/transience_exp/")
#source("make_trait_grid.R")
source("model_functions.R")
require(mvtnorm)
require(deSolve)
######

#seed rain options are "0", "5", "10", "50", "100", or "1000": i.e. no rain, 5%, 10%, 50%, 100%, or 1000%)
niter<-64
simmat<-data.frame(ntot=rep(c(1,2,4,8,16,32,64),each=(6*niter)),
                   type=rep(c("0", "5", "10", "50", "100", "1000"),(7*niter)),
                   rep=rep(1:(7*niter), each=6))
samplemat<-unique(simmat[,c("ntot", "rep")])

#model parameters
rho   = 10                  # CN ratio
N0    = 0.07                # N concentration in soil at beginning of season
ind0  = 1000                # total number of individuals
tmax = 140                  # timesteps per year
yrmax = 100                 # number of years

set.seed(230519)
#Make pool of 100 species based on real species traits
#traits_out0<-data.frame(Vi=seq(1000,10000,length=64))
traits_out0<-data.frame(Vi=rnorm(64, (10000-1000)/2, (10000-1000)/4))
traits_out0$thetai<-(18000/(traits_out0$Vi+3000))
traits_out0$Vi = traits_out0$Vi + rnorm(nrow(traits_out0), 0, sd = mean(traits_out0$Vi)/5)
plot(traits_out0$Vi, traits_out0$thetai)
traits_out0$Vi[traits_out0$Vi<0] = 0

#matrix for saving results
datout<-NULL

write.csv(data.frame(SpeciesID=1:64, traits_out0), "lindsaymod_trans_exp_speciesdata_variable.csv", row.names=F)

#make species list
sppslst<-list()
for(i in 1:nrow(samplemat)) {
  if(samplemat$ntot[i]==1) {
    sppslst[[i]]<-samplemat$rep[i]
  } else {
    sppslst[[i]]<-sample(64, samplemat[i,"ntot"])
  }
}

for(ii in 1:nrow(simmat)) {
  n=simmat$ntot[ii] #number of species
  smpps<-which(samplemat$ntot==simmat$ntot[ii] & samplemat$rep==simmat$rep[ii])
  indi0  = rep(ind0/n,n)      # individuals per species
  M0i = rep(2/n, n)/indi0     # Starting biomass per individual - one per species  
  
  #if(n==1) { #if monoculture, then run each species, one at a time
  #  ps<-simmat$rep[[ii]]
  #} else { #if not a monoculture, then choose random species
     ps<-sppslst[[smpps]]
  #}
  
  traits_out<-traits_out0[ps,]
  
  #run model
  Mout<-array(dim=c(n, tmax, yrmax))
  Mout[,1,1]<-M0i
  
  Iout<-array(dim=c(n, yrmax))
  Iout[,1]<-indi0
  
  Mout_TOT<-array(dim=c(n, tmax, yrmax))
  Mout_TOT[,1,1]<-M0i*indi0
  
  indi<-indi0
  
  if(as.character(simmat$type[[ii]])=="0") {
    fracrep<-0
  } else if(as.character(simmat$type[[ii]])=="5") {
    fracrep<-0.05
  } else if(as.character(simmat$type[[ii]])=="10") {
    fracrep<-0.1
  } else if(as.character(simmat$type[[ii]])=="50") {
    fracrep<-0.5
  } else if(as.character(simmat$type[[ii]])=="100") {
    fracrep<-1
  } else if(as.character(simmat$type[[ii]])=="1000") {
    fracrep<-10
  }
  
  for(i in 1:yrmax) {
    Iout[,i]<-indi
    Mout[,1,i]<-M0i
    
    for(j in 2:tmax) {
      Itmp<-traits_out$thetai*Mout[,(j-1),i]*N0*
        ((N0*traits_out$Vi*rho-sum(indi*Mout[,(j-1),i]))/(N0*traits_out$Vi*rho))
      Itmp[Itmp<0]<-0
      
      Mout[,j,i]=Mout[,(j-1),i]+Itmp
    }
    Mout_TOT[,,i]<-Mout[,,i]*indi
    
    #redistribute individuals into new seeds
    indi<-round(((Mout_TOT[,j,i]/sum(Mout_TOT[,j,i]))+(1/n)*fracrep)*ind0)
  }
  
  tmp<-Mout_TOT[,140,,drop="FALSE"] #rows are species, cols are years
  
  
  #repeat after turning off rain
  if(FALSE) {
    #Potentially wrong? Need to start with end of simulation's biomass
    
    Mout<-array(dim=c(n, tmax, yrmax))
    Mout[,1,1]<-M0i
    
    Iout<-array(dim=c(n, yrmax))
    Iout[,1]<-indi0
    
    Mout_TOT<-array(dim=c(n, tmax, yrmax))
    Mout_TOT[,1,1]<-M0i*indi0
    
    indi<-indi0
  } else {
    Mout<-array(dim=c(n, tmax, yrmax))
    Mout[,1,1]<-M0i
    
    Iout<-array(dim=c(n, yrmax))
    Iout[,1]<-indi
    
    Mout_TOT<-array(dim=c(n, tmax, yrmax))
    Mout_TOT[,1,1]<-M0i*indi
    
    indi<-indi
  }
  
  for(i in 1:yrmax) {
    Iout[,i]<-indi
    Mout[,1,i]<-M0i
    
    for(j in 2:tmax) {
      Itmp<-traits_out$thetai*Mout[,(j-1),i]*N0*
        ((N0*traits_out$Vi*rho-sum(indi*Mout[,(j-1),i]))/(N0*traits_out$Vi*rho))
      Itmp[Itmp<0]<-0
      
      Mout[,j,i]=Mout[,(j-1),i]+Itmp
    }
    Mout_TOT[,,i]<-Mout[,,i]*indi
    
    #redistribute individuals into new seeds
    indi<-round(Mout_TOT[,j,i]/sum(Mout_TOT[,j,i])*ind0)
  }
  
  tmp2<-Mout_TOT[,140,,drop="FALSE"] #rows are species, cols are years
  
  
  if(FALSE) {
    #plotting
    par(mfrow=c(1,2), mar=c(4,4,2,2))
    matplot(1:nrow(t((tmp[,1,]))), t((tmp[,1,])), type="l", xlab="time", ylab="biomass"); abline(h=0, lty=3)
    matplot(1:nrow(t((tmp2[,1,]))), t((tmp2[,1,])), type="l", xlab="time", ylab="biomass"); abline(h=0, lty=3)
  }
    
  #save results
  datout<-rbind(datout,
                data.frame(Model="lindsay",
                           Ninitial=n,
                           Rep=simmat$rep[[ii]],
                           SeedRain=as.character(simmat$type[ii]),
                           SpeciesID=ps,
                           Stage = rep(c("assembly", "disassembly"), each=n*(dim(tmp)[3])),
                           Year = rep(1:((dim(tmp)[3])*2), each=n),
                           Biomass = NA,
                           Productivity = c(c(tmp[,1,]), c(tmp2[,1,]))))
  
  if(ii/10==floor(ii/10)) {
    print(round(ii/nrow(simmat),3))
    trtable<-save(list = c("datout"), file = "lindsaymod_trans_exp_out_HPC_variable.rda")
  }
}
trtable<-save(list = c("datout"), file = "lindsaymod_trans_exp_out_HPC_variable.rda")
