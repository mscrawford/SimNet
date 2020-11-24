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

trmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    ord <- order(no3)
    dN <- numeric(length(State))
    for(j in 1:length(State)) {
      if(State[ord[j]]>0) {
        dN[ord[j]]<-r*pmax(0,State[ord[j]])*(abm[ord[j]]-
                         sum((pN[ord[1:j]]*pmax(0,State[ord[1:j]]))/pN[ord[j]]))/abm[ord[j]]
      } else {
        dN[ord[j]]<-0
      }
    }
    
    return(list(c(dN)))
  })
}

#load data
dat<-read.table("tableS1_paramtertable.csv", header=T, sep=";")[,c("Updated.Species.Name", "Measured.Soil.Nitrate.mg.kg.1", "Measured.Aboveground.Biomass.g.m.2", "Measured.Proportion.Tissue.N")]
colnames(dat)<-c("species", "no3i", "abmi", "pNi")

#Make pool of 100 species based on real species traits
trdat<-data.frame(lno3=log(dat$no3i),labmi=log(dat$abmi),lpNi=logit(dat$pNi))
cv<-cov(trdat)
mn<-colMeans(trdat)

set.seed(230519)
trsmp<-rmvnorm(n = 64, mean = mn, sigma = cv)
dat2<-data.frame(no3i=exp(trsmp[,1]), abmi=exp(trsmp[,2]), pNi=ilogit(trsmp[,3]))

pca_result<-traits_to_pca(dat2)
traits_out0<-pca_to_traits(pcaout=pca_result$pcaout, prout=pca_result$prout, trdat=pca_result$trdat)

#matrix for saving results
datout<-NULL

#mean recruits
seed_input<-mean(traits_out0$abmi)*0.14

write.csv(data.frame(SpeciesID=1:64, traits_out0), "adammod_trans_exp_speciesdata.csv", row.names=F)

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
  
  #if(n==1) { #if monoculture, then run each species, one at a time
  #  ps<-simmat$rep[[ii]]
  #} else { #if not a monoculture, then choose random species
     ps<-sppslst[[smpps]]
  #}
  
  traits_out<-traits_out0[ps,]
  
  pars  <- list(no3=traits_out$no3i,
                abm=traits_out$abmi,
                pN=traits_out$pNi,
                r=0.2)
  
  yini  <- rep(0.1,length(pars$no3))
  times <- seq(0, 1, by = 0.01)
  
  #matrix for saving temporary simulation results
  datout_FULL<-matrix(nrow=100, ncol=n)
  
  #define seed rain magnitude
  if(simmat$type[ii]=="0") {
    db<-seed_input*0/n
  } else if(simmat$type[ii]=="5") {
    db<-seed_input*0.05/n
  } else if(simmat$type[ii]=="10") {
    db<-seed_input*0.1/n
  } else if(simmat$type[ii]=="50") {
    db<-seed_input*0.5/n
  } else if(simmat$type[ii]=="100") {
    db<-seed_input*1/n
  } else if(simmat$type[ii]=="1000") {
    db<-seed_input*10/n
  }
  
  #run experiment while seed rain in coming in
  datout_FULL[1,]<-yini
  y0<-yini
  for(im in 2:100) {
    if(im!=2) {
      y0<-y0+db
    }
    out   <- ode(y = y0, times = times, func = trmod, parms = pars, method = "ode45")
    out[out<0]<-0
    y0<-out[nrow(out),-1]
    datout_FULL[im,]<-y0
    #matplot(out[,1], out[,-1], type="l")
  }
  
  #repeat after turning off rain
  y0<-datout_FULL[nrow(datout_FULL),]
  times2<-seq(0, 99, by=1)
  out   <- ode(y = y0, times = times2, func = trmod, parms = pars, method = "ode45")
  
  if(FALSE) {
    #plotting
    par(mfrow=c(1,2), mar=c(4,4,2,2))
    matplot(1:nrow(datout_FULL), datout_FULL, type="l", xlab="time", ylab="biomass"); abline(h=0, lty=3)
    matplot(out[,1], out[,-1], type="l", xlab="time", ylab="biomass"); abline(h=0, lty=3)
  }
    
  #save results
  datout<-rbind(datout,
                data.frame(Model="adam",
                           Ninitial=n,
                           Rep=simmat$rep[[ii]],
                           SeedRain=as.character(simmat$type[ii]),
                           SpeciesID=ps,
                           Stage = rep(c("assembly", "disassembly"), each=n*nrow(out)),
                           Year = rep(1:(nrow(out)*2), each=n),
                           Biomass = NA,
                           Productivity = c(t(datout_FULL), t(out[,-1]))))
  
  if(ii/10==floor(ii/10)) {
    print(round(ii/nrow(simmat),3))
    trtable<-save(list = c("datout"), file = "adammod_trans_exp_out_HPC.rda")
  }
}
trtable<-save(list = c("datout"), file = "adammod_trans_exp_out_HPC.rda")
