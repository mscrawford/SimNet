#Run model given no3i, pNi, and abmi for each species
getbmest<-function(no3i, pNi, abmi) {
  abm_esti<-numeric(length(abmi))
  
  subsi<-no3i==min(no3i, na.rm=T)
  abm_esti[subsi]<-abmi[subsi]
  
  pNuptake<-pNi
  
  ord<-order(no3i)
  ord<-ord[!(ord%in%which(subsi))]
  
  if(sum(ord)>0) {
    for(j in 1:length(ord)) {
      abm_esti[ord[j]]<-abmi[ord[j]]-sum((pNuptake*abm_esti)/pNi[ord[j]])
      abm_esti[!is.finite(abm_esti)]<-0
      if(abm_esti[ord[j]]<0) {
        abm_esti[ord[j]]<-0
      }
    }
  }
  return(abm_esti)
}

#logit functions
logit<-function(x) {
  suppressWarnings(res<-(-log(1/x-1)))
  res[!is.finite(res)]<-NA
  res
}

ilogit<-function(x) {
  1/(1+exp(-x))
}

#transform traits to PCA space
traits_to_pca<-function(dat, dimuse=2) {
  #transformed data
  trdat<-data.frame(lno3i=log(dat$no3i), labmi=log(dat$abmi), lpNi=logit(dat$pNi))
  #standardized transformed data
  stddat<-apply(trdat, 2, function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))
  
  #run PCA
  prout<-prcomp(trdat, center=TRUE, scale=TRUE)
  varexp<-round(prout$sdev^2/sum(prout$sdev^2),2)
  psuse<-rev(order(varexp))[1:dimuse]
  
  #project data to pca, retaining dimuse dimensions
  pcaout<-as.matrix(stddat)%*%as.matrix(prout$rotation[,psuse])
  
  return(list(pcaout=pcaout, trdat=trdat, stddat=stddat, prout=prout, varexp=varexp))
}

#back-transform PCA variables to traits
pca_to_traits<-function(pcaout=pcaout, prout=prout, trdat=trdat, dimuse=2) {
  varexp<-round(prout$sdev^2/sum(prout$sdev^2),2)
  psuse<-rev(order(varexp))[1:dimuse,drop=FALSE]
  
  stddat<-as.matrix(pcaout)%*%t(as.matrix(prout$rotation[,psuse,drop=FALSE]))
  trdat_new<-data.frame(sapply(1:ncol(stddat), function(i) {stddat[,i,drop=FALSE]*sd(trdat[,i],na.rm=T)+mean(trdat[,i],na.rm=T)}))
  if(dim(trdat_new)[1]!=dim(stddat)[1]) {
    trdat_new<-data.frame(t(trdat_new)); row.names(trdat_new)<-NULL
  }
  
  colnames(trdat_new)<-colnames(trdat)
  dat_new<-data.frame(no3i=exp(trdat_new$lno3i), abmi=exp(trdat_new$labmi), pNi=ilogit(trdat_new$lpNi))
  
  return(dat_new)
}



#examples:
if(FALSE) {
  #load data
  dat<-read.table("~/Dropbox/Projects/047_Barry_BEF/SIMNET_meeting/src/resource_tradeoff_model/tableS1_paramtertable.csv", header=T, sep=";")[,c("Updated.Species.Name", "Measured.Soil.Nitrate.mg.kg.1", "Measured.Aboveground.Biomass.g.m.2", "Measured.Proportion.Tissue.N")]
  colnames(dat)<-c("species", "no3i", "abmi", "pNi")
  modout<-getbmest(no3i = dat$no3i, pNi=dat$pNi, abmi = dat$abmi)
  
  pca_result<-traits_to_pca(dat)
  traits_out<-pca_to_traits(pcaout=pca_result$pcaout, prout=pca_result$prout, trdat=pca_result$trdat)
  
  plot(traits_out$no3i, dat$no3i, xlab="tradeoff", ylab="raw", main="no3"); abline(a=0, b=1, lty=3)
  plot(traits_out$abmi, dat$abmi, xlab="tradeoff", ylab="raw", main="abm"); abline(a=0, b=1, lty=3)
  plot(traits_out$pNi, dat$pN, xlab="tradeoff", ylab="raw", main="pN"); abline(a=0, b=1, lty=3)
  
  modout_tradeoff<-getbmest(no3i = traits_out$no3i, abmi = traits_out$abmi, pNi = traits_out$pNi)
  modout_raw<-getbmest(no3i = dat$no3i, abmi = dat$abmi, pNi =dat$pNi)
  
  plot(modout_tradeoff, modout_raw, xlab="tradeoff", ylab="raw", main="predicted abundance"); abline(a=0, b=1, h=0, lty=3)
}
