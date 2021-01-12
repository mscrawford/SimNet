##################################################################################
##################################################################################
############################## SimNet simulations ################################
##################################################################################
##################################################################################

##### 24th June 2019




##################################################################################
######################  defining the pool of 64 species ##########################
##################################################################################

setwd("/Users/marechaux/OneDrive/SimNet/TROLL")
par_file="input_v241abs.txt"
numsp=190
species_par=read.table(par_file, header=TRUE,dec=".", sep="", skip=57, row.names=1)
#### n=190 species 
summary(species_par)
nrow(species_par)

library(Hmisc)
corr_all<-rcorr(species_par[,1:7], type="pearson")
plot(species_par[,1:7])
species_par$hmax
length(which(species_par$hmax==50.92))
species_par$ah
length(which(species_par$ah==0.308))

plot(species_par[,1:5])

############# making 64 clusters based on species traits ###########

trait_scale <- scale(species_par[,1:5]) ### scaling the data, so that the results are not biased because of different units.
set.seed(123)
cluster <- kmeans(trait_scale, 64, nstart=1000) #### paritioning clustering with k-means
cluster$cluster
summary(as.factor(cluster$cluster))
summary(as.factor(cluster$size))
species_par[which(cluster$cluster==14),1:7]

############ pick one species per cluster at random ############

species_64=NULL

for (c in 1:64) {
  if (c==1)  {
    species_64 <- c(species_64, which(cluster$cluster==c))}
  else { 
    species_64 <- c(species_64, sample(which(cluster$cluster==c), 1))}
}

species_par_64 <- species_par[species_64,]

summary(species_par_64)

##################################################################################
############ create parameter input files ##################
##################################################################################

setwd("/Users/marechaux/OneDrive/SimNet/TROLL/Simulations")

### n=64 species
write.table(species_par_64, file="par_sp64.txt", append=TRUE, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)

seed_rain=c(0,5,10,50,100,1000)
### n=1 species (64 simulations)

for (N in seed_rain) {
name_par=paste(paste("par", N, sep="_"), "sp1", sep="_")
for (s in 1:64) {
  file.copy(paste(name_par, ".txt", sep=""), paste(paste(name_par, as.character(s), sep="_"), ".txt", sep=""))
  write.table(as.data.frame(species_par_64[s,]), file=paste(paste(name_par, as.character(s), sep="_"), ".txt", sep=""), append=TRUE, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
}
}

### n=2 species (64 simulations)

for (N in seed_rain) {
  name_par=paste(paste("par", N, sep="_"), "sp2", sep="_")
for (i in 1:64) {
  file.copy(paste(name_par, ".txt", sep=""), paste(paste(name_par, as.character(i), sep="_"), ".txt", sep=""))
  write.table(as.data.frame(species_par_64[sample(1:64,2),]), file=paste(paste(name_par, as.character(i), sep="_"), ".txt", sep=""), append=TRUE, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
}
}

### n=4 species (64 simulations)


for (N in seed_rain) {
  name_par=paste(paste("par", N, sep="_"), "sp4", sep="_")
  for (i in 1:64) {
    file.copy(paste(name_par, ".txt", sep=""), paste(paste(name_par, as.character(i), sep="_"), ".txt", sep=""))
    write.table(as.data.frame(species_par_64[sample(1:64,4),]), file=paste(paste(name_par, as.character(i), sep="_"), ".txt", sep=""), append=TRUE, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
  }
}

### n=8 species (64 simulations)

for (N in seed_rain) {
  name_par=paste(paste("par", N, sep="_"), "sp8", sep="_")
  for (i in 1:64) {
    file.copy(paste(name_par, ".txt", sep=""), paste(paste(name_par, as.character(i), sep="_"), ".txt", sep=""))
    write.table(as.data.frame(species_par_64[sample(1:64,8),]), file=paste(paste(name_par, as.character(i), sep="_"), ".txt", sep=""), append=TRUE, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
  }
}

### n=16 species (64 simulations)

for (N in seed_rain) {
  name_par=paste(paste("par", N, sep="_"), "sp16", sep="_")
  for (i in 1:64) {
    file.copy(paste(name_par, ".txt", sep=""), paste(paste(name_par, as.character(i), sep="_"), ".txt", sep=""))
    write.table(as.data.frame(species_par_64[sample(1:64,16),]), file=paste(paste(name_par, as.character(i), sep="_"), ".txt", sep=""), append=TRUE, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
  }
}

### n=32 species (64 simulations)

for (N in seed_rain) {
  name_par=paste(paste("par", N, sep="_"), "sp32", sep="_")
  for (i in 1:64) {
    file.copy(paste(name_par, ".txt", sep=""), paste(paste(name_par, as.character(i), sep="_"), ".txt", sep=""))
    write.table(as.data.frame(species_par_64[sample(1:64,32),]), file=paste(paste(name_par, as.character(i), sep="_"), ".txt", sep=""), append=TRUE, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)
  }
}
