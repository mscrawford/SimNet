###################################################################################################
################## TROLL SimNet -- simulation set-up and outputs handling  #########################
###################################################################################################

### IM june 2019

setwd("/Users/marechaux/OneDrive/SimNet/TROLL/Simulations")


####### First step: define how long each simulation has to be run in each phase.
## to do so, run a simulation for the richest mixture (n=64 species) with intense seed rain.

Table1=read.table(paste(name, "_Table1.txt", sep=""), header=F,dec=".", sep="")
assign(paste(name,  "_Table1", sep=""), Table1, pos=.GlobalEnv)
colnames(test64_5_Table1) <- c("Model", "Ninitial", "Rep", "SeedRain", "SpeciesID", "Stage", "Year", "Biomass", "Productivity")
summary(test64_5_Table1)

plot(1:200, test64_5_Table1$Biomass[which(test64_5_Table1$SpeciesID==levels(test64_5_Table1$SpeciesID)[1])], ylim=c(0,40), type="l", xlab="Time", ylab="Biomass")
for (s in 2:length(levels(test64_5_Table1$SpeciesID))) {
  par(new=TRUE)
  plot(1:200, test64_5_Table1$Biomass[which(test64_5_Table1$SpeciesID==levels(test64_5_Table1$SpeciesID)[s])], ylim=c(0,40), type="l", xlab="", ylab="")
}

plot(1:200, test64_5_Table1$Productivity[which(test64_5_Table1$SpeciesID==levels(test64_5_Table1$SpeciesID)[1])], ylim=c(0,2.3), type="l", xlab="Time", ylab="Productivity")
for (s in 2:length(levels(test64_5_Table1$SpeciesID))) {
  par(new=TRUE)
  plot(1:200, test64_5_Table1$Productivity[which(test64_5_Table1$SpeciesID==levels(test64_5_Table1$SpeciesID)[s])], ylim=c(0,2.3), type="l", xlab="", ylab="")
}

## to be conservative ==> 500 years for each phase (assembly and disassembly), that is to say 6000 iterations each, and 12000 in total.
## with outputs every 12000/200=60 timesteps (i.e. every 5 years).

####### Second step: define the level of seed rain that corresponds to 100%
### to do so, run all monoculture simulations with no seed rain and take the average of seed production at final state.

name="sp1_1" #this is the simulation name
par_file="par_sp1_1.txt"   #this is the name of the input file parameter
Table1=read.table(paste(name, "_Table1.txt", sep=""), header=F,dec=".", sep="")
assign(paste(name,  "_Table1", sep=""), Table1, pos=.GlobalEnv)
colnames(sp1_1_Table1) <- c("Model", "Ninitial", "Rep", "SeedRain", "SpeciesID", "Stage", "Year", "Biomass", "Productivity", "Nseed")
summary(sp1_1_Table1)
plot(1:200, sp1_1_Table1$Nseed)
plot(1:200, sp1_1_Table1$Biomass)
plot(1:200, sp1_1_Table1$Productivity)

Nseed=NULL ## Nseed contains the number of seeds produces per one hectare of monoculture and per timestep, average across the 50 last outputs.
Biomass=NULL
Productivity=NULL
for (i in 1:64) {
  name=paste("sp1", i, sep="_")
  Table1=read.table(paste(name, "_Table1.txt", sep=""), header=F,dec=".", sep="")
  colnames(Table1) <- c("Model", "Ninitial", "Rep", "SeedRain", "SpeciesID", "Stage", "Year", "Biomass", "Productivity", "Nseed")
  assign(paste(name,  "_Table1", sep=""), Table1, pos=.GlobalEnv)
  #plot(1:200, Table1$Nseed)
  Nseed <- c(Nseed, mean(Table1$Nseed[150:200]))
  Biomass <- c(Biomass, mean(Table1$Biomass[150:200]))
  Productivity <- c(Productivity, mean(Table1$Productivity[150:200]))
}

mean(Nseed)
# = 3208.868 ~ 3210
summary(Productivity)
summary(Biomass)

####### Create Table 1, by merging all simulations' file.

Initial_species_richness=c(1,2,4,8,16,32,64)
Seed_rain=c(0,5,10,50, 100,1000)

TROLL_Table1=NULL
for (r in Initial_species_richness) {
  for (s in Seed_rain) {
    if (r==64) { 
      i=1
      name <- paste(paste(paste(paste("sp", r, sep=""), s, sep="_"), i, sep="_"), "Table1.txt", sep="_")
      file <- read.table(name, header=F,dec=".", sep="")
      TROLL_Table1 <- rbind(TROLL_Table1, file)
    }
    else {
      for (i in 1:64) {
        name <- paste(paste(paste(paste("sp", r, sep=""), s, sep="_"), i, sep="_"), "Table1.txt", sep="_")
        file <- read.table(name, header=F,dec=".", sep="")
        TROLL_Table1 <- rbind(TROLL_Table1, file)
      }
    }
  }
}

colnames(TROLL_Table1) <- c("Model", "Ninitial", "Rep", "SeedRain", "SpeciesID", "Stage", "Year", "Biomass", "Productivity")

write.table(TROLL_Table1, file = "TROLL_Table1.txt", quote = FALSE, dec=".", sep = "\t", na = "-", row.names=T)


test=read.table("TROLL_Table1.txt", header=T,dec=".", sep="\t")
summary(test)
TROLL_Table1 <-test

TROLL_Table2=read.table("TROLL_Table2.txt", header=T,dec=".", sep="\t")
summary(TROLL_Table2)
nrow(TROLL_Table2)
TROLL_Table2 <- TROLL_Table2[,-c(9,10)]
write.table(TROLL_Table2, file = "TROLL_Table2.txt", quote = FALSE, dec=".", sep = "\t", na = "-", row.names=F)




#dev.new()
s=0
I_prod=NULL
S_prod=NULL
E_prod=NULL
I_mass=NULL
S_mass=NULL
E_mass=NULL
I_rich=NULL
S_rich=NULL
E_rich=NULL
par(mfrow=c(1,2))
for (r in Initial_species_richness) {
  if (r==64) {
      sim=TROLL_Table1[which(TROLL_Table1$SeedRain==s & TROLL_Table1$Ninitial==r & TROLL_Table1$Rep==1),]
      #sr=NULL
      #prod=NULL
      #mass=NULL
      #for (t in 1:200) {
        #sr <- c(sr, length(which(sim$Year==t & sim$Productivity>0)))
        #prod <- c(prod, sum(sim$Productivity[which(sim$Year==t)]))
        #mass <- c(mass, sum(sim$Biomass[which(sim$Year==t)]))
        #if (t==1) { 
      t=1
          I_prod <- c(I_prod,sum(sim$Productivity[which(sim$Year==t)]))
          I_mass <- c(I_mass,sum(sim$Biomass[which(sim$Year==t)]))
          I_rich <- c(I_rich,length(which(sim$Year==t & sim$Productivity>0)))
          
        #}
        #if (t==100) { 
          t=100
          S_prod <- c(S_prod,sum(sim$Productivity[which(sim$Year==t)]))
          S_mass <- c(S_mass,sum(sim$Biomass[which(sim$Year==t)]))
          S_rich <- c(S_rich,length(which(sim$Year==t & sim$Productivity>0)))
          
        #}
       # if (t==200) { 
          t=200
          E_prod <- c(E_prod,sum(sim$Productivity[which(sim$Year==t)]))
          E_mass <- c(E_mass,sum(sim$Biomass[which(sim$Year==t)]))
          E_rich <- c(E_rich,length(which(sim$Year==t & sim$Productivity>0)))
          
        #}
      #}

        #par(new=TRUE, mfg=c(1,1))
       # plot(sr, prod, type="p", pch=20, cex=0.5, col=palette, xlab="", ylab="", main="", yaxt="n", xaxt="n", ylim=c(0,30), xlim=c(0,65))
        #par(new=TRUE, mfg=c(1,2))
        #plot(sr, mass, type="p", pch=20, cex=0.5, col=palette, xlab="", ylab="", main="", yaxt="n", xaxt="n", ylim=c(0,600), xlim=c(0,65))
  }
  else {
  for (i in 1:64) {
    sim=TROLL_Table1[which(TROLL_Table1$SeedRain==s & TROLL_Table1$Ninitial==r & TROLL_Table1$Rep==i),]
    #sr=NULL
    #prod=NULL
    #mass=NULL
    #for (t in 1:200) {
      #sr <- c(sr, length(which(sim$Year==t & sim$Productivity>0)))
      #prod <- c(prod, sum(sim$Productivity[which(sim$Year==t)]))
      #mass <- c(mass, sum(sim$Biomass[which(sim$Year==t)]))
      #if (t==1) { 
    t=1
        I_prod <- c(I_prod,sum(sim$Productivity[which(sim$Year==t)]))
        I_mass <- c(I_mass,sum(sim$Biomass[which(sim$Year==t)]))
        I_rich <- c(I_rich,length(which(sim$Year==t & sim$Productivity>0)))
        
      #}
      #if (t==100) { 
        t=100
        S_prod <- c(S_prod,sum(sim$Productivity[which(sim$Year==t)]))
        S_mass <- c(S_mass,sum(sim$Biomass[which(sim$Year==t)]))
        S_rich <- c(S_rich,length(which(sim$Year==t & sim$Productivity>0)))
        
      #}
      #if (t==200) { 
        t=200
        E_prod <- c(E_prod,sum(sim$Productivity[which(sim$Year==t)]))
        E_mass <- c(E_mass,sum(sim$Biomass[which(sim$Year==t)]))
        E_rich <- c(E_rich,length(which(sim$Year==t & sim$Productivity>0)))
        
      #}
    }
    #if (i==1) {
      #par( mfg=c(1,1))
      #plot(sr, prod, type="p", pch=20, cex=0.5, col=palette, xlab="Species richness", ylab="Productivity",  ylim=c(0,30), xlim=c(0,65), main="")
      #par(mfg=c(1,2))
      #plot(sr, mass, type="p", pch=20, cex=0.5, col=palette, xlab="Species richness", ylab="Biomass", ylim=c(0,800), xlim=c(0,65), main="")
    #}
    #else {
      #par(new=TRUE, mfg=c(1,1))
      #plot(sr, prod, type="p", pch=20, cex=0.5, col=palette, xlab="", ylab="", main="", yaxt="n", xaxt="n", ylim=c(0,30), xlim=c(0,65))
      #par(new=TRUE, mfg=c(1,2))
      #plot(sr, mass, type="p", pch=20, cex=0.5, col=palette, xlab="", ylab="", main="", yaxt="n", xaxt="n", ylim=c(0,800), xlim=c(0,65))
    #}
  }
  }


