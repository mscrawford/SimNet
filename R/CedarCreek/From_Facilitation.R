library(tidyverse)
library(party)

set.seed(1987)

BigBio<-read.csv("CedarCreek-BigBio.csv", sep=";")
#print(dim(BigBio))
#print(str(BigBio))
#print(summary(BigBio))
species.names<-read.csv("BigBio-SpeciesNames.csv", sep=";")
species.traits<-read.csv("CedarCreek_traits.csv", sep=",")
#print(species.traits)
species.traits<-species.traits %>%
	select(-Year)
BigBio<-BigBio[BigBio$NumSp>0,]
BigBio$Biomass..g.m2.<-as.numeric(BigBio$Biomass..g.m2.)
BigBioSpMatrix<-BigBio[,c(4,18:35)]
BigBioSpMatrix<-BigBioSpMatrix[!duplicated(BigBioSpMatrix$Plot),]
BigBio<-BigBio %>%
  group_by(Year, Plot, Species, NumSp)%>%
  summarise(species_biomass_m2=mean(Biomass..g.m2.))%>%
  ungroup()

BigBio<-merge(BigBio, BigBioSpMatrix, by="Plot")

BigBio$Species<-if_else(BigBio$Species=="Achillea millefolium(lanulosa)", "Achillea millefolium", BigBio$Species)

BigBio<-BigBio%>%
  filter(Species %in% species.names$Species.names) %>%
  select(Plot, Year, Species, NumSp, species_biomass_m2)
print("##############################    orig    ##########################")
print(str(BigBio))

BigBio <- merge(BigBio,species.traits, by.x="Species", by.y="Species")

# Add an 'id' column to facilitate cforest analysis
BigBio$id <- seq_along(BigBio[,1])

print("##############################    Merged    ##########################")
print(str(BigBio))

#bigbio.mono<-BigBio[BigBio$Plot==151,] #%>% 
bigbio.mono<-BigBio[BigBio$NumSp==1,] #%>% 

bigbio.mix<-BigBio[BigBio$NumSp==16,]
print("##############################    Monoculture    ##########################")
#print(bigbio.mono)
print(str(bigbio.mono))
print("##############################    Mix    ##########################")
print(str(bigbio.mix))

#BigBio<-BigBio%>%
#  arrange(Species)%>%
#  pivot_wider(names_from=Species, values_from=species_biomass_m2)
#print("##############################    arranged    ##########################")
#print(head(BigBio))
#
#BigBio<-BigBio%>%
#  select(-c("Petalostemum candidum","Petalostemum villosum", "Solidago rigida"))
#BigBio<-as.data.frame(BigBio)
#for(i in 22:39){
#  BigBio[,i]<-if_else(BigBio[,i-18]==1&is.na(BigBio[,i]),0,BigBio[,i])
#  BigBio[,i]<-ifelse(BigBio[,i-18]==0 & is.na(BigBio[,i]) == FALSE,NA,BigBio[,i])
#}
#
#BigBio<-BigBio%>%
#  select(-c(Achmi:Sornu))%>%
#  pivot_longer(cols=c(4:21), names_to="Species", values_to="species_biomass_m2")
#
#count.table<-aggregate(data = BigBio,                # Applying aggregate
#                          Species ~ Plot,
#                          function(Species) length(unique(Species)))
#
#BigBio<-na.omit(BigBio)
#
##print(dim(BigBio))
##print(str(BigBio))
##print(summary(BigBio))
#bigbio.mono<-BigBio[BigBio$NumSp==1,]
#
#bigbio.mono<-bigbio.mono%>%
#  select(-c(Plot, NumSp))
#
#names(bigbio.mono)[3]<-"mono_biomass_m2"
##print(summary(bigbio.mono))
#
#bigbio.mix<-BigBio[BigBio$NumSp!=1,]
##print(summary(bigbio.mix))
#
#bigbio.all<-left_join(bigbio.mono, bigbio.mix)
#
#bigbio.all<-bigbio.all[!duplicated(bigbio.all[c("Year","Plot","Species","species_biomass_m2")]),]
#bigbio.all<-na.omit(bigbio.all)
#bigbio.all$experiment<-"Cedar Creek"
#bigbio.all$site<-"BigBio"
#bigbio.all$block<-"A"
#names(bigbio.all)[1:7]<-c("year","species", "mono_biomass_m2", "plot", "sowndiv", "species_biomass_m2", "experiment")
#
#bigbio.all$year<-bigbio.all$year-2000
#range(bigbio.all$mono_biomass_m2)
##print(summary(bigbio.all))
##print(head(bigbio.all))
##data.all<-rbind(jena.main.all, bigbio.all)


###########################################
############### cforest ###################
###########################################

fx_cforest_single_condition <- function(model,NoSpp){
	model <- model %>%
            filter(NumSp == NoSpp) %>%
            select(-Species, -Plot, -Year, -NumSp)
        
	# Cforest analysis
	set.seed(1987)
        train <- model %>% sample_frac(.70)
        test <- anti_join(model, train, by = 'id')
        
        train <- train %>% select(-id)
        test <- test %>% select(-id)
        
        rf <- cforest(species_biomass_m2 ~ .,
                      data = as.data.frame(train),
		      control = cforest_unbiased(mtry = 2, ntree = 501))
		      
	print(rf)
        pred <- data.frame(pred = predict(rf, newdata=test,OOB=TRUE))
	test_res = test$species_biomass_m2
	corr <- cor(pred, test_res)[[1]] #, method = "kendall" #pearson (default)
	var <- round((corr^2), 2) #, method = "kendall" #pearson (default)
        title = paste("Correlation: ", corr,sep = "")

        ##Conditional permutation importance:
	# We used conditional = TRUE if predictor variables are correlated
        CPI <- varimp(rf, conditional = TRUE)
        #CPI <- varimp(rf, conditional = FALSE)
	print(head(CPI))

        #Prepare data frame for plotting
        rf_df <- as.data.frame(CPI)
        rf_df$varnames <- rownames(rf_df)
        rownames(rf_df) <- NULL
#	if(grepl("^Grass3_PCA",modelName)){
#	rf_df$type <- ifelse(rf_df$varnames == "PC2score","Mixed","Resource related")}
#	else{rf_df$type <- ifelse(rf_df$varnames %in% sizeT,"Size related","Resource related")}
        rf_df <- arrange(rf_df,desc(CPI))
        rf_df$var_categ <- c(1: dim(rf_df)[1])
	cond <- NoSpp
	rf_df$condition <- rep(cond,times=dim(rf_df)[1])
	rf_df <- rf_df %>% 
		mutate(CPI = replace(CPI, which(CPI<0), NA) ) %>%
		mutate(sCPI = (CPI*var/sum(CPI,na.rm=TRUE)))
		#mutate(sCPI = (CPI*corr/sum(CPI,na.rm=TRUE)))
        #print(rf_df)
	#print(paste0("var=",var,"; sum of scaled importance=",sum(rf_df$CPI,na.rm=TRUE)))
	#print(paste0("corr=",corr,"; sum of scaled importance=",sum(rf_df$CPI,na.rm=TRUE)))
        return(rf_df)
}

cforest_mono <- fx_cforest_single_condition(BigBio,1)
write.csv(cforest_mono, "~/Documents/SimNet/R/CedarCreek/cforest_mono.csv", row.names=FALSE)

cforest_mix <- fx_cforest_single_condition(BigBio,16)
write.csv(cforest_mix, "~/Documents/SimNet/R/CedarCreek/cforest_mix.csv", row.names=FALSE)

#Results:
#NULL
#
#	 Random Forest using Conditional Inference Trees
#
#Number of trees:  501
#
#Response:  species_biomass_m2
#Inputs:  Area_of_leaf_blade_cm2, P.A, P_A.L, SLA_cm2.g, Seed_weight_.g., height_.m.
#Number of observations:  1162
#
#Area_of_leaf_blade_cm2                    P.A                  P_A.L
#              84.03091               67.55114               70.53475
#             SLA_cm2.g        Seed_weight_.g.             height_.m.
#              44.53847              149.39597              120.18692
#
#	 Random Forest using Conditional Inference Trees
#
#Number of trees:  501
#
#Response:  species_biomass_m2
#Inputs:  Area_of_leaf_blade_cm2, P.A, P_A.L, SLA_cm2.g, Seed_weight_.g., height_.m.
#Number of observations:  5095
#
#Area_of_leaf_blade_cm2                    P.A                  P_A.L
#              18.83767               35.96540               23.91159
#             SLA_cm2.g        Seed_weight_.g.             height_.m.
#              36.70641               45.60213               22.15710

