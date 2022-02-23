library(ggbiplot)
library("party")
library("pls")
set.seed(42)

meta = seq(80, 100)
iso = seq(180, 200)

fx_PCA_all_models <- function(modelName,model){
	print(paste0("PCA: ",modelName))
        model <- model %>%
            mutate_if(is.character, as.factor) %>%
            select(-Productivity,-SpeciesID, -Ninitial, -Stage, -Rep, -Year, -Biomass) 

        model.pca <- prcomp(model, center=TRUE, scale.=TRUE)
        #print(summary(model.pca))
        #print(str(model.pca))
        #print(model.pca$center)
        #print(model.pca$rotation)
        
#	ggbiplot(mtcars.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,  labels=rownames(mtcars), groups=mtcars.country) +
        ggbiplot(model.pca, scale=0) +
		scale_colour_manual(name="Origin", values= c("forest green", "red3", "dark blue")) +
		ggtitle(paste0("PCA of ", modelName)) +
	       	theme_minimal() +
	       	theme(legend.position = "bottom")

        ggsave(paste0(tmp_dir,"/PCA/PCA_",modelName,".pdf"))
        while (!is.null(dev.list()))  dev.off()
        return(model.pca)
}

##Principal components regression
#fx_PCR <- function(modelName,model){
#	print(paste0("PCR: ",modelName))


