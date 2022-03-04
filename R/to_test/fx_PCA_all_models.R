library(ggbiplot)
library("pls")
library("factoextra")
library(FactoMineR)
library(corrplot)
set.seed(1987)

meta = seq(80, 100)
iso = seq(180, 200)

fx_PCA <- function(modelName,model){
# Perform a Principal Components Analysis for a given model and produce a biplot
	print(paste0("PCA: ",modelName))
        model <- model %>%
            mutate_if(is.character, as.factor) %>%
            select(-Productivity,-SpeciesID, -Ninitial, -Stage, -Rep, -Year, -Biomass) 

        model.pca <- prcomp(model, center=TRUE, scale.=TRUE)
        #print(summary(model.pca))
        #print(str(model.pca))
        #print(model.pca$center)
        #print(model.pca$rotation)
        
	#ggbiplot(mtcars.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,  labels=rownames(mtcars), groups=mtcars.country) +
        ggbiplot(model.pca, scale=0) +
		scale_colour_manual(name="Origin", values= c("forest green", "red3", "dark blue")) +
		ggtitle(paste0("PCA of ", modelName)) +
	       	theme_minimal() +
	       	theme(text = element_text(size = 16))

        ggsave(paste0(tmp_dir,"/PCA/PCA_",modelName,".pdf"))
        #ggsave(paste0(tmp_dir,"/PCA/PCA_",modelName,".pdf"), width=9, height=7, dpi=300)
        while (!is.null(dev.list()))  dev.off()
        return(model.pca)
}

fx_get_coordinates <- function(modelID,PCAresult,fileName){
	ind <- get_pca_ind(PCAresult)
	df <- ind$coord
	modelID$PC1score <- df[,1]
	modelID$PC2score <- df[,2]
	modelID$PC3score <- df[,3]
	saveRDS(modelID,file=fileName)
	return()
}

fx_correlation_plot <- function(modelName,PCAresult){
#Correlation plot to highlight most contributing variables to each dimension
	var <- get_pca_var(PCAresult)
	pdf(file = paste0(tmp_dir,"/PCA/Corr_",modelName,".pdf"), family = "sans")
	p <- corrplot(var$contrib, is.corr=FALSE,col=colorRampPalette(c("#B2182B","#FDDBC7","#2166AC"))(100),
		 title = "Contribution of variables to each dimension", 
		 mar = c(0,0,1,0), tl.cex = 1.5, , font.main = 1, cex.axis = 1.5, cex.main = 1.7, number.digits = 1)
	dev.off()
	return(p)
}
