library(ggbiplot)
library("pls")
library("factoextra")
library(FactoMineR)
library(corrplot)
library(ggcorrplot)
set.seed(1987)

meta = seq(80, 100)
iso = seq(180, 200)

fx_cor_plot <- function(modelName,model){
# Plot the correlation between all the traits of a model 
	set.seed(1987)
        model <- model %>%
            mutate_if(is.character, as.factor) %>%
            select(-Productivity,-SpeciesID, -Ninitial, -Stage, -Rep, -Year, -Biomass) 
        # Compute a correlation matrix
        corr <- round(cor(model), 1)
	# Compute a matrix of correlation p-values
	p.mat <- cor_pmat(model)
	ggcorrplot(corr)
	ggcorrplot(corr, method = "circle", hc.order = TRUE, type = "upper",
	ggtheme = ggplot2::theme_minimal, lab = TRUE,
	p.mat = p.mat, insig = "blank") + # Leave blank on no significant coefficient
	ggtitle(paste0("Correlation between traits of ", modelName)) + 
	theme(text = element_text(size = 16))
        ggsave(paste0(tmp_dir,"/PCA/corPlot_",modelName,".pdf"))
        while (!is.null(dev.list()))  dev.off()
        return()
}

fx_PCA <- function(modelName,model){
# Perform a Principal Components Analysis for a given model and produce a biplot
	set.seed(1987)
	print(paste0("PCA: ",modelName))
        model <- model %>%
            mutate_if(is.character, as.factor) %>%
            select(-Productivity,-SpeciesID, -Ninitial, -Stage, -Rep, -Year, -Biomass) 

        model.pca <- prcomp(model, center=TRUE, scale.=TRUE)
        
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
#Correlation plot: correlations variables - dimensions
	print("var explained")
	var_exp <- (summary(PCAresult))$importance[2,]
	print(var_exp) 
	var <- get_pca_var(PCAresult)
	df <-var$cor
	#we only select the first 3 dimensions (PCscores)
	df <- df[,1:3]
	num_dim <- dim(df)[2]
	# Proportion of variation explained by each component
	var_exp <- round(var_exp[1:num_dim], 2)
	colnames(df)[colnames(df) == c("Dim.1","Dim.2","Dim.3")] <- 
		case_when(grepl("^G",modelName) ~c('LES1','Size/Growth','Spacing'), 
			  grepl("^F",modelName) ~c('LES2','maxHeight','woodDensity'))
	colnames(df) <- paste0(colnames(df),"\n (",var_exp,")")  
	
	pdf(file = paste0(tmp_dir,"/PCA/Corr_",modelName,".pdf"), family = "sans",height = 7, width = 5)
	corrplot(df, col=colorRampPalette(c("#B2182B","#FDDBC7","#2166AC"))(6), 
		 addCoef.col = "black", details=FALSE, mar = c(0,0,0,0), type = "lower",
		 #title = "Correlation of variables with each dimension", 
		 cex.col=2.5, cex.var=3.5, digits=1, cl.cex = 1.25)
	dev.off()
	return()
}
