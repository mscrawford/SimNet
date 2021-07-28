library(ggbiplot)
library("party")
set.seed(42)

meta = seq(80, 100)
iso = seq(180, 200)

fx_PCA_all_models <- function(modelName,model,NoSpp,stage){
	print(modelName)
        model <- model %>%
            filter(Ninitial == NoSpp,
                   Year %in% stage) %>%
            mutate_if(is.character, as.factor) %>%
            select(-Productivity,-SpeciesID, -Ninitial, -Stage, -Rep, -Year) 

        model.pca <- prcomp(model, center = TRUE,scale. = TRUE)
        print(summary(model.pca))
        
        ggbiplot(model.pca)
        ggsave(paste0(tmp_dir,"/PCA/",modelName,".pdf"))
        while (!is.null(dev.list()))  dev.off()
        return(model.pca)
}
