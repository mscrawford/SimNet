library(ggbiplot)
library(aroma.light)
set.seed(42)

meta = seq(80, 100)
iso = seq(180, 200)

fx_wPCA_all_models <- function(modelName,model,NoSpp,stage){
	print(modelName)
        model <- model %>%
            filter(Ninitial == NoSpp,
                   Year %in% stage) %>%
            mutate_if(is.character, as.factor) %>%
            select(-Productivity,-SpeciesID, -Ninitial, -Stage, -Rep, -Year)
        
        w <- model$Biomass
        model <- model %>%
          select(-Biomass)
        
        model.pca <- wpca(as.matrix(model), w=w, center=TRUE, scale.=TRUE)
        print(model.pca$d)
        
        # ggbiplot(model.pca)
        # ggsave(paste0(tmp_dir,"/wPCA/",modelName,".pdf"))
        # while (!is.null(dev.list()))  dev.off()
        return(model.pca)
}
