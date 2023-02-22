library(ggplot2)
library(viridis) 
set.seed(1987)

meta = seq(80, 100) #100th year: last step of Metacommunity stage 
iso = seq(180, 200) #200th year: last step of Isolation stage

fx_traits_vs_biomass <- function(plotName,model,NoSpp,stage,g_by,trait1,trait2,xlab,ylab){
      model <- model %>%
      filter(Ninitial == NoSpp,
             Year %in% stage) %>% 
      group_by_at(g_by) %>%
      summarise(Biomass = mean(Biomass))

      model.32 <- model %>%
      filter(Biomass > 0) %>%
      mutate(Biomass = log(Biomass))

      model.032 <- model %>%
      filter(Biomass == 0)
      
  p1 <- ggplot() +
    geom_point(data = model.32,
               aes_string(x = trait1,
                          y = trait2,
                          color = "Biomass", alpha = 0.7,
                          size = "Biomass")) +
    geom_point(data = model.032, shape = 4,
               aes_string(x = trait1,
                          y = trait2)) +
    scale_color_viridis() +
    labs(size = "Log mean biomass",
	 color = "Log mean biomass",
         x = xlab,
         y = ylab) +
    theme_bw() +
    theme(aspect.ratio = 1)
  
  filename <- paste0(plotName,".pdf")
  path <- paste0(tmp_dir,"/traits_vs_biomass/")
  ggsave(filename = filename
         ,path = path
         ,plot = p1
         ,height = 12
         ,width = 19
         ,units = "cm")
  return(p1)
}

fx_traits_vs_biomass_jitter <- function(plotName,model,NoSpp,stage,driver,response,xlab,ylab){
	print(head(model))
	is_productivity = grepl("_P$",plotName)
	response <- if(is_productivity){"Productivity"}else{"Biomass"}
	model <- model %>%
	filter(Ninitial == NoSpp,
	       Year %in% stage) %>%
	select(-Rep, -Ninitial)
	model <- aggregate(.~SpeciesID,data=model, mean)

	model.32 <- model %>%
	filter(Biomass > 0) %>%
	mutate(Biomass = log(Biomass)) %>%
	filter(Productivity > 0) %>%
	mutate(Productivity = log(Productivity))
	
	model.032 <- model %>%
	filter(Biomass == 0) %>%
	filter(Productivity == 0)
	
	jf_1 <- unlist(model[driver])
	jf1 <- median(jf_1)*0.06
	jf_2 <- unlist(model[response])
	jf2 <- median(jf_2)*0.06

	modelName <- case_when(grepl("^G1",plotName) ~'Grass 1',
			       grepl("^G2",plotName) ~'Grass 2',
			       grepl("^G3",plotName) ~'Grass 3',
			       grepl("^F1",plotName) ~'Forest 1',
			       grepl("^F2",plotName) ~'Forest 2',
			       grepl("^D",plotName) ~'Dryland')
	print(modelName)
  p1 <- ggplot() +
    geom_point(data = model.32,
               aes_string(y = response,
                          x = driver),
                          #color = response, alpha = 0.8,
                          #size = response
               position=position_jitter(h=jf2, w=jf1)) +
    ggtitle(modelName) + 
    geom_point(data = model.032, shape = 4,
               aes_string(x = driver,
                          y = response),
               position=position_jitter(h=jf2, w=jf1)) +
    #scale_colour_viridis() + #direction = -1) +
    labs(x = xlab, y = ylab)+
	 #size = paste0("Log mean \n",response),
	 #color = paste0("Log mean \n",response)) +
    theme_bw() +
    theme_classic() +
    theme(aspect.ratio = 0.5,text = element_text(size = 30),
    plot.title=element_text(hjust=0.5, vjust=0.5), legend.position = c(0.8, 0.57), 
    legend.text = element_text(size=20), legend.title = element_text(size=20)) +
	    if(grepl("G1C1|DC1_P",plotName)){guides(alpha="none", size="none")}else{guides(alpha="none", size="none", color="none")} 
  filename <- paste0(plotName,".pdf")
  path <- paste0(tmp_dir,"/traits_vs_biomass/")
  ggsave(filename = filename, path = path, plot = p1
         ,height = 13.5, width = 22, units = "cm")
  return(p1)
}
