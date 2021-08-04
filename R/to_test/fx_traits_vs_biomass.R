library(ggplot2)
library(viridis) 

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
      
      jf1 <- unlist(model[trait1])
      jf1 <- median(jf1)*0.003
      jf2 <- unlist(model[trait2])
      jf2 <- median(jf2)*0.003

  p1 <- ggplot() +
    geom_point(data = model.32,
               aes_string(x = trait1,
                          y = trait2,
                          color = "Biomass", alpha = 0.7,
                          size = "Biomass"),
               position=position_jitter(h=jf1, w=jf2)) +
    geom_point(data = model.032, shape = 4,
               aes_string(x = trait1,
                          y = trait2),
               position=position_jitter(h=jf1, w=jf2)) +
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
