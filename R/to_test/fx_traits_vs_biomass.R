library(ggplot2)
library(viridis) 
library(reshape2)
library(data.table)

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

#fx_traits_vs_biomass_jitter <- function(plotName,df,NoSpp,stage,driver,response,xlab,ylab){
fx_traits_vs_biomass_jitter <- function(plotName,df,NoSpp,stage,trait1,trait2,xlab,ylab){
	is_productivity = grepl("_P$",plotName)
	response <- if(is_productivity){"Productivity"}else{"Biomass"}
	df <- df %>%
	filter(Ninitial == NoSpp,
	       Year %in% stage) %>%
	select(all_of(-Rep,-Ninitial,-Year,-Stage))
	df <- aggregate(. ~ SpeciesID, data = df, mean) #Calculate mean response per specie

	df.32 <- df %>%
	#mutate(Biomass = scales::rescale(Biomass, to = c(0, 100))) %>%
	filter(Biomass > 0) %>%
	mutate(Biomass = log(Biomass)) %>%
	#mutate(Productivity = scales::rescale(Productivity, to = c(0, 100))) %>%
	filter(Productivity > 0) %>%
	mutate(Productivity = log(Productivity))
	
	df.032 <- df %>%
	filter(Biomass == 0) %>%
	filter(Productivity == 0)
	
	jf_1 <- unlist(df[trait1])
	jf1 <- median(jf_1)*0.06
	jf_2 <- unlist(df.32[trait2])
	jf2 <- median(jf_2)*0.06
	xmax <- max(jf_1) + (max(jf_1) *0.06)
	ymin <- min(jf_2) #*0.06
	ymax <- max(jf_2) #*0.06
	ymean <- mean(jf_2)#(ymax+ymin)/2

	modelName <- case_when(grepl("^G1",plotName) ~'Grass 1',
			       grepl("^G2",plotName) ~'Grass 2',
			       grepl("^G3",plotName) ~'Grass 3',
			       grepl("^F1",plotName) ~'Forest 1',
			       grepl("^F2",plotName) ~'Forest 2',
			       grepl("^D",plotName) ~'Dryland')
	print(modelName)
  p1 <- ggplot() +
    geom_point(data = df.32,
               aes_string(x = trait1,
                          y = trait2,
                          color = response, alpha = 0.8,
                          size = response),
               position=position_jitter(h=jf2, w=jf1)) +
    ylim(ymin, (ymax+(ymax*0.06))) +
    xlim(NA, (xmax+(xmax*0.06))) +
    ggtitle(modelName) + 
    geom_point(data = df.032, shape = 4,
               aes_string(x = trait1,
                          y = trait2),
               position=position_jitter(h=jf2, w=jf1)) +
    scale_size_continuous(name = "Prop.",
                          breaks = c(ymin,ymean,ymax),#c(1,50,100),
#                          limits = c(ymin-0.06, ymax*0.06),
                          labels = c("Low","Mean","High"),
                          range = c(0, 6)) +
    scale_colour_viridis(breaks = c(ymin,ymean,ymax),labels = c("Low","Mean","High")
			 #,limits = c(ymin*0.06, ymax*0.06)
			 ) + #direction = -1) +
    labs(x = xlab, y = ylab#)+
	 #size = paste0("Log mean \n",response),
	 , color = paste0("Log mean \n",response)) +
    theme_bw() +
    theme_classic() +
    theme(aspect.ratio = 0.5,text = element_text(size = 30),
    plot.title=element_text(hjust=0.5, vjust=0.5), legend.position = c(0.8, 0.57), 
    legend.text = element_text(size=20), legend.title = element_text(size=20)) +
	   if(grepl("G1C1|DC1_P",plotName)){guides(alpha="none", size="none")}else{guides(alpha="none", size="none", color="none")} 
  filename <- paste0(plotName,".png")
  path <- paste0(tmp_dir,"/traits_vs_biomass/")
  ggsave(filename = filename, path = path, plot = p1
         ,height = 13.5, width = 22, units = "cm")
  return(p1)
}

fx_traits_vs_biomass_jitter_I <- function(plotName,df,NoSpp,stage,trait1,trait2,xlab,ylab){
#	source("http://www.openintro.org/stat/data/arbuthnot.R")
#	print(head(arbuthnot))
	is_productivity = grepl("_P$",plotName)
	response <- if(is_productivity){"Productivity"}else{"Biomass"}
	df <- df %>%
	filter(Ninitial == NoSpp,
	       Year %in% stage) %>%
	select(all_of(-Rep,-Ninitial,-Year,-Stage,-id))
	df <- aggregate(. ~ SpeciesID, data = df, mean) #Calculate mean response per specie

	df.32 <- df %>%
	#mutate(Biomass = scales::rescale(Biomass, to = c(0, 100))) %>%
	filter(Biomass > 0) %>%
	mutate(Biomass = log(Biomass)) %>%
	#mutate(Productivity = scales::rescale(Productivity, to = c(0, 100))) %>%
	filter(Productivity > 0) %>%
	mutate(Productivity = log(Productivity)) %>%
	select(all_of(response,trait1,trait2))

print(head(df.32))
        dfmelt <- melt(df.32, id.vars = trait1, variable.name = 'Traits', value.name = response)
print(head(dfmelt))
	
	df.032 <- df %>%
	filter(Biomass == 0) %>%
	filter(Productivity == 0)
	
	jf_1 <- unlist(df[trait1])
	jf1 <- median(jf_1)*0.06
	jf_2 <- unlist(df.32[response])
	jf2 <- median(jf_2)*0.06
	xmax <- max(jf_1) + (max(jf_1) *0.06)
	ymin <- min(jf_2) #*0.06
	ymax <- max(jf_2) #*0.06
	ymean <- mean(jf_2)#(ymax+ymin)/2

	modelName <- case_when(grepl("^G1",plotName) ~'Grass 1',
			       grepl("^G2",plotName) ~'Grass 2',
			       grepl("^G3",plotName) ~'Grass 3',
			       grepl("^F1",plotName) ~'Forest 1',
			       grepl("^F2",plotName) ~'Forest 2',
			       grepl("^D",plotName) ~'Dryland')
	print(modelName)
  p1 <- ggplot() +
    geom_point(data = df.32,
               aes_string(x = trait1,
                          y = response,
                          color = "blue", alpha = 0.8),
               position=position_jitter(h=jf2, w=jf1)) +
    geom_point(data = df.32,
               aes_string(x = trait2,
                          y = response,
                          color = "red", alpha = 0.8),
               position=position_jitter(h=jf2, w=jf1)) +
    ylim(ymin, (ymax+(ymax*0.06))) +
#    xlim(NA, (xmax+(xmax*0.06))) +
    ggtitle(modelName) + 
    geom_point(data = df.032, shape = 4,
               aes_string(x = trait1,
                          y = response),
               position=position_jitter(h=jf2, w=jf1)) +
#    scale_size_continuous(name = "Prop.",
#                          breaks = c(ymin,ymean,ymax),#c(1,50,100),
##                          limits = c(ymin-0.06, ymax*0.06),
#                          labels = c("Low","Mean","High"),
#                          range = c(0, 6)) +
#    scale_colour_viridis(breaks = c(ymin,ymean,ymax),labels = c("Low","Mean","High")
#			 #,limits = c(ymin*0.06, ymax*0.06)
#			 ) + #direction = -1) +
    labs(x = "Trait", y = paste0("Log mean \n",response)#)+
	 #size = paste0("Log mean \n",response),
	 , color = "Trait") +
    theme_bw() +
    theme_classic() +
    theme(aspect.ratio = 0.5,text = element_text(size = 30),
    plot.title=element_text(hjust=0.5, vjust=0.5), legend.position = c(0.8, 0.57), 
    legend.text = element_text(size=20), legend.title = element_text(size=20)) +
	   guides(alpha="none", size="none")
  filename <- paste0(plotName,".png")
  path <- paste0(tmp_dir,"/traits_vs_biomass/")
  ggsave(filename = filename, path = path, plot = p1
         ,height = 13.5, width = 22, units = "cm")
  return(p1)
}

fx_prepare_df  <- function(modelName,data,trait1,trait2,lab1,lab2){
	dfmono <- data %>%
	filter(Ninitial == 1, Year %in% meta) %>%
	select(SpeciesID, Biomass, trait1, trait2)
    dfmono$type_bm <- rep("Monoculture",times=dim(dfmono)[1])

	dfmix <- data %>%
	filter(Ninitial == 32, Year %in% meta) %>%
	select(SpeciesID, Biomass, trait1, trait2) #%>%
	dfmix <- aggregate(. ~ SpeciesID, data = dfmix, mean) #Calculate mean response per specie
	dfmix$type_bm <- rep("Mixture",times=dim(dfmix)[1])

	df <- rbind(dfmono,dfmix)  %>%
	setnames(old=c(trait1,trait2), new=c(lab1,lab2))

	df$Model <- rep(modelName,times=dim(df)[1])
#        df <- reshape2::melt(df, id="SpeciesID")
#	print("#######     No melt      #############")
	return(df)
}
#"Monoculture \n Biomass""Mixture \n Biomass"

#To do, after melt
#fx_plot_all <- function(df){
#	print(paste0("####################",spec(df)))
#	p <- ggplot(df, aes(x=SpeciesID, y=value)) + geom_point() +
#	facet_wrap(~variable, scales = "free")
#	theme(aspect.ratio = 0.5,text = element_text(size = 30),
#    plot.title=element_text(hjust=0.5, vjust=0.5), legend.position = c(0.8, 0.57), 
#    legend.text = element_text(size=20), legend.title = element_text(size=20)) +
#        guides(alpha="none", size="none")
#        path <- paste0(tmp_dir,"/traits_vs_biomass/")
#        ggsave(filename = "Figure_3.png", path = path, plot = p
#        ,height = (13.5*6), width = (22*4), units = "cm")
#return(p)
#}

fx_biomass_vs_traits <- function(df){

	df.32 <- df %>%
	#mutate(Biomass = scales::rescale(Biomass, to = c(0, 100))) %>%
	filter(Biomass > 0) %>%
	mutate(Biomass = log(Biomass))
print(df.32)
	df.032 <- df %>%
	filter(Biomass == 0) 
	
#	jf_1 <- unlist(df[trait1])
#	jf1 <- median(jf_1)*0.06
#	jf_2 <- unlist(df.32[Biomass])
#	jf2 <- median(jf_2)*0.06
#	xmax <- max(jf_1) + (max(jf_1) *0.06)
#	ymin <- min(jf_2) #*0.06
#	ymax <- max(jf_2) #*0.06
#	ymean <- mean(jf_2)#(ymax+ymin)/2

  p1 <- ggplot() +
	  geom_point(data = df.32, 
		     aes(x = t1,
			 y = Biomass
			 )) +
	  geom_point(data = df.32, 
		     aes(x = t2,
			 y = Biomass
			 )) +
	  facet_grid(Ninitial+lab1~Model, scales="free") +
	  stat_smooth(color="red", method="loess") +

#    geom_point(data = df.32,
#               aes_string(x = trait2,
#                          y = Biomass,
#                          color = "red", alpha = 0.8),
#               position=position_jitter(h=jf2, w=jf1)) +
#    ylim(ymin, (ymax+(ymax*0.06))) +
##    xlim(NA, (xmax+(xmax*0.06))) +
#    ggtitle(modelName) + 
    geom_point(data = df.032, shape = 4, aes(x = t1, y = Biomass)) +
##    scale_size_continuous(name = "Prop.",
##                          breaks = c(ymin,ymean,ymax),#c(1,50,100),
###                          limits = c(ymin-0.06, ymax*0.06),
##                          labels = c("Low","Mean","High"),
##                          range = c(0, 6)) +
##    scale_colour_viridis(breaks = c(ymin,ymean,ymax),labels = c("Low","Mean","High")
##			 #,limits = c(ymin*0.06, ymax*0.06)
##			 ) + #direction = -1) +
#    labs(x = "Trait", y = "Log mean \n Biomass"#)+
#	 #size = "Log mean \n Biomass",
#	 , color = "Trait") +
    theme_bw() +
    theme_classic() +
    theme(aspect.ratio = 0.5,text = element_text(size = 30),
    plot.title=element_text(hjust=0.5, vjust=0.5), legend.position = c(0.8, 0.57), 
    legend.text = element_text(size=20), legend.title = element_text(size=20)) +
	   guides(alpha="none", size="none")
  path <- paste0(tmp_dir,"/traits_vs_biomass/")
  ggsave(filename = "Figure_3.png", path = path, plot = p1
         ,height = (13.5*6), width = (22*4), units = "cm")
  return(p1)
}

#fx_single_plot  <- function(modelName,df,t1,t2,lab1,lab2){
#	df <- data.frame(df) %>%
#	filter(Ninitial == 1 | Ninitial == 32,
#	       Year %in% meta) %>%
#	select(-Rep, -Stage,-Year,-id) %>%
#	mutate(Ninitial = recode(Ninitial, "1" = "Monoculture",
#				 "32" = "Mixture"))
#       # melt(df, id=SpeciesID)
#	print(dim(df))
#	aggregate(. ~ SpeciesID, data = df, mean) #Calculate mean response per specie
#	print(dim(df))
#	#df <- df %>%
#	#select(Ninitial, Biomass, t1, t2) 
##	df$Model <- rep(modelName,times=dim(df)[1])
##	setnames(df, old=c(trait1,trait2), new=c("t1","t2"))
##	df$lab1 <- rep(lab1,times=dim(df)[1])
##	df$lab2 <- rep(lab2,times=dim(df)[1])
#	  facet_grid(.~Ninitial, scales="free") +
#  p1 <- ggplot() +
#	  geom_point(data = df.32, 
#		     aes(x = .data[[lab1]],
#			 y = Biomass_mono,
#			 )) +
#	  stat_smooth(color="red", method="loess") +
#    plot_vars
#
#  p2 <- ggplot() +
#	  geom_point(data = df.32, 
#		     aes(x = .data[[lab1]],
#			 y = Biomass_mix,
#			 )) +
#	  stat_smooth(color="red", method="loess") +
#    geom_point(data = df.032, shape = 4, aes(x = .data[[lab1]], y = Biomass_mix)) +
#    plot_vars
#
#  p3 <- ggplot() +
#	  geom_point(data = df.32, 
#		     aes(x = .data[[lab2]],
#			 y = Biomass_mono,
#			 )) +
#	  stat_smooth(color="red", method="loess") +
#    plot_vars
#

fx_scatter_plot <- function(condition,df,lab,plot_name){
    df <- if(condition == "mono"){df[df$type_bm=="Monoculture",]} else {df[df$type_bm=="Mixture",]}
    df.32 <- df %>%
	#mutate(Biomass = scales::rescale(Biomass, to = c(0, 100))) %>%
        filter(Biomass > 0) %>%
        mutate(Biomass = log(Biomass))
	#df.032 <- df %>%
	#filter(Biomass == 0) 
	#print(head(df.32$Model))
  p <- ggplot() +
	  geom_point(data = df.32, 
		     aes(x = .data[[lab]],
			 y = Biomass,
			 )) + {
	  if (df$Model[1] != "Dryland") stat_smooth(aes(df.32[[lab]],df.32$Biomass), fullrange = TRUE, color="red", method="loess", se=FALSE)
          } +
    #geom_point(data = df.032, shape = 4, aes(x = .data[[lab]], y = Biomass)) +
    labs(#size = "Log mean biomass", color = "Log mean biomass",
      y = "log biomass") +
    ggtitle(df$Model[1]) + 
    theme_bw() +
    theme_classic() +
    theme(aspect.ratio = 0.5,text = element_text(size = 30),
    plot.title=element_text(hjust=0.5, vjust=0.5), legend.position = c(0.8, 0.57), 
    legend.text = element_text(size=20), legend.title = element_text(size=20)) +
    guides(alpha="none", size="none")
  path <- paste0(tmp_dir,"/traits_vs_biomass/")
  ggsave(filename = paste0(plot_name,".png"), path = path, plot = p
         ,height = (13.5), width = (22), units = "cm")
	return(p)
}

