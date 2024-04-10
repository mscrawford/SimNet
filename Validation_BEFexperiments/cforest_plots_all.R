library(party)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(ggpubr)
library(ggrepel)

set.seed(1987)
#Rcolorbrewer colors:
red <- "#e41a1cff"
blue <- "#377eb8ff"
purple <- "#984ea3ff"

base_dir          <- setwd("../")
val_dir          <- paste0(base_dir, "/Validation_BEFexperiments/")
tmp_dir           <- paste0(val_dir, "tmp/")
cache_dir           <- paste0(tmp_dir, "cache/")

#jena_mono <-read.csv(paste0(cache_dir,"cforest_Jena_mono.csv"), sep=",")
#jena_mix <-read.csv(paste0(cache_dir,"cforest_Jena_mix.csv"), sep=",")
#jena_Fons_mix <-read.csv(paste0(cache_dir,"cforest_Jena_Fons_mix.csv"), sep=",")
jena_Fons_mix_2t <-read.csv(paste0(cache_dir,"cforest_Jena_Fons_mix_2t.csv"), sep=",")
jena_Fons_mono_2t <-read.csv(paste0(cache_dir,"cforest_Jena_Fons_mono_2t.csv"), sep=",")
#jena_dom_mono <-read.csv(paste0(cache_dir,"cforest_Jena_Dom_mono.csv"), sep=",")
#jena_dom_mix <-read.csv(paste0(cache_dir,"cforest_Jena_Dom_mix.csv"), sep=",")
cedar_mono <-read.csv(paste0(cache_dir,"cforest_cedar_mono.csv"), sep=",")
cedar_mix <-read.csv(paste0(cache_dir,"cforest_cedar_mix.csv"), sep=",")
cedar_2t_mono <-read.csv(paste0(cache_dir,"cforest_cedar_2t_mono.csv"), sep=",")
cedar_2t_mix <-read.csv(paste0(cache_dir,"cforest_cedar_2t_mix.csv"), sep=",")
cedar_mono_mean <-read.csv(paste0(cache_dir,"cforest_cedar_mono_mean.csv"), sep=",")
cedar_mix_mean <-read.csv(paste0(cache_dir,"cforest_cedar_mix_mean.csv"), sep=",")
cedar_2t_mono_mean <-read.csv(paste0(cache_dir,"cforest_cedar_2t_mono_mean.csv"), sep=",")
cedar_2t_mix_mean <-read.csv(paste0(cache_dir,"cforest_cedar_2t_mix_mean.csv"), sep=",")
#sardinilla_mono <-read.csv(paste0(cache_dir,"cforest_sardinilla_mono.csv"), sep=",")
#sardinilla_mix <-read.csv(paste0(cache_dir,"cforest_sardinilla_mix.csv"), sep=",")
#sardinilla2_mono <-read.csv(paste0(cache_dir,"cforest_sardinilla2_mono.csv"), sep=",")
#sardinilla2_mix <-read.csv(paste0(cache_dir,"cforest_sardinilla2_mix.csv"), sep=",")
#sardinilla2_2012_mono <-read.csv(paste0(cache_dir,"cforest_sardinilla2_2012_mono.csv"), sep=",")
#sardinilla2_2012_mix <-read.csv(paste0(cache_dir,"cforest_sardinilla2_2012_mix.csv"), sep=",")
sardinilla_mono <-read.csv(paste0(cache_dir,"cforest_sardinilla_mono_.csv"), sep=",")
sardinilla_mix <-read.csv(paste0(cache_dir,"cforest_sardinilla_mix_.csv"), sep=",")
sardinilla2_mono <-read.csv(paste0(cache_dir,"cforest_sardinilla2_mono_.csv"), sep=",")
sardinilla2_mix <-read.csv(paste0(cache_dir,"cforest_sardinilla2_mix_.csv"), sep=",")
sardinilla2_2012_mono <-read.csv(paste0(cache_dir,"cforest_sardinilla2_2012_mono_.csv"), sep=",")
sardinilla2_2012_mix <-read.csv(paste0(cache_dir,"cforest_sardinilla2_2012_mix_.csv"), sep=",")

#jena_mix$mName <- rep("Jena",times=dim(jena_mix)[1])
#jena_Fons_mix$mName <- rep("Jena_Fons",times=dim(jena_Fons_mix)[1])
jena_Fons_mix_2t$mName <- rep("Jena_Fons_2t",times=dim(jena_Fons_mix_2t)[1])
jena_Fons_mono_2t$mName <- rep("Jena_Fons_2t",times=dim(jena_Fons_mono_2t)[1])
#jena_dom_mix$mName <- rep("Jena_dom",times=dim(jena_dom_mix)[1])
cedar_mono_mean$mName <- rep("Cedar_mean",times=dim(cedar_mono_mean)[1])
cedar_mix_mean$mName <- rep("Cedar_mean",times=dim(cedar_mono_mean)[1])
cedar_2t_mono_mean$mName <- rep("Cedar_2t_mean",times=dim(cedar_2t_mono_mean)[1])
cedar_2t_mix_mean$mName <- rep("Cedar_2t_mean",times=dim(cedar_2t_mix_mean)[1])
cedar_mono$mName <- rep("Cedar",times=dim(cedar_mono)[1])
cedar_mix$mName <- rep("Cedar",times=dim(cedar_mono)[1])
cedar_2t_mono$mName <- rep("Cedar_2t",times=dim(cedar_2t_mono)[1])
cedar_2t_mix$mName <- rep("Cedar_2t",times=dim(cedar_2t_mix)[1])
sardinilla_mono$mName <- rep("Sardinilla",times=dim(sardinilla_mono)[1])
sardinilla_mix$mName <- rep("Sardinilla",times=dim(sardinilla_mix)[1])
sardinilla2_mono$mName <- rep("Sardinilla st",times=dim(sardinilla2_mono)[1])
sardinilla2_mix$mName <- rep("Sardinilla st",times=dim(sardinilla2_mix)[1])
sardinilla2_2012_mono$mName <- rep("Sardinilla st >2012",times=dim(sardinilla2_2012_mono)[1])
sardinilla2_2012_mix$mName <- rep("Sardinilla st >2012",times=dim(sardinilla2_2012_mix)[1])

#all_d <- rbind(jena_Fons_mix)
#all_d <- rbind(jena_Fons_mix_2t, jena_mix, jena_Fons_mix)
#all_d <- rbind(jena_Fons_mix_2t, jena_Fons_mono_2t)
#all_d <- rbind(cedar_2t_mono, cedar_2t_mix, cedar_mono, cedar_mix)
#all_d <- rbind(jena_mix, jena_Fons_mix, jena_dom_mix, cedar_mono, cedar_mix)
#print(head(all_d))
resvar = "Biomass"
fx_plot_all <- function(df,resvar,plot_name){	
# Function to plot the random forest results for all four conditions (for all models)
	# include function-dominance correlation in model name
	df$mNameFDC <-  paste0(df$mName,'\n (',if(resvar=="Biomass"){df$funcdom}else{df$funcdom_p},')')
	p <- ggplot(df, aes(x=varnames, y=sCPI, fill=type)) +
	#p <- ggplot(df, aes(x=reorder(varnames,typen), y=sCPI, fill=type)) +
	    geom_bar(position='dodge',stat='identity') +
	#            geom_text(aes(label=scientific(CPI, digits = 2),size=10)
	#    		  ,position = position_dodge(width = 1)
	#                      ,hjust=0,vjust=0.2,color="black",size=2
	#                      ,show.legend = FALSE,angle = 90) +
	    #ggtitle(paste0("Trait importance in determining ",resvar)) +
	    ylab("Proportion of variance explained") +
	    #ylab("Conditional permutation importance, \n scaled to r^2") +
	    #ylab("Conditional permutation importance, \n scaled to correlation") +
	    xlab("Traits") +
	    scale_fill_manual(name = "Trait type",
			      values = c("Resource related" = red,
					  "Size related" = blue,
					  "Mixed" = purple),
			      labels = c("Mixed","Resource", "Size")) +
	    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5)) + 
	    facet_grid(condition ~ mName, scales = "free_x") +
	    #facet_grid(reorder(condition,modeln) ~ reorder(mNameFDC, -if(resvar=="Biomass"){funcdom}else{funcdom_p}), scales = "free_x") +
	    theme_bw() +
	    theme(text = element_text(size = 26),legend.position = "bottom",
	    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
	    panel.grid.minor = element_blank(),
	    panel.grid.major = element_blank(),
	    plot.title = element_text(hjust = 0.5))

	ggsave(file=paste0(tmp_dir,plot_name,".png")
	       , width=13, height=9, dpi=300
	)
	while (!is.null(dev.list()))  dev.off()
	return(p)
}

#plot_name <- "cforest_grid_JenaFons"
all_d <- rbind(cedar_mono, cedar_mix, cedar_2t_mono, cedar_2t_mix, jena_Fons_mix_2t, jena_Fons_mono_2t, sardinilla2_mono, sardinilla2_mix, sardinilla2_2012_mono, sardinilla2_2012_mix)
plot_name <- "cforest_grid_experiments_"
fx_plot_all(all_d,resvar,plot_name)

cedars <- rbind(cedar_mono, cedar_mix, cedar_2t_mono, cedar_2t_mix, cedar_mono_mean, cedar_mix_mean, cedar_2t_mono_mean, cedar_2t_mix_mean)
plot_name <- "cforest_cedars"
fx_plot_all(cedars,resvar,plot_name)

sard <- rbind(cedar_2t_mono, cedar_2t_mix, jena_Fons_mix_2t, jena_Fons_mono_2t, sardinilla_mono, sardinilla_mix, sardinilla2_mono, sardinilla2_mix)
plot_name <- "cforest_sard"
fx_plot_all(sard,resvar,plot_name)
