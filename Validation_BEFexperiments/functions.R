library(ggplot2)

fx_plot_trait_Vs_biomass <- function(df, plot_name){
df$Year <- as.character(df$Year)
    df %>%
        gather(-Species, -Plot, -Year, -id, -NumSp, -species_biomass_m2, key= "var", value = "value") %>%
        ggplot(aes(x=value, y=log(species_biomass_m2), color=Year, shape=Year)) +
        geom_point() +
        labs(y = "log biomass") +
        geom_smooth(method="lm", fill=NA) +
    	#stat_smooth(fullrange = TRUE, color="red", method="loess", se=FALSE) +
        facet_wrap(~var, scales ="free") +
        theme_bw() +
        theme(text = element_text(size = 25), strip.text = element_text(size=25),
              legend.text = element_text(size=20), legend.title = element_text(size=20))
    ggsave(file=paste0(tmp_dir, plot_name)
           , width=18, height=14, dpi=300)
    while (!is.null(dev.list()))  dev.off()
    return()
}
