library(ggplot2)
library(forcats)
data<-read.csv("./R/Ceballos-Nunez.et.al-DataForPredictionPlot.csv")

red <- "#e41a1cff"
blue <- "#377eb8ff"

data$Monoculture.Mixture<-factor(data$Monoculture.Mixture, levels=c("Monoculture", "Mixture"))
pred.plot<-ggplot(data, aes(x=Traits, y=Importance))+
  facet_grid(Monoculture.Mixture~Strong.Weak)+
  geom_bar(stat="identity", color="black", aes(fill=Trait.Type))+
  scale_fill_manual(values=c(red, blue), guide="none")+
  ylim(0,0.9)+
  theme_bw() +
  theme(text = element_text(size = 30),legend.position = "none")
pred.plot

ggsave("PredictionPlot-11.1.2022.png", dpi=300)
