#https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html
#https://jkzorz.github.io/2019/07/02/boxplots.html

library(tidyr)
library(tidyverse)
library(indicspecies)
library(dplyr)
library(phyloseq)
library(microbiome)
library(ggplot2)
groups = c("Spring", "Summer", "Autumn","Winter")

l4_indval_ps = prune_taxa(taxa_sums(l4_fungi_genus) > 0, l4_fungi_genus) 
L4_indval<-psmelt(l4_indval_ps)

L4_indval_1 = L4_indval %>% unite(season_year, season, year, sep = "_", remove = FALSE)
L4_indval_1 <- L4_indval_1[, c("Genus", "Abundance", "season","season_year","sample")]
L4_indval_1<-aggregate(Abundance ~ season+season_year+Genus, data = L4_indval_1, mean)
L4_indval_1.wide <- pivot_wider(L4_indval_1, names_from = Genus, values_from = Abundance) 
abund = L4_indval_1.wide[,3:ncol(L4_indval_1.wide)]
season = L4_indval_1.wide$season
inv = multipatt(abund, season,  control = how(nperm=999),duleg = TRUE,func = "indval.g")
summary(inv)

summarize_phyloseq(l4_fungi)

indval_winter <- subset_taxa(l4_fungi,Genus == "Clitocybe" |
                        Genus == "Phlebia" |
                          Genus == "Phlebiopsis" |
                        Genus == "Fomitopsis" | 
                        Genus == "Vuilleminia" |
                          Genus == "Diatrypella" |
                          Genus == "Hyphodontia" |
                          Genus == "Flammulina" |
                          Genus == "Tarzetta" |
                          Genus == "Tetracladium" |
                          Genus == "Myrmecridium" |
                          Genus == "Otidea" |
                          Genus == "Sporopachydermia" |
                          Genus == "Proliferodiscus")

#extract abundances
indval_winter_genus<-tax_glom(indval_winter, taxrank = "Genus")
indval_winter_genus.melt<-psmelt(indval_winter_genus)
indval_winter_genus.melt <- aggregate(indval_winter_genus.melt$Abundance,
                                      by=list(indval_winter_genus.melt$Genus), FUN=sum, na.rm=TRUE)

indval_winter <- filter_taxa(indval_winter, function(x){sum(x > 0) > 20}, prune = TRUE)                  
indval_winter.melt<-psmelt(indval_winter)
indval_winter_plot<-ggplot(indval_winter.melt, aes(x = week_str, y = log(Abundance),fill="#20639B")) + 
  geom_boxplot(stat = "boxplot") + 
   facet_wrap("Genus", ncol=4) + 
  theme_ro()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(face = 'italic'),legend.position = "none") +
  xlab("Weeks, January - December (2002-2018)") +
  ylab("Log10 Abundance (reads)") +theme(aspect.ratio = 0.5)+ 
  scale_fill_manual(values=c("#20639B","#20639B")) + 
  ggtitle("Winter")+ ylim(0, 11.8)


indval_spring <- subset_taxa(l4_fungi,Genus == "Trametes" |
                               Genus == "Metschnikowia" |
                               Genus == "Lichina" |
                               Genus == "Bjerkandera" | 
                               Genus == "Pycnoporellus" |
                               Genus == "Crocicreas" |
                               Genus == "Scoliciosporum" |
                               Genus == "Bradymyces" |
                               Genus == "Fibulomyces" |
                               Genus == "Coprinopsis" |
                               Genus == "Phomatodes" |
                               Genus == "Ypsilina"|
                               Genus == "Bagliettoa")

#extract abundances
indval_spring_genus<-tax_glom(indval_spring, taxrank = "Genus")
indval_spring_genus.melt<-psmelt(indval_spring_genus)
indval_spring_genus.melt <- aggregate(indval_spring_genus.melt$Abundance,
                                      by=list(indval_spring_genus.melt$Genus), FUN=sum, na.rm=TRUE)

indval_spring <- filter_taxa(indval_spring, function(x){sum(x > 0) > 20}, prune = TRUE)                  

indval_spring.melt<-psmelt(indval_spring)
indval_spring_plot<-ggplot(indval_spring.melt, aes(x = week_str, y = log(Abundance),fill=Phylum)) + 
  geom_boxplot(stat = "boxplot") + 
  facet_wrap("Genus", ncol=4) + 
  theme_ro()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(face = 'italic'),legend.position = "none") +
  xlab("Weeks, January - December (2002-2018)") +
  ylab("Log10 Abundance (reads)") +theme(aspect.ratio = 0.5)+ 
  scale_fill_manual(values=c("#3CAEA3","#3CAEA3")) + 
  ggtitle("Spring")+ ylim(0, 11.8)


indval_summer <- subset_taxa(l4_fungi,Genus == "Alternaria" |
                               Genus == "Cladosporium" |
                               Genus == "Limonomyces" | 
                               Genus == "Leptosphaerulina" |
                               Genus == "Claviceps" |
                               Genus == "Podosphaera"|
                               Genus == "Venturia")
#extract abundances
indval_summer_genus<-tax_glom(indval_summer, taxrank = "Genus")
indval_summer_genus.melt<-psmelt(indval_summer_genus)
indval_summer_genus.melt <- aggregate(indval_summer_genus.melt$Abundance,
                                           by=list(indval_summer_genus.melt$Genus), FUN=sum, na.rm=TRUE)

indval_summer <- filter_taxa(indval_summer, function(x){sum(x > 0) > 20}, prune = TRUE)                  
indval_summer.melt<-psmelt(indval_summer)
indval_summer_plot<-ggplot(indval_summer.melt, aes(x = week_str, y = log(Abundance),fill=Phylum)) + 
  geom_boxplot(stat = "boxplot") + 
  facet_wrap("Genus", ncol=4,nrow=3) + 
  theme_ro()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(face = 'italic'),legend.position = "none") +
  xlab("Weeks, January - December (2002-2018)") +
  ylab("Log10 Abundance (reads)") +theme(aspect.ratio = 0.5)+ 
  scale_fill_manual(values=c("#F6D55C","#F6D55C")) + 
  ggtitle("Summer")+ ylim(0, 11.8)

indval_autumn <- subset_taxa(l4_fungi,Genus == "Hypholoma" |
                               Genus == "Lepista" |
                               Genus == "Armillaria" | 
                               Genus == "Paramicrothyrium" |
                               Genus == "Catenulostroma" |
                               Genus == "Kluyveromyces" |
                               Genus == "Curvularia" |
                               Genus == "Gloiothele")

#extract abundances
indval_autumn_genus<-tax_glom(indval_autumn, taxrank = "Genus")
indval_autumn_genus.melt<-psmelt(indval_autumn_genus)
indval_autumn_genus.melt <- aggregate(indval_autumn_genus.melt$Abundance,
                                      by=list(indval_autumn_genus.melt$Genus), FUN=sum, na.rm=TRUE)

indval_autumn <- filter_taxa(indval_autumn, function(x){sum(x > 0) > 20}, prune = TRUE)                  

indval_autumn.melt<-psmelt(indval_autumn)
indval_autumn_plot<-ggplot(indval_autumn.melt, aes(x = week_str, y = log(Abundance),fill=Phylum)) + 
  geom_boxplot(stat = "boxplot") + 
  facet_wrap("Genus", ncol=4,nrow=3) + 
  theme_ro()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(face = 'italic'),legend.position = "none") +
  xlab("Weeks, January - December (2002-2018)") +
  ylab("Log10 Abundance (reads)") +theme(aspect.ratio = 0.5)+ 
  scale_fill_manual(values=c("#ED553B","#ED553B")) + 
  ggtitle("Autumn")+ ylim(0, 11.8)


ggplot(indval_winter.melt, aes(x = Genus, y = log(Abundance), fill = season)) + 
  geom_boxplot(colour = "black", position = position_dodge(0.8)) +
  theme(legend.title = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 10, face = "bold"), legend.position = "right", 
        axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 12, colour = "black"), 
        axis.title.y = element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.key=element_blank()) + 
  labs(x= "", y = "Relative Abundance (%)", fill = "Season") + 
  scale_fill_manual(values = c("white", "lightgray","darkgray", "black"))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#indval examples
indval_examples <- subset_taxa(l4_fungi,Genus == "Phlebia" |
                               Genus == "Lichina" |
                               Genus == "Alternaria" |
                               Genus == "Hypholoma" )

#extract abundances
indval_examples_genus<-tax_glom(indval_examples, taxrank = "Genus")
indval_examples_genus.melt<-psmelt(indval_examples_genus)
indval_examples_genus.melt <- aggregate(indval_examples_genus.melt$Abundance,
                                      by=list(indval_examples_genus.melt$Genus), FUN=sum, na.rm=TRUE)

indval_examples <- filter_taxa(indval_examples, function(x){sum(x > 0) > 20}, prune = TRUE)                  
indval_examples.melt<-psmelt(indval_examples)
indval_examples_plot<-ggplot(indval_examples.melt, aes(x = week_str, y = log(Abundance),fill=season)) + 
  geom_boxplot(stat = "boxplot") + 
  facet_wrap("Genus", ncol=2) + 
  theme_ro()+
  theme(strip.text = element_text(face = 'italic'),legend.position = "none") +
  xlab("Weeks, January - December (2002-2018)") +
  ylab("Log10 Abundance (reads)") +theme(aspect.ratio = 0.5)+ 
  scale_fill_manual(values=c("#ED553B","#3CAEA3","#F6D55C","#20639B")) + 
  ggtitle("examples")+ ylim(0, 11.8)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
