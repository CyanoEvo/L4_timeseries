fig_3_taxa <- subset_taxa(l4_fungi,Genus == "Metschnikowia" |
                               Genus == "Epicoccum" |
                               Genus == "Cladosporium" |
                               Genus == "Symmetrospora" | 
                               Genus == "Rhodotorula" |
                               Genus == "Penicillium")

fig_3_taxa.melt<-psmelt(fig_3_taxa)
fig_3_taxa_plot<-ggplot(fig_3_taxa.melt, aes(x = week_str, y = log(Abundance),fill=Phylum)) + 
  geom_boxplot(stat = "boxplot") + 
  facet_wrap("Genus", ncol=1) + 
  theme_ro()+
  theme(
  #      axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank(),
        strip.text = element_text(face = 'italic'),legend.position = "none") +
  xlab("Weeks, January - December (2002-2018)") +
  ylab("Log10 Abundance (reads)") +theme(aspect.ratio = 0.5)+ 
  #geom_smooth(data = fig_3_taxa.melt, method = "gam", aes(week, log(Abundance), group =1), colour="red") +
  #  geom_vline(xintercept = 18, linetype = "longdash") +
  #  geom_vline(xintercept = 22, linetype = "longdash")+
  scale_fill_manual(values=c("white","white")) + 
  ggtitle("")+ ylim(0, 11.8)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))