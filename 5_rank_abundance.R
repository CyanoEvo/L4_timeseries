library(dplyr)

# for all genus based anlysis we futher rationalise the dataset removing all uncertain genera
l4_fungi_genus<-l4_all_fungi
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Unknown Ascomycota" )
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Unknown Basidiomycota")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Dothideomycetes incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Cystobasidiomycetes incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Helotiales incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Agaricomycetes incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Pleosporales incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Capnodiales incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Bionectriaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Unknown")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Agaricales incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Amylocorticiaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Cucurbitariaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Boletaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Chaetothyriales incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Dictyosporiaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Cystobasidiaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Didymosphaeriaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Erythrobasidiaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Glomeraceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Helotiaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Herpotrichiellaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Hypocreales incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Lasiosphaeriaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Lecanorales incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Lecanoromycetes incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Malasseziales incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Microbotryomycetes incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Mortierellales incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Peniophoraceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Phaeosphaeriaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Roccellaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Sordariomycetes incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Sporidiobolaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Sporormiaceae incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Trechisporales incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Tremellales incertae sedis")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Unknown Glomeromycota")
l4_fungi_genus<-subset_taxa(l4_fungi_genus, !Genus == "Unknown Mucoromycota")

# Abundance plots
# get the most abundant taxa
#make sure we are using the unrarefied data

l4_fungi_rank_abundance.melt<-psmelt(l4_fungi_genus)

l4_fungi_rank_abundance.melt1 <- l4_fungi_rank_abundance.melt[,c("Abundance", "Genus","Phylum")]
l4_fungi_rank_abundance.melt1 <- aggregate(l4_fungi_rank_abundance.melt1$Abundance,
                                           by=list(l4_fungi_rank_abundance.melt$Genus,l4_fungi_rank_abundance.melt$Phylum), FUN=sum, na.rm=TRUE)
colnames(l4_fungi_rank_abundance.melt1)<-c("Genus", "Phylum","Abundance")
l4_fungi_rank_abundance.melt1<-arrange(l4_fungi_rank_abundance.melt1, -Abundance)
l4_fungi_rank_abundance.melt1$Genus <- factor(l4_fungi_rank_abundance.melt1$Genus, levels = l4_fungi_rank_abundance.melt1$Genus[order(l4_fungi_rank_abundance.melt1$Abundance)])
#choose this wisely - here we use the TOTAL unfiltered dataset (all years)
Total_Abundance<-sum(l4_all_fungi.melt$Abundance)
l4_fungi_rank_abundance.melt1$Proportion_total<-l4_fungi_rank_abundance.melt1$Abundance/Total_Abundance*100
l4_fungi_toptaxa<-subset(l4_fungi_rank_abundance.melt1, Proportion_total>1)

all_plot<-ggplot(l4_fungi_rank_abundance.melt1,aes(x=Abundance, y=Genus, fill = Phylum)) +
  geom_bar(stat="identity", position = "stack", colour="black") +
  theme(legend.position = "none") +
  theme_ro() + 
  theme(legend.position = "none",
        axis.text.y=element_text(face = 'italic'))+
  xlab("Total reads 2002-2019") +
  ylab("All genera") +
  theme(aspect.ratio=2) 

mostabundant_plot<-ggplot(l4_fungi_toptaxa,aes(x=Abundance, y=Genus, fill = Phylum)) +
  geom_bar(stat="identity", position = "stack", colour="black") +
  theme(legend.position = "none") +
  theme_ro() + 
  theme(legend.position = "none",
        axis.text.y=element_text(face = 'italic'))+
  xlab("Total reads 2002-2019") +
  ylab("Genera comprising > 1 % Total fungal reads") +
  theme(aspect.ratio=2) +
#  scale_fill_manual(values=c("#20639b","#ed553b")) +
  theme(legend.position = c(.6, .5))
