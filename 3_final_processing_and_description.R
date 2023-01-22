library(phyloseq)
library(ggplot2)

# To get to this point:
# - PS object made after main dada2 pipeline
# - PS object seperatated out into separate files
# - taxonomy table manually curated and saved as .csv
# - metadata manually curated 
# - additional metadata added using metadata analysis script and saved

# Load curated data ------------------------------------------------------------
setwd("~/Projects/L4_timeseries/L4_Nathan/tidy_analysis")
l4_TAX<-read.csv("l4_TAX.csv",header = T, row.names=1)
l4_OTU<-read.csv("l4_OTU_1.csv",header = T, row.names=1)
l4_META_sd<-read.csv("L4_SD.csv",header = T, row.names=22)

# MAKE SURE SAMPLE NAMES MATCH FOR META AND OTU
a1<-row.names(l4_OTU)
a2<-row.names(l4_META_sd)
setdiff(a1,a2)

L4_ps<- phyloseq(tax_table(as.matrix(l4_TAX)), 
                 otu_table(l4_OTU, taxa_are_rows = FALSE),
                 sample_data(l4_META_sd))

#Remove incomplete years
L4_ps<-subset_samples(L4_ps, !year=="2019")
L4_ps<-subset_samples(L4_ps, !year=="2001")

# basic dataset descriptiojn
L4_ps.melt1<-psmelt(L4_ps)
L4_ps.totals <- aggregate(L4_ps.melt1$Abundance,
                          by=list(L4_ps.melt1$Kingdom,
                                  L4_ps.melt1$Phylum,
                                  L4_ps.melt1$OTU), 
                          FUN=sum, na.rm=TRUE)
kingdom_totals<-data.frame(table(L4_ps.totals$Group.1))

kingdom_plot<-ggplot(kingdom_totals, aes(x=Var1, y=Freq, fill= Var1)) + 
  geom_bar(stat="identity", position = "stack", colour="black") +
  scale_y_continuous(limits = c(0, 1500), labels = scales::comma, breaks = seq(0, 1500, by = 250)) +
  xlab("") + ylab("n ASVs") +
  theme_ro() + theme(aspect.ratio=1.5) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = c(.5, .5),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  scale_fill_manual(values=c("#ed553b","white","grey"))+ labs(fill='Kingdom') 

#Focus on the fungi
l4_all_fungi <- subset_taxa(L4_ps, Kingdom == "Fungi")
l4_all_fungi.melt<-psmelt(l4_all_fungi)
l4_all_fungi.totals <- aggregate(l4_all_fungi.melt$Abundance,
                                 by=list(l4_all_fungi.melt$Kingdom,
                                         l4_all_fungi.melt$Phylum,
                                         l4_all_fungi.melt$OTU),
                                 FUN=sum, na.rm=TRUE)
fungi_totals<-data.frame(table(l4_all_fungi.totals$Group.2))

fungi_plot<-ggplot(fungi_totals, aes(x=Var1, y=Freq, fill= Var1)) + 
  geom_bar(stat="identity", position = "stack", colour="black") +
  scale_y_continuous(limits = c(0, 662), labels = scales::comma, breaks = seq(0, 700, by = 100)) +
  xlab("") + ylab("n ASVs") +
  theme_ro() + theme(aspect.ratio=1.5) +
  theme(axis.text.x = element_blank()) + 
  theme(legend.position = c(.5, .5),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) +
  scale_fill_manual(values=c("#173f5f","#20639b","#3caea3","#f6d55c","#ed553b","grey"))+ labs(fill='Phylum')

plot_grid(kingdom_plot,fungi_plot, ncol=2, align ="hv")

# For the diversity analyses we want to use rarefied counts
# We only want to use real fungi for this
l4_fungi <- subset_taxa(l4_all_fungi, !Phylum == "Unknown")
l4_fungi.melt<-psmelt(l4_fungi)

rarecurve(t(otu_table(l4_fungi)), step=50, cex=0.5)
# Take out small samples or rarefaction will be to extreme

#rarefied to 1008 reads
#61samples lost
#21 OTUS lost