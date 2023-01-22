#https://www.rdocumentation.org/packages/RVAideMemoire/versions/0.9-80/topics/pairwise.perm.manova

library(tidyr)
library(vegan)
library(reshape2)
library(phyloseq)
library(devtools)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mixOmics")
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(RVAideMemoire)
library(pairwiseAdonis)
library(pheatmap)
set.seed(12345)

l4_fungip <- l4_fungi
sl4_fungip <- prune_samples(sample_sums(l4_fungip) > 1000, l4_fungip)

l4_fungi_rare<- rarefy_even_depth(l4_fungip, sample.size = min(sample_sums(l4_fungip)),
                                  rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)


#are seasons different from each other?
ad_metadata <- as(sample_data(l4_fungi_rare), "data.frame")
dist<-distance(l4_fungi_rare, method="bray", permutations = 999)
adonis(dist ~ season,
       data = ad_metadata)
pairwise.adonis(dist,ad_metadata$season,p.adjust.m = "none")
#here I manually constructed a matrix season_pw from the output of pairwise.adonis
season_pw_mat<-as.matrix(season_pw)

pheatmap(season_pw_mat,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = TRUE, 
         color = c("#d53e4f","#f7f7f7"),
         breaks = c(0, 0.05, 1),  # distances 0 to 3 are red, 3 to 9 black
         main = 'Pairwise adonis (p-value)')

adonis(dist ~ month_str,
       data = ad_metadata)
pairwise.adonis(dist,ad_metadata$month_str,p.adjust.m = "none")
#here I manually constructed a matrix month_pw from the output of pairwise.adonis
month_pw_mat<-as.matrix(month_pw)


pheatmap(month_pw_mat,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = TRUE, 
         color = c("#d53e4f","#f7f7f7"),
         breaks = c(0, 0.05, 1),  # distances 0 to 3 are red, 3 to 9 black
         main = 'Pairwise adonis (p-value)')


adonis(dist ~ week_str,
       data = ad_metadata)
pairwise.adonis(dist,ad_metadata$week_str,p.adjust.m = "none")
#here I manually constructed a matrix pairwise_adonis from the output of pairwise.adonis
pwad<-as.matrix(pairwise_adonis)

pheatmap(pwad,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = TRUE, 
         color = c("#d53e4f","#f7f7f7"),
         breaks = c(0, 0.05, 1),  # distances 0 to 3 are red, 3 to 9 black
         main = 'Pairwise adonis (p-value)')


# alpha diversity ----------------------------------------------------------------
## metric calculation
l4_alpha <- estimate_richness(l4_fungi_rare, measures = c("Observed", "Chao1", "Shannon")) 
alpha_tib <- as_tibble(cbind(l4_alpha, sample_data(l4_fungi_rare)))
alpha_tib$Pielou <- alpha_tib$Shannon/log(alpha_tib$Observed)

## weekly plots - if arrival of new data makes trends clear at weekly scale, this is preferable to monthly
obs_plot<-ggplot(alpha_tib, aes(x = week, y = Observed, fill = as.factor(season))) +
  #  geom_boxplot()+
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.25, width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", shape = 21, alpha = 1, stroke = 0.3,  size = 3) +
  #  geom_smooth(data = alpha_tib, method = "gam", aes(week, Observed, group =1), colour="red") +
  xlab("") + 
  theme_ro() +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("#ED553B","#3CAEA3","#F6D55C","#20639B"))+ theme(aspect.ratio=0.18)+ 
  scale_x_continuous(breaks = seq(0, 52, by = 4))

shannon_plot<-ggplot(alpha_tib, aes(x = week, y = Shannon, fill = as.factor(season))) +
  #  geom_boxplot()+
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.25, width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", shape = 21, alpha = 1, stroke = 0.3,  size = 3) +
  #  geom_smooth(data = alpha_tib, method = "gam", aes(week, Observed, group =1), colour="red") +
  xlab("") + 
  theme_ro() +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("#ED553B","#3CAEA3","#F6D55C","#20639B"))+ theme(aspect.ratio=0.18)+ 
  scale_x_continuous(breaks = seq(0, 52, by = 4))

chao1_plot<-ggplot(alpha_tib, aes(x = week, y = Chao1, fill = as.factor(season))) +
  #  geom_boxplot()+
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.25, width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", shape = 21, alpha = 1, stroke = 0.3,  size = 3) +
  #  geom_smooth(data = alpha_tib, method = "gam", aes(week, Observed, group =1), colour="red") +
  xlab("") + 
  theme_ro() +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("#ED553B","#3CAEA3","#F6D55C","#20639B"))+ theme(aspect.ratio=0.18)+ 
  scale_x_continuous(breaks = seq(0, 52, by = 4))


pie_plot <- ggplot(alpha_tib, aes(x = week, y = Pielou, fill = as.factor(season))) +
  #  geom_boxplot()+
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.25, width = 0.1) +
  stat_summary(fun.y = mean, geom = "point", shape = 21, alpha = 1, stroke = 0.3,  size = 3) +
 #    geom_smooth(data = alpha_tib, method = "gam", aes(week, Pielou, group =1), colour="red") +
  xlab("") +
  theme_ro()+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("#ED553B","#3CAEA3","#F6D55C","#20639B"))+ theme(aspect.ratio=0.18)+ 
  scale_x_continuous(breaks = seq(0, 52, by = 4))

plot_grid(obs_plot,shannon_plot,chao1_plot,pie_plot,ncol=1,align = "vh")
kruskal.test(Chao1 ~ season, data = alpha_tib)
pairwise.wilcox.test(alpha_tib$Chao1, alpha_tib$season,
                     p.adjust.method = "none")

# beta diversity ----------------------------------------------------------
l4_otu <- as_tibble(as.data.frame(otu_table(l4_fungi_rare)))
rownames(l4_otu) <- rownames(otu_table(l4_fungi_rare))
l4_bray <- vegdist(l4_otu, method = "bray")

# CCA ---------------------------------------------------------------------

## remove variables with NA values (i.e. missing CTD data) this is necessary for CCA to run. Note that the corresponding rows in the otu_table will also have to be removed
l4_otu_cca <- as.data.frame(l4_otu)
sd_cca <-as(sample_data(l4_fungi_rare), "data.frame") %>% select(Tv290C, Sal00, total_gdf, season, mean_nirate_nitrite, mean_phosphate, mean_silicate,month_str)
rownames(l4_otu_cca) <- rownames(sd_cca)
setdiff(rownames(sd_cca), rownames(sd_cca[complete.cases(sd_cca), ])) # rownames which will also need to be removed from the otu_table note that many sites are lost at this stage.
l4_otu_cca <- l4_otu_cca[!rownames(l4_otu_cca) %in% c(setdiff(rownames(sd_cca), rownames(sd_cca[complete.cases(sd_cca), ]))),]
sd_cca <- sd_cca[complete.cases(sd_cca), ]

## z score transform environmental variables
sd_cca$Tv290C <- (sd_cca$Tv290C - mean(sd_cca$Tv290C))/sd(sd_cca$Tv290C)
sd_cca$Sal00 <- (sd_cca$Sal00 - mean(sd_cca$Sal00))/sd(sd_cca$Sal00)
sd_cca$total_gdf <- (sd_cca$total_gdf - mean(sd_cca$total_gdf))/sd(sd_cca$total_gdf)
sd_cca$mean_nirate_nitrite <- (sd_cca$mean_nirate_nitrite - mean(sd_cca$mean_nirate_nitrite))/sd(sd_cca$mean_nirate_nitrite)
sd_cca$mean_phosphate <- (sd_cca$mean_phosphate - mean(sd_cca$mean_phosphate))/sd(sd_cca$mean_phosphate)
sd_cca$mean_silicate <- (sd_cca$mean_silicate - mean(sd_cca$mean_silicate))/sd(sd_cca$mean_silicate)

# testing the cca model for variance inflation (if VIF over 10, highest scoring variable removed and test conducted again)
test1 <- cca(l4_otu_cca ~ Tv290C + Sal00 +  total_gdf + mean_nirate_nitrite + mean_phosphate + mean_silicate, data = sd_cca, scaling = "sites")
vif.cca(test1)
test1 <- cca(l4_otu_cca ~ Tv290C + Sal00 + total_gdf  + mean_nirate_nitrite + mean_phosphate + mean_silicate, data = sd_cca, scaling = "sites")
vif.cca(test1)
anova(test1, by = "terms", scaling = 1) # tests the significance of each term in teh model
summary(test1, scaling = 1) # temperature, salinity, and nitrate_nitrite all contribute significantly to constraining variance. Note that only 8% of variance is constrains (v. low) due to many zeros problem.
plot(test1, scaling = 1)

# plotting cca 
smry1 <- summary(test1, scaling = 1)
scrs1 <- as.data.frame(scores(test1, scaling = 1)$sites)
scrs1 <- cbind(scrs1, sd_cca)
biplot1 <- as.data.frame(smry1$biplot[c(1,3,4),]) # only select signficant environmental parameters for plotting

l4_cca <- ggplot(scrs1, aes(x=CCA1, y=CCA2)) + 
  geom_point(size = 3, shape = 21, stroke = 0.2, aes(fill = as.factor(month_str))) +
  geom_segment(data=biplot1, aes(x=0, xend=CCA1, y=0, yend=CCA2), 
               color="black", arrow=arrow(length=unit(0.01,"npc"))) +
  geom_label(data=biplot1, aes(x=CCA1,y=CCA2,label=rownames(biplot1),
                               hjust=0.5*(1-sign(CCA1)),vjust=0.5*(1-sign(CCA2))), 
             color="black", size=3) +
  theme_ro() +
  scale_fill_manual(values=c("#1d91c0","#225ea8",
                             "#78c679","#41ab5d","#238443",
                             "#ffffcc","#ffeda0","#fed976",
                             "#ec7014","#cc4c02","#993404",
                             "#41b6c4","red","green","blue","orange")) +
  theme(aspect.ratio=1)+ stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=season))


# beta diversity time-series -------------------------------------------------------------
## bind difference in weeks between samples with bray-curtis dissimilarity between samples 
sample_data(l4_fungi_rare)$week_diff  <- round(difftime(sample_data(l4_fungi_rare)$date, "2000-12-28", units = "weeks"))
l4_bray_ts <- as_tibble(cbind(as.vector(dist(sample_data(l4_fungi_rare)$week_diff, upper = F)),as.vector(l4_bray)))
colnames(l4_bray_ts) <- c("week_dist", "bray")

## time-series plot
l4_bray_ts_plot <- ggplot(l4_bray_ts, aes(x = week_dist, y = bray)) + 
  #  geom_rect(aes(xmin = 441,xmax = 495),fill="grey80") +
  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.25, width = 0.1, colour = "grey40") +
  stat_summary(fun.y = mean, geom = "point", shape = 21, alpha = 1, stroke = 0.3,  size = 1.5, fill = "white", colour = "grey40") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 1+(20*2)), colour = "firebrick1", fill = NA, size = 0.6) +
  geom_smooth(method = "lm", colour = "#4E73AB", fill = NA, size = 0.6) +
  theme_ro() + 
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0, 52*17, 52)) +
  scale_y_reverse() +
  labs(y = "Bray-Curtis Dissimilarity", x = "∆ Weeks")+
  geom_hline(yintercept = 0.92, colour = "grey", size = 0.3,alpha=0.5) +
  theme(legend.position = "none")+ theme(aspect.ratio=0.2)
