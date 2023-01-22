# recurrence index --------------------------------------------------------
# See Recurrence Index repository https://github.com/CaterinaRG/Recurrence-Index 
# Function:
seasonality.test<-function(comm.tab,n=1000,probs=c(0.025, 0.975),lag.max=884,na.action=na.pass) # lag max should be total no. weeks (i.e. 905)
{
  require(vegan)
  
  season.index<-function(x){
    acf.all<-apply(as.matrix(x),2,acf,plot=F,lag.max=lag.max,na.action=na.pass)
    acf.all<-sapply(acf.all,"[[",1)
    apply(acf.all,2,function(x) sum(abs(x)))
  }
  
  n<-n
  season.index.real<-season.index(comm.tab)
  
  names(season.index.real) <- colnames(comm.tab)
  season.index.simul<-matrix(NA, ncol = dim(comm.tab)[2],nrow = n)
  for (i in 1:n) {
    season.index.simul[i, ]<-season.index(comm.tab[sample(1:nrow(comm.tab)),])
  }
  colnames(season.index.simul) <- colnames(comm.tab)
  season.index.simul <- as.data.frame(season.index.simul)
  media <- apply(season.index.simul, 2, mean, na.rm=TRUE) #na.rm pq funcioni amb NA
  ci <- apply(season.index.simul, 2, quantile, probs = probs, na.rm=TRUE) #na.rm pq funcioni amb NA
  resultats <- data.frame(observed = season.index.real, mean.simulated = media,lowCI = ci[1, ], uppCI = ci[2, ], sign = NA)
  for (j in 1:dim(resultats)[1]) {
    if (resultats$observed[j] > resultats$uppCI[j]) 
      resultats$sign[j] <- "SIG. HIGHER"
    if (resultats$observed[j] < resultats$lowCI[j]) 
      resultats$sign[j] <- "SIG. LOWER"
    if (resultats$observed[j] >= resultats$lowCI[j] & resultats$observed[j] <= 
        resultats$uppCI[j]) 
      resultats$sign[j] <- "NON SIG."
  }
  resultats$sign <- as.factor(resultats$sign)
  resultats
}

l4_fungi = prune_taxa(taxa_sums(l4_fungi) > 0, l4_fungi) 
####################Exophiala####################
Exophiala_stats <- subset_taxa(l4_fungi,Genus == "Exophiala")
Exophiala_stats.melt<-psmelt(Exophiala_stats)
Exophiala_stats_otu<-data.frame(otu_table(Exophiala_stats))
Exophiala_stats.melt$week_diff <- round(difftime(Exophiala_stats.melt$date, "2002-01-01", units = "weeks"))
Exophiala_stats.melt$week_diff <- as.numeric(Exophiala_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Exophiala_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Exophiala_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Exophiala_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Exophiala_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                lag.max = 884, na.action = na.pass)$lag, 
                            acf(l4_otu_dummy %>% select(totals), 
                                lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Exophiala_acf_plot <- ggplot(Exophiala_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Exophiala_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Exophiala####################

####################Cadophora####################
Cadophora_stats <- subset_taxa(l4_fungi,Genus == "Cadophora")
Cadophora_stats.melt<-psmelt(Cadophora_stats)
Cadophora_stats_otu<-data.frame(otu_table(Cadophora_stats))
Cadophora_stats.melt$week_diff <- round(difftime(Cadophora_stats.melt$date, "2002-01-01", units = "weeks"))
Cadophora_stats.melt$week_diff <- as.numeric(Cadophora_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Cadophora_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Cadophora_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Cadophora_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Cadophora_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                   lag.max = 884, na.action = na.pass)$lag, 
                               acf(l4_otu_dummy %>% select(totals), 
                                   lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Cadophora_acf_plot <- ggplot(Cadophora_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_hline(yintercept = 0.08, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Cadophora_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Cadophora####################

####################Metschnikowia####################
Metschnikowia_stats <- subset_taxa(l4_fungi,Genus == "Metschnikowia")
Metschnikowia_stats.melt<-psmelt(Metschnikowia_stats)
Metschnikowia_stats_otu<-data.frame(otu_table(Metschnikowia_stats))
Metschnikowia_stats.melt$week_diff <- round(difftime(Metschnikowia_stats.melt$date, "2002-01-01", units = "weeks"))
Metschnikowia_stats.melt$week_diff <- as.numeric(Metschnikowia_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Metschnikowia_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Metschnikowia_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Metschnikowia_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Metschnikowia_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                    lag.max = 884, na.action = na.pass)$lag, 
                acf(l4_otu_dummy %>% select(totals), 
                    lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Metschnikowia_acf_plot <- ggplot(Metschnikowia_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_hline(yintercept = 0.05, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Metschnikowia_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Metschnikowia####################

####################Cladosporium####################
Cladosporium_stats <- subset_taxa(l4_fungi,Genus == "Cladosporium")
Cladosporium_stats.melt<-psmelt(Cladosporium_stats)
Cladosporium_stats_otu<-data.frame(otu_table(Cladosporium_stats))
Cladosporium_stats.melt$week_diff <- round(difftime(Cladosporium_stats.melt$date, "2002-01-01", units = "weeks"))
Cladosporium_stats.melt$week_diff <- as.numeric(Cladosporium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Cladosporium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Cladosporium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Cladosporium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Cladosporium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$lag, 
                                   acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Cladosporium_acf_plot <- ggplot(Cladosporium_acf, aes(x = V1, y = V2)) + 
 # geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_hline(yintercept = 0.05, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Cladosporium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Cladosporium####################

####################Rhodotorula####################
Rhodotorula_stats <- subset_taxa(l4_fungi,Genus == "Rhodotorula")
Rhodotorula_stats.melt<-psmelt(Rhodotorula_stats)
Rhodotorula_stats_otu<-data.frame(otu_table(Rhodotorula_stats))
Rhodotorula_stats.melt$week_diff <- round(difftime(Rhodotorula_stats.melt$date, "2002-01-01", units = "weeks"))
Rhodotorula_stats.melt$week_diff <- as.numeric(Rhodotorula_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Rhodotorula_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Rhodotorula_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Rhodotorula_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Rhodotorula_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                      lag.max = 884, na.action = na.pass)$lag, 
                                  acf(l4_otu_dummy %>% select(totals), 
                                      lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Rhodotorula_acf_plot <- ggplot(Rhodotorula_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Rhodotorula_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Rhodotorula####################

####################Parengyodontium####################
Parengyodontium_stats <- subset_taxa(l4_fungi,Genus == "Parengyodontium")
Parengyodontium_stats.melt<-psmelt(Parengyodontium_stats)
Parengyodontium_stats<-subset_taxa(l4_fungi, Genus == "Parengyodontium")
Parengyodontium_stats = prune_taxa(c("asv_19","asv_41","asv_2904"), Parengyodontium_stats)
Parengyodontium_stats_otu<-data.frame(otu_table(Parengyodontium_stats))
Parengyodontium_stats.melt$week_diff <- round(difftime(Parengyodontium_stats.melt$date, "2002-01-01", units = "weeks"))
Parengyodontium_stats.melt$week_diff <- as.numeric(Parengyodontium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Parengyodontium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Parengyodontium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Parengyodontium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Parengyodontium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$lag, 
                                 acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Parengyodontium_acf_plot <- ggplot(Parengyodontium_acf, aes(x = V1, y = V2)) + 
 #geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Parengyodontium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Parengyodontium####################

####################Engyodontium####################
Engyodontium_stats <- subset_taxa(l4_fungi,Genus == "Engyodontium")
Engyodontium_stats.melt<-psmelt(Engyodontium_stats)
Engyodontium_stats_otu<-data.frame(otu_table(Engyodontium_stats))
Engyodontium_stats.melt$week_diff <- round(difftime(Engyodontium_stats.melt$date, "2002-01-01", units = "weeks"))
Engyodontium_stats.melt$week_diff <- as.numeric(Engyodontium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Engyodontium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Parengyodontium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Engyodontium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Engyodontium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                         lag.max = 884, na.action = na.pass)$lag, 
                                     acf(l4_otu_dummy %>% select(totals), 
                                         lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Engyodontium_acf_plot <- ggplot(Engyodontium_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Engyodontium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Engyodontium####################

####################Penicillium####################
Penicillium_stats <- subset_taxa(l4_fungi,Genus == "Penicillium")
Penicillium_stats.melt<-psmelt(Penicillium_stats)
Penicillium_stats_otu<-data.frame(otu_table(Penicillium_stats))
Penicillium_stats.melt$week_diff <- round(difftime(Penicillium_stats.melt$date, "2002-01-01", units = "weeks"))
Penicillium_stats.melt$week_diff <- as.numeric(Penicillium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Penicillium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Penicillium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Penicillium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Penicillium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                      lag.max = 884, na.action = na.pass)$lag, 
                                  acf(l4_otu_dummy %>% select(totals), 
                                      lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Penicillium_acf_plot <- ggplot(Penicillium_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Penicillium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Penicillium####################

####################Sarocladium####################
Sarocladium_stats <- subset_taxa(l4_fungi,Genus == "Sarocladium")
Sarocladium_stats.melt<-psmelt(Sarocladium_stats)
Sarocladium_stats_otu<-data.frame(otu_table(Sarocladium_stats))
Sarocladium_stats.melt$week_diff <- round(difftime(Sarocladium_stats.melt$date, "2002-01-01", units = "weeks"))
Sarocladium_stats.melt$week_diff <- as.numeric(Sarocladium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Sarocladium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Sarocladium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Sarocladium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Sarocladium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$lag, 
                                 acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Sarocladium_acf_plot <- ggplot(Sarocladium_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Sarocladium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Sarocladium####################

####################Debaryomyces####################
Debaryomyces_stats <- subset_taxa(l4_fungi,Genus == "Debaryomyces")
Debaryomyces_stats.melt<-psmelt(Debaryomyces_stats)
Debaryomyces_stats = prune_taxa(c("asv_29","asv_1131","asv_4346"), Debaryomyces_stats)
Debaryomyces_stats_otu<-data.frame(otu_table(Debaryomyces_stats))
Debaryomyces_stats.melt$week_diff <- round(difftime(Debaryomyces_stats.melt$date, "2002-01-01", units = "weeks"))
Debaryomyces_stats.melt$week_diff <- as.numeric(Debaryomyces_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Debaryomyces_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Debaryomyces_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Debaryomyces_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Debaryomyces_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$lag, 
                                 acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Debaryomyces_acf_plot <- ggplot(Debaryomyces_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Debaryomyces_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Debaryomyces####################

####################Kondoa####################
Kondoa_stats <- subset_taxa(l4_fungi,Genus == "Kondoa")
Kondoa_stats.melt<-psmelt(Kondoa_stats)
Kondoa_stats_otu<-data.frame(otu_table(Kondoa_stats))
Kondoa_stats.melt$week_diff <- round(difftime(Kondoa_stats.melt$date, "2002-01-01", units = "weeks"))
Kondoa_stats.melt$week_diff <- as.numeric(Kondoa_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Kondoa_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Kondoa_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Kondoa_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Kondoa_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                lag.max = 884, na.action = na.pass)$lag, 
                            acf(l4_otu_dummy %>% select(totals), 
                                lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Kondoa_acf_plot <- ggplot(Kondoa_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Kondoa_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))+
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Kondoa####################

####################Epcicoccum####################
Epicoccum_stats <- subset_taxa(l4_fungi,Genus == "Epicoccum")
Epicoccum_stats.melt<-psmelt(Epicoccum_stats)
Epicoccum_stats_otu<-data.frame(otu_table(Epicoccum_stats))
Epicoccum_stats.melt$week_diff <- round(difftime(Epicoccum_stats.melt$date, "2002-01-01", units = "weeks"))
Epicoccum_stats.melt$week_diff <- as.numeric(Epicoccum_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Epicoccum_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Epicoccum_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Epicoccum_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Epicoccum_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                      lag.max = 884, na.action = na.pass)$lag, 
                                  acf(l4_otu_dummy %>% select(totals), 
                                      lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Epicoccum_acf_plot <- ggplot(Epicoccum_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Epcicoccum_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Epcicoccum####################

####################Filobasidium####################
Filobasidium_stats <- subset_taxa(l4_fungi,Genus == "Filobasidium")
Filobasidium_stats.melt<-psmelt(Filobasidium_stats)
Filobasidium_stats_otu<-data.frame(otu_table(Filobasidium_stats))
Filobasidium_stats.melt$week_diff <- round(difftime(Filobasidium_stats.melt$date, "2002-01-01", units = "weeks"))
Filobasidium_stats.melt$week_diff <- as.numeric(Filobasidium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Filobasidium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Filobasidium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Filobasidium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Filobasidium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                   lag.max = 884, na.action = na.pass)$lag, 
                               acf(l4_otu_dummy %>% select(totals), 
                                   lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Filobasidium_acf_plot <- ggplot(Filobasidium_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("ACF") + 
  ggtitle("", subtitle = "Filobasidium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Filobasidium####################

####################Sporidiobolus####################
Sporidiobolus_stats <- subset_taxa(l4_fungi,Genus == "Sporidiobolus")
Sporidiobolus_stats.melt<-psmelt(Sporidiobolus_stats)
Sporidiobolus_stats_otu<-data.frame(otu_table(Sporidiobolus_stats))
Sporidiobolus_stats.melt$week_diff <- round(difftime(Sporidiobolus_stats.melt$date, "2002-01-01", units = "weeks"))
Sporidiobolus_stats.melt$week_diff <- as.numeric(Sporidiobolus_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Sporidiobolus_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Sporidiobolus_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Sporidiobolus_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Sporidiobolus_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                  lag.max = 884, na.action = na.pass)$lag, 
                              acf(l4_otu_dummy %>% select(totals), 
                                  lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Sporidiobolus_acf_plot <- ggplot(Sporidiobolus_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Sporidiobolus_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Sporidiobolus####################

####################Symmetrospora####################
Symmetrospora_stats <- subset_taxa(l4_fungi,Genus == "Symmetrospora")
Symmetrospora_stats.melt<-psmelt(Symmetrospora_stats)
Symmetrospora_stats_otu<-data.frame(otu_table(Symmetrospora_stats))
Symmetrospora_stats.melt$week_diff <- round(difftime(Symmetrospora_stats.melt$date, "2002-01-01", units = "weeks"))
Symmetrospora_stats.melt$week_diff <- as.numeric(Symmetrospora_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Symmetrospora_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Symmetrospora_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Symmetrospora_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Symmetrospora_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$lag, 
                                   acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Symmetrospora_acf_plot <- ggplot(Symmetrospora_acf, aes(x = V1, y = V2)) + 
 # geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Symmetrospora_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Symmetrospora####################

####################Simplicillium####################
Simplicillium_stats <- subset_taxa(l4_fungi,Genus == "Simplicillium")
Simplicillium_stats.melt<-psmelt(Simplicillium_stats)
Simplicillium_stats_otu<-data.frame(otu_table(Simplicillium_stats))
Simplicillium_stats.melt$week_diff <- round(difftime(Simplicillium_stats.melt$date, "2002-01-01", units = "weeks"))
Simplicillium_stats.melt$week_diff <- as.numeric(Simplicillium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Simplicillium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Simplicillium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Simplicillium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Simplicillium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$lag, 
                                   acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Simplicillium_acf_plot <- ggplot(Simplicillium_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Simplicillium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Simplicillium####################

####################Leucosporidium####################
Leucosporidium_stats <- subset_taxa(l4_fungi,Genus == "Leucosporidium")
Leucosporidium_stats.melt<-psmelt(Leucosporidium_stats)
Leucosporidium_stats_otu<-data.frame(otu_table(Leucosporidium_stats))
Leucosporidium_stats.melt$week_diff <- round(difftime(Leucosporidium_stats.melt$date, "2002-01-01", units = "weeks"))
Leucosporidium_stats.melt$week_diff <- as.numeric(Leucosporidium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Leucosporidium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Leucosporidium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Leucosporidium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Leucosporidium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$lag, 
                                   acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Leucosporidium_acf_plot <- ggplot(Leucosporidium_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Leucosporidium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Leucosporidium####################

####################Aureobasidium####################
Aureobasidium_stats <- subset_taxa(l4_fungi,Genus == "Aureobasidium")
Aureobasidium_stats.melt<-psmelt(Aureobasidium_stats)
Aureobasidium_stats_otu<-data.frame(otu_table(Aureobasidium_stats))
Aureobasidium_stats.melt$week_diff <- round(difftime(Aureobasidium_stats.melt$date, "2002-01-01", units = "weeks"))
Aureobasidium_stats.melt$week_diff <- as.numeric(Aureobasidium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Aureobasidium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Aureobasidium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Aureobasidium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Aureobasidium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                        lag.max = 884, na.action = na.pass)$lag, 
                                    acf(l4_otu_dummy %>% select(totals), 
                                        lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Aureobasidium_acf_plot <- ggplot(Aureobasidium_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Aureobasidium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Aureobasidium####################


######################END############################


####################Aspergillus####################
Aspergillus_stats_otu<-data.frame(otu_table(Aspergillus_stats))
Aspergillus_stats.melt$week_diff <- round(difftime(Aspergillus_stats.melt$date, "2002-01-01", units = "weeks"))
Aspergillus_stats.melt$week_diff <- as.numeric(Aspergillus_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Aspergillus_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Aspergillus_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Aspergillus_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Aspergillus_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$lag, 
                                   acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Aspergillus_acf_plot <- ggplot(Aspergillus_acf, aes(x = V1, y = V2)) + 
 # geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Aspergillus_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Aspergillus####################

####################Wallemia####################
Wallemia_stats_otu<-data.frame(otu_table(Wallemia_stats))
Wallemia_stats.melt$week_diff <- round(difftime(Wallemia_stats.melt$date, "2002-01-01", units = "weeks"))
Wallemia_stats.melt$week_diff <- as.numeric(Wallemia_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Wallemia_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Wallemia_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Wallemia_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=1, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Wallemia_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$lag, 
                                 acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Wallemia_acf_plot <- ggplot(Wallemia_acf, aes(x = V1, y = V2)) + 
 # geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Wallemia_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Wallemia####################

####################Cystobasidium####################
Cystobasidium_stats_otu<-data.frame(otu_table(Cystobasidium_stats))
Cystobasidium_stats.melt$week_diff <- round(difftime(Cystobasidium_stats.melt$date, "2002-01-01", units = "weeks"))
Cystobasidium_stats.melt$week_diff <- as.numeric(Cystobasidium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Cystobasidium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Cystobasidium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Cystobasidium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Cystobasidium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                  lag.max = 884, na.action = na.pass)$lag, 
                              acf(l4_otu_dummy %>% select(totals), 
                                  lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Cystobasidium_acf_plot <- ggplot(Cystobasidium_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Cystobasidium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Cystobasidium####################

####################Acremonium####################
Acremonium_stats_otu<-data.frame(otu_table(Acremonium_stats))
Acremonium_stats.melt$week_diff <- round(difftime(Acremonium_stats.melt$date, "2002-01-01", units = "weeks"))
Acremonium_stats.melt$week_diff <- as.numeric(Acremonium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Acremonium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Acremonium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Acremonium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Acremonium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$lag, 
                                   acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Acremonium_acf_plot <- ggplot(Acremonium_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Acremonium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Acremonium####################

####################Neoascochyta####################
Neoascochyta_stats_otu<-data.frame(otu_table(Neoascochyta_stats))
Neoascochyta_stats.melt$week_diff <- round(difftime(Neoascochyta_stats.melt$date, "2002-01-01", units = "weeks"))
Neoascochyta_stats.melt$week_diff <- as.numeric(Neoascochyta_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Neoascochyta_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Neoascochyta_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Neoascochyta_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Neoascochyta_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                    lag.max = 884, na.action = na.pass)$lag, 
                                acf(l4_otu_dummy %>% select(totals), 
                                    lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Neoascochyta_acf_plot <- ggplot(Neoascochyta_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Neoascochyta_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Neoascochyta####################

####################Peniophora####################
Peniophora_stats_otu<-data.frame(otu_table(Peniophora_stats))
Peniophora_stats.melt$week_diff <- round(difftime(Peniophora_stats.melt$date, "2002-01-01", units = "weeks"))
Peniophora_stats.melt$week_diff <- as.numeric(Peniophora_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Peniophora_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Peniophora_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Peniophora_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Peniophora_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                      lag.max = 884, na.action = na.pass)$lag, 
                                  acf(l4_otu_dummy %>% select(totals), 
                                      lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Peniophora_acf_plot <- ggplot(Peniophora_acf, aes(x = V1, y = V2)) + 
 # geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Peniophora_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Peniophora####################

####################Botrytis####################
Botrytis_stats_otu<-data.frame(otu_table(Botrytis_stats))
Botrytis_stats.melt$week_diff <- round(difftime(Botrytis_stats.melt$date, "2002-01-01", units = "weeks"))
Botrytis_stats.melt$week_diff <- as.numeric(Botrytis_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Botrytis_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Botrytis_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Botrytis_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Botrytis_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                    lag.max = 884, na.action = na.pass)$lag, 
                                acf(l4_otu_dummy %>% select(totals), 
                                    lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Botrytis_acf_plot <- ggplot(Botrytis_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Botrytis_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Botrytis####################

####################Trametes####################
Trametes_stats_otu<-data.frame(otu_table(Trametes_stats))
Trametes_stats.melt$week_diff <- round(difftime(Trametes_stats.melt$date, "2002-01-01", units = "weeks"))
Trametes_stats.melt$week_diff <- as.numeric(Trametes_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Trametes_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Trametes_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Trametes_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Trametes_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                      lag.max = 884, na.action = na.pass)$lag, 
                                  acf(l4_otu_dummy %>% select(totals), 
                                      lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Trametes_acf_plot <- ggplot(Trametes_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Trametes_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Trametes####################

####################Monocillium####################
Monocillium_stats_otu<-data.frame(otu_table(Monocillium_stats))
Monocillium_stats.melt$week_diff <- round(difftime(Monocillium_stats.melt$date, "2002-01-01", units = "weeks"))
Monocillium_stats.melt$week_diff <- as.numeric(Monocillium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Monocillium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Monocillium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Monocillium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Monocillium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$lag, 
                                 acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Monocillium_acf_plot <- ggplot(Monocillium_acf, aes(x = V1, y = V2)) + 
 # geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Monocillium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Monocillium####################

####################Vishniacozyma####################
Vishniacozyma_stats_otu<-data.frame(otu_table(Vishniacozyma_stats))
Vishniacozyma_stats.melt$week_diff <- round(difftime(Vishniacozyma_stats.melt$date, "2002-01-01", units = "weeks"))
Vishniacozyma_stats.melt$week_diff <- as.numeric(Vishniacozyma_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Vishniacozyma_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Vishniacozyma_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Vishniacozyma_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Vishniacozyma_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$lag, 
                                 acf(l4_otu_dummy %>% select(totals), 
                                     lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Vishniacozyma_acf_plot <- ggplot(Vishniacozyma_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("∆ Weeks") + 
  ylab("") + 
  ggtitle("", subtitle = "Vishniacozyma_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Vishniacozyma####################

####################Protomycesa####################
Protomyces_stats_otu<-data.frame(otu_table(Protomyces_stats))
Protomyces_stats.melt$week_diff <- round(difftime(Protomyces_stats.melt$date, "2002-01-01", units = "weeks"))
Protomyces_stats.melt$week_diff <- as.numeric(Protomyces_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Protomyces_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Protomyces_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Protomyces_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Protomyces_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$lag, 
                                   acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Protomyces_acf_plot <- ggplot(Protomyces_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Protomyces_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)
####################Protomyces####################

####################Erythrobasidium####################
Erythrobasidium_stats_otu<-data.frame(otu_table(Erythrobasidium_stats))
Erythrobasidium_stats.melt$week_diff <- round(difftime(Erythrobasidium_stats.melt$date, "2002-01-01", units = "weeks"))
Erythrobasidium_stats.melt$week_diff <- as.numeric(Erythrobasidium_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Erythrobasidium_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Erythrobasidium_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Erythrobasidium_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Erythrobasidium_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                    lag.max = 884, na.action = na.pass)$lag, 
                                acf(l4_otu_dummy %>% select(totals), 
                                    lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Erythrobasidium_acf_plot <- ggplot(Erythrobasidium_acf, aes(x = V1, y = V2)) + 
#  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Erythrobasidium_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 1))
####################Erythrobasidiums####################

####################Armillaria####################
Armillaria_stats<-subset_taxa(L4_fungi, Genus == "Armillaria")
Armillaria_stats_otu<-data.frame(otu_table(Armillaria_stats))
Armillaria_stats.melt<-psmelt(Armillaria_stats)
Armillaria_stats.melt$week_diff <- round(difftime(Armillaria_stats.melt$date, "2002-01-01", units = "weeks"))
Armillaria_stats.melt$week_diff <- as.numeric(Armillaria_stats.melt$week_diff)
## setup otu_table with dummy weeks
l4_otu_dummy <- Armillaria_stats_otu %>% mutate(week_diff = as.numeric(round(difftime(sample_data(Armillaria_stats)$date, "2002-01-01", units = "weeks"))))
dummy.weeks <- data.frame(week_diff = c(1:max(Armillaria_stats.melt$week_diff)))
## week_diff column must also be unselected or it will be treated as a species erroneously
l4_otu_dummy <- as_tibble(left_join(dummy.weeks, l4_otu_dummy, by = "week_diff")) %>% select(-c("week_diff"))
l4_otu_dummy$totals<-rowSums(l4_otu_dummy)
l4_otu_dummy[is.na(l4_otu_dummy)] <- 0
## the recurrence index calculation fails when the maximum lag is too high. This interacts with NA values, but currently cannot successfully be set beyond n = 843 (i.e. 843 months). 
SI_test <- seasonality.test(as.matrix(l4_otu_dummy),n=884, na.action = na.pass)
## tibble of 'recurrent' ASVs Check totals for whole taxon
SI_tib <- as_tibble(SI_test) %>% mutate(asv = rownames(SI_test)) %>% filter(sign == "SIGNIFICANTLY HIGHER")

Armillaria_acf<-as_tibble(cbind(acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$lag, 
                                   acf(l4_otu_dummy %>% select(totals), 
                                       lag.max = 884, na.action = na.pass)$acf))
## Plotting seasonal recurrence
Armillaria_acf_plot <- ggplot(Armillaria_acf, aes(x = V1, y = V2)) + 
  #  geom_vline(xintercept = seq(0, 52*17, 52), size = 0.2, colour = "black", linetype = 3) + 
  geom_hline(yintercept = 0, colour = "black", size = 0.2) +
  geom_bar(stat = "identity", colour = "black") + 
  xlab("") + 
  ylab("") + 
  ggtitle("", subtitle = "Armillaria_acf") +
  scale_x_continuous(limits = c(0,(52*17)), expand = c(0,0), breaks = seq(0,52*17, 52)) +
  theme_ro()   +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  ylim(-0.05,0.5)+theme(aspect.ratio = 0.5)+ theme(plot.subtitle = element_text(hjust = 0))
####################Armillaria####################

#9*7"
acf_point5<-plot_grid(Exophiala_acf_plot,
                      Kondoa_acf_plot,
                      Symmetrospora_acf_plot,
                      Epicoccum_acf_plot,
                      Sarocladium_acf_plot,
                      Metschnikowia_acf_plot,
                      Parengyodontium_acf_plot,
                      Cladosporium_acf_plot,
                      Simplicillium_acf_plot,
                      Engyodontium_acf_plot,
                      Filobasidium_acf_plot,
                      Rhodotorula_acf_plot,
                      Cadophora_acf_plot, 
                      Sporidiobolus_acf_plot,
                      Penicillium_acf_plot,
                      Debaryomyces_acf_plot,ncol=4)

acf_fig_3<-plot_grid(Metschnikowia_acf_plot, 
                     Epicoccum_acf_plot,
                      Cladosporium_acf_plot,
                     Symmetrospora_acf_plot,
                      Rhodotorula_acf_plot,
                      Penicillium_acf_plot,ncol=1)

labels=c("2019","2018","2017","2016","2015","2014","2013","2012","2011","2010","2009","2008","2007","2006","2005","2004","2003","2002")