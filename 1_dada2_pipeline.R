# Pipeline based on only forward ITS reads, according to recommendations of Pauvert et al. (2019; Fungal Ecology; method Si5)

# Dependencies ---------------------------------------------------
install.packages("devtools")
library(devtools)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2", version = "3.9")
BiocManager::install("phyloseq", version = "3.9")
BiocManager::install("ShortRead", version = "3.9")
BiocManager::install("Biostrings", version = "3.9")
library(dada2)
library(phyloseq)
library(ShortRead)
library(Biostrings)
library(seqinr)
library(dplyr)
library(tibbletime)
library(lubridate)

install.packages("ggplot2")
install.packages("cowplot")
library(ggplot2)
library(cowplot)


# File path ------------------------------------------------------------
## set path to fastq location
path <- ("~/DATA/L4_fungi/files")

## retain only forward reads
fnFs <- sort(list.files(path, pattern = "R1_001", full.names = TRUE))

# Primer removal ----------------------------------------------------------
## primer sequences
ITS1F <- "CTTGGTCATTTAGAGGAAGTAA"  
ITS2 <- "GCTGCGTTCTTCATCGATGC" 

## primer orientations
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)  # convert primer sequence in DNAString for Biostrings package handling
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # convert back to character vector
}
ITS1F.orients <- allOrients(ITS1F)
ITS2.orients <- allOrients(ITS2)

## primer removal 
cutadapt <- "/usr/bin/cutadapt" # path to cutadapt in system
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))

ITS1F.RC <- dada2:::rc(ITS1F)
ITS2.RC <- dada2:::rc(ITS2)

R1.flags <- paste("-g", ITS1F, "-a", ITS2.RC) 

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], # output files
                             fnFs[i])) # input files
} #for more information see arguments for cutadapt

# check the successful removal of primers (0 values should be returned on all counts)
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(ITS1F.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ForwardReads = sapply(ITS2.orients, primerHits, fn = fnFs.cut[[1]]))

## forward fastq filenames already processd by cutadapt
cutFs <- sort(list.files(path.cut, pattern = ".fastq.gz", full.names = TRUE))

## extract sample names for downstream use
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))

# Quality profiles ------------------------------------------------
## plot forward read quality scores
plotQualityProfile(cutFs[16]) 

# Filter and trim Sequences -----------------------------------------------
## create file path for filtered and trimmed reads to be written
filtFs <- file.path(path.cut, "filtered", basename(cutFs))

## Si5 pipeline maxEE = 2, minLen = 50, trimRight = 100 (50 selected as 100 is excessive with the quality of our sequences)
out_100 <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2, truncQ = 2, minLen = 50, trimRight = 100, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  #

# Learn error rates -------------------------------------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

# Infer ASVs --------------------------------------------------------------
## note: some files with no reads don't 'exist' at this stage so must not be called by filtFs or an error will be produced.
file.exists(filtFs)
filtFs.x <- filtFs[-c(559, 562,566, 568)]
file.exists(filtFs.x) ## because of this situation, we also need to drop these from sample names etc
sample.names.x <- sample.names[-c(559, 562,566, 568)]

## dereplicated forward reads
derepFs <- derepFastq(filtFs.x, verbose = TRUE)

## name the derepilicated objects by the sample names
names(derepFs) <- sample.names.x

## infer amplicon sequence variants using the DADA2 algorithm
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)

#start from here 2/12/21
## construct sequence table
seqtab.r1.100 <- makeSequenceTable(dadaFs) 

# Remove chimeras ---------------------------------------------------------
seqtab.nochim.r1.100 <- removeBimeraDenovo(seqtab.r1.100, method = "consensus", multithread = TRUE, verbose = TRUE)

# Pipeline retention -----------------------------------------------------
## inspect the length distribution of sequences
table(nchar(getSequences(seqtab.nochim.r1.100)))

## evaluate loss of sequences throughout pipeline stage
getN <- function(x) sum(getUniques(x))
track.r1.100 <- cbind(out_100[-c(559, 562,566, 568),], sapply(dadaFs, getN), rowSums(seqtab.nochim.r1.100))
colnames(track.r1.100) <- c("input", "filtered", "denoisedF", 
                        "nonchim")
rownames(track.r1.100) <- sample.names.x
track.r1.100

# Taxonomic Assignment ----------------------------------------------------
## path to reference database
unite.ref <- "~/DBs/sh_general_release_04.02.2020/sh_general_release_dynamic_04.02.2020.fasta"

## taxonomic assignment, with a minimum bootstrap confidence of 50 
taxa.r1.100 <- assignTaxonomy(seqtab.nochim.r1.100, unite.ref, multithread = TRUE, tryRC = TRUE, minBoot = 50, outputBootstraps = TRUE, verbose = TRUE)

## extract complete sequence list for posterity
sequence_saving_100 <- rownames(taxa.r1.100$tax)

## remove sequence names 
taxa.print.r1.100 <- taxa.r1.100$tax  
rownames(taxa.print.r1.100) <- NULL

# Plot taxonomic resolution -----------------------------------------------
taxa.print.r1.100 <- as.data.frame(taxa.print.r1.100)

## vector of assigned reads at each rank
reads.as.r1.100 <- as.vector(c(sum(!is.na(taxa.print.r1.100$Kingdom)), sum(!is.na(taxa.print.r1.100$Phylum)),sum(!is.na(taxa.print.r1.100$Class)),sum(!is.na(taxa.print.r1.100$Order)),sum(!is.na(taxa.print.r1.100$Family)),sum(!is.na(taxa.print.r1.100$Genus)),sum(!is.na(taxa.print.r1.100$Species))))

## vector of unassgined reads at each rank
reads.un.r1.100 <- as.vector(c(sum(is.na(taxa.print.r1.100$Kingdom)), sum(is.na(taxa.print.r1.100$Phylum)),sum(is.na(taxa.print.r1.100$Class)),sum(is.na(taxa.print.r1.100$Order)),sum(is.na(taxa.print.r1.100$Family)),sum(is.na(taxa.print.r1.100$Genus)),sum(is.na(taxa.print.r1.100$Species))))

## create plotting dataframe
rank <- as.data.frame(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
bind.as.r1.100 <- cbind(reads.as.r1.100, as.vector(rank))
bind.as.r1.100$type <- "Assigned"
colnames(bind.as.r1.100) <- c("value", "rank", "type")
bind.un.r1.100 <- cbind(reads.un.r1.100, as.vector(rank))
bind.un.r1.100$type <- "Unassigned"
colnames(bind.un.r1.100) <- c("value", "rank", "type")
df.pipe.r1.100 <- rbind(bind.as.r1.100, bind.un.r1.100)

## level fators for plotting
df.pipe.r1.100$rank <- factor(df.pipe.r1.100$rank, levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
df.pipe.r1.100$type <- factor(df.pipe.r1.100$type, levels = c("Unassigned", "Assigned"))

## create and save plot
p.pipe.assign.100 <- ggplot(df.pipe.r1.100, aes(x=rank, y=value, fill=type)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("gray80","firebrick1")) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("") +
  ylab("ASVs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

save_plot("l4_pipeline_assignment_pauvert_100.png", p.pipe.assign.100, base_width = 8, base_aspect_ratio = 1.3)



# Sample metadata ---------------------------------------------------------
## constructed metadata
sd.l4 <- as_tibble(substr(rownames(seqtab.nochim.r1.100),1,regexpr("_",rownames(seqtab.nochim.r1.100))-1))
colnames(sd.l4) <- "sample"
sd.l4$no_mba <- gsub("S","",sd.l4$sample)

## imported metadata
l4_meta <- readRDS("sample_ctd_meta.rds")
l4_meta$no_mba <- as.character(l4_meta$no_mba)

## join metadata objects
l4_meta_join <- left_join(sd.l4, l4_meta, by = "no_mba") %>% mutate(month = month(date)) %>% mutate(week = week(date))

## create a vector for sample names
s.vec <- as.vector(1:length(rownames(l4_meta_join)))  
s.nam <- cbind("sample_", s.vec)
s.nam <- as.data.frame(s.nam)
s.names <- paste0(s.nam$V1, s.nam$s.vec)
s.names <- as.data.frame(s.names)

## apply sample names to metadata 
row.names(l4_meta_join) <- s.names$s.names
l4_meta_join <- as.data.frame(l4_meta_join)

## apply sample names to sequence table
row.names(seqtab.nochim.r1.100) <- s.names$s.names

# ASV names ---------------------------------------------------------------
dim(taxa.print.r1.100)
dim(seqtab.nochim.r1.100)
a.vec <- as.vector(1:length(colnames(seqtab.nochim.r1.100)))  #number should reflect your total ASVs
a.nam <- cbind("asv_", a.vec)
a.nam <- as.data.frame(a.nam)
asv.names <- paste0(a.nam$V1, a.nam$a.vec)
asv.names <- as.data.frame(asv.names)

colnames(seqtab.nochim.r1.100) = asv.names$asv.names
taxa.l4.100 <- tax_table(taxa.print.r1.100)
rownames(taxa.l4.100) <- asv.names$asv.names
colnames(taxa.l4.100) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Phyloseq object ---------------------------------------------------------
l4.phy.100 <- phyloseq(tax_table(taxa.l4.100), otu_table(seqtab.nochim.r1.100, taxa_are_rows = FALSE), sample_data(l4_meta_join))

# Optimal taxonomy --------------------------------------------------------
# Restructuring tax_table to include 'best' classification
library(zoo)
library(tibble)
bc.t = t(as.data.frame(tax_table(l4.phy.100)))
bc.t[bc.t=="k__"] <- ""
bc.t[bc.t=="p__"] <- ""
bc.t[bc.t=="c__"] <- ""
bc.t[bc.t=="o__"] <- ""
bc.t[bc.t=="f__"] <- ""
bc.t[bc.t=="g__"] <- ""
bc.t[bc.t=="s__"] <- ""
bc.t[bc.t==""] <- NA
bc.fill <- na.locf(bc.t, na.rm = TRUE)
t.bc.fill <- as.data.frame(t(bc.fill))
head(t.bc.fill)
rnc.bc <- rownames_to_column(t.bc.fill, "ASV")
rnc.bc <- as.data.frame(rnc.bc)

## creates a column with the best classification and the ASV
rnc.bc$taxa_ASV <- paste(rnc.bc$Genus,rnc.bc$Species,rnc.bc$ASV)

## bind this column back onto the original tax_table 
safe.bc <- as.data.frame(tax_table(l4.phy.100))
safe.bc$taxa_ASV <- paste(rnc.bc$taxa_ASV)

# setup object as tax_table
bc.tax <- tax_table(safe.bc)
colnames(bc.tax) <- colnames(safe.bc)
rownames(bc.tax) <- rownames(safe.bc)

## update phyloseq object with new table
identical(bc.tax[1:length(colnames(seqtab.nochim.r1.100)),1:7], tax_table(l4.phy.100)) #should be true
tax_table(l4.phy.100) <- bc.tax


# Phyloseq processing ----------------------------------------------------
## remove taxa which are unassigned at the Phylum level
l4.phy.fun.100 <- subset_taxa(l4.phy.100, !is.na(Phylum))

# Read count histogram ----------------------------------------------------
hist.df.100 <- as.data.frame(rbind(cbind(rowSums(otu_table(l4.phy.100)), "total"),
                               cbind(rowSums(otu_table(l4.phy.fun.100)), "fungal")))
colnames(hist.df.100) <- c("reads", "type")
hist.df.100$reads <- as.numeric(as.character(hist.df.100$reads))

read.plot.100 <- ggplot(hist.df.100, aes(x=reads, fill = type)) + 
  geom_histogram(binwidth=5000) + theme_minimal() + xlab("Read Count") + ylab("Frequency") + scale_fill_manual(values = c("firebrick1","grey80")) +
  scale_y_continuous(expand = c(0,0)) 
  
save_plot("l4_100_pipeline_reads.png", read.plot.100, base_width = 8, base_aspect_ratio = 1.3)

## remove bad samples with <100 sequecnes
l4.phy.fun.100.GOOD <- prune_samples(sample_sums(l4.phy.fun.100) > 100, l4.phy.fun.100)

# ASV sequences -----------------------------------------------------------
length(sequence_saving_100)

seq_fungi_100 <- sequence_saving_100[as.numeric(gsub("asv_","",taxa_names(l4.phy.fun.100.GOOD)))]

## saving sequences of retained ASVs
write.fasta(as.list(seq_fungi_100), taxa_names(l4.phy.fun.100.GOOD), file = "l4_100_fungi_seqs.fasta")

## with 0.01% (rare) taxa removed
FSr  = transform_sample_counts(l4.phy.fun.100.GOOD, function(x) x / sum(x) )
FSfr = filter_taxa(FSr, function(x) sum(x) > .01, TRUE)
l4_clean <-FSfr

seq_fungi_100_0.01 <- sequence_saving_100[as.numeric(gsub("asv_","",taxa_names(l4_clean)))]
write.fasta(as.list(seq_fungi_100_0.01), taxa_names(l4_clean), file = "l4_100_fungi_seqs_001_removed.fasta")


write.csv(otu_table(l4_clean),"asv_abundance.csv")
