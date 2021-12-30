library(ggplot2)
# n50 for 1,534 Germany Pseudomonas genomes
n50 <- read.table("./pseudomonas_n50/n50.txt",sep = " ")

# subset to only n50 values
n50_df <- n50[c(1,19)]

# rename columns
names(n50_df)[names(n50_df) == "V19"] <- "n50_score"
names(n50_df)[names(n50_df) == "V1"] <- "genome"

# n50 is the sequence length of the shortest contig at 50% of the total genome length

mean(n50_df$n50_score) # 416652.5
min(n50_df$n50_score) # 977
max(n50_df$n50_score) # 1451350

# look at distribution of all n50 scores
hist(n50_df$n50_score)

# look at n50 of all "genomes with no phages" and "OTU5 without a cluster2 phage"
nophages <- subset(n50_df, genome == c("plate22.D2.pilon.contigs_renamed.fasta","plate22.D5.pilon.contigs_renamed.fasta",
                                     "plate22.D7.pilon.contigs_renamed.fasta","plate23.D8.pilon.contigs_renamed.fasta",
                                     "plate26.F6.pilon.contigs_renamed.fasta","plate26.F7.pilon.contigs_renamed.fasta",
                                     "plate2.H2.pilon.contigs_renamed.fasta","plate6.H5.pilon.contigs_renamed.fasta",
                                     "plate7.G9.pilon.contigs_renamed.fasta","plate9.H9.pilon.contigs_renamed.fasta",
                                     "plate22.D8.pilon.contigs_renamed.fasta"))

nophages <- subset(n50_df, genome == "plate22.D2.pilon.contigs_renamed.fasta" | genome == "plate22.D5.pilon.contigs_renamed.fasta" | 
                     genome == "plate22.D7.pilon.contigs_renamed.fasta" | genome == "plate23.D8.pilon.contigs_renamed.fasta" | genome ==
                   "plate26.F6.pilon.contigs_renamed.fasta" | genome == "plate26.F7.pilon.contigs_renamed.fasta" | genome ==
                   "plate2.H2.pilon.contigs_renamed.fasta" | genome == "plate6.H5.pilon.contigs_renamed.fasta" | genome ==
                   "plate7.G9.pilon.contigs_renamed.fasta" | genome == "plate9.H9.pilon.contigs_renamed.fasta" | genome ==
                   "plate22.D8.pilon.contigs_renamed.fasta")

hist(nophages$n50_score)


nocluster2phage <- subset(n50_df, genome == "plate23.D7.pilon.contigs_renamed.fasta")
hist(nocluster2phage$n50_score)

ggplot(nocluster2phage, aes(x=genome, y=n50_score)) + geom_point()
