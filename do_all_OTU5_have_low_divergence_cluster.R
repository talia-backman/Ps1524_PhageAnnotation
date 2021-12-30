library(cluster)
library(NbClust)
library(clustertend)
library(factoextra)
library(stats)
library(ape)
library(phytools)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(dplyr)
library(BBmisc)


# set the seed
set.seed(123)

# load and tidy phage data from mash 
tbl <- read.table("./lysogenic_fna/lysogenic_fna_final_results",header=FALSE,row.names=1)
colnames(tbl) <- rownames(tbl)

# plot k-means using 2 clusters (from elbow method)
tbl2 <- as.matrix(tbl)
#p2 <- kmeans(tbl2, 2)
#fviz_cluster(p2, data = tbl2, geom="point")

# plot k-means using 4 clusters (from silhouette method)
p4 <- kmeans(tbl2, 4,nstart=20)
fviz_cluster(p4, data = tbl2, geom="point")

# load in the pseudomonas OTU5 phylogeny 
treex2 <- read.tree(file="./ps_1524_uncollapsed_5_2018.nwk")

# cleaning the cluster data
# make p1 into a data frame to see if we can map the values to a phylogeny
clusters <- data.frame(p4$cluster)
fragmented_cluster_info <- data.frame(p4$cluster)
clusters$genome_ID <- row.names(clusters)
clusters[] <- lapply(clusters, gsub, pattern='late', replacement='')
clusters$genome_ID <- gsub("^.{0,4}", "", clusters$genome_ID)
clusters[] <- lapply(clusters, gsub, pattern='.pilon.contigs_renamed.phages_lysogenic.fna', replacement='')
rownames(clusters) <- NULL
clusters <- clusters[,c(2,1)]

# Add a column to clusters: Phage #
clusters['phages'] <- NA
df <- clusters %>% group_by(genome_ID) %>% dplyr::mutate(phages = row_number())

# cleaning the cluster data
# subset df based off phages column, rename their cluster column
df1 <- subset(df, phages == 1) #1582
names(df1)[2] <- "phage1"
df2 <- subset(df, phages == 2) # 884
names(df2)[2] <- "phage2"
df3 <- subset(df, phages == 3) # 443
names(df3)[2] <- "phage3"
df4 <- subset(df, phages == 4) # 149
names(df4)[2] <- "phage4"
df5 <- subset(df, phages == 5) # 40
names(df5)[2] <- "phage5"
df6 <- subset(df, phages == 6) # 4
names(df6)[2] <- "phage6"
df7 <- subset(df, phages == 7) # 1
names(df7)[2] <- "phage7"
df8 <- subset(df, phages == 8) # 1
names(df8)[2] <- "phage8"

# see how many times we see each genome_ID
#table(df$genome_ID)

uneeded <- which(!(treex2$tip.label) %in% df1$genome_ID)
length(uneeded)
tree <- drop.tip(treex2, uneeded)

# make comprehensive phage clustering dataframe with the corresponding genome as rownames
# make all 15 dataframes 1510 observations long but with NA in places that don't have values
temp <- merge(df1, df2, by="genome_ID", all=TRUE)
temp <- merge(temp,df3, by="genome_ID", all=TRUE)
temp <- merge(temp,df4, by="genome_ID", all=TRUE)
temp <- merge(temp,df5, by="genome_ID", all=TRUE)
temp <- merge(temp,df6, by="genome_ID", all=TRUE)
temp <- merge(temp,df7, by="genome_ID", all=TRUE)
all_df <- merge(temp,df8, by="genome_ID", all=TRUE)

all_df <- all_df[ ,c(1,2,4,6,8,10,12,14,16)]

# convert to matrix 
all_df_matrix <- as.matrix(all_df)
rownames(all_df_matrix) <- all_df_matrix[,1]
all_df_matrix <- all_df_matrix[ ,c(2:9)]





############# IS THE LOWEST DIVERGENT CLUSTER PHAGE PRESENT IN ALL OTU5?
all_df_matrix[all_df_matrix=="2"] <- "YES"
all_df_matrix[all_df_matrix=="3"] <- "NO"
all_df_matrix[all_df_matrix=="4"] <- "NO"
all_df_matrix[all_df_matrix=="1"] <- "NO"

all_df_matrix <- as.data.frame(all_df_matrix)
has_lowest_divergent_cluster_phage1 <- (which(all_df_matrix$phage1=="YES"))
has_lowest_divergent_cluster_phage2 <- (which(all_df_matrix$phage2=="YES"))
has_lowest_divergent_cluster_phage3 <- (which(all_df_matrix$phage3=="YES"))
has_lowest_divergent_cluster_phage4 <- (which(all_df_matrix$phage4=="YES"))
has_lowest_divergent_cluster_phage5 <- (which(all_df_matrix$phage5=="YES"))



has_lowest_divergent_cluster <- (which(all_df_matrix=="YES"))

which(all_df_matrix[]=="YES")
has_lowest_divergent_cluster
OTU5 <- 106:1513
#OTU5 <- tree$tip.label[106:1513]

hasit <- which(OTU5 %in% has_lowest_divergent_cluster)
doesnthaveit <- which(!OTU5 %in% has_lowest_divergent_cluster)





has_lowest_divergent_cluster <- apply(all_df_matrix, 1, function(r) any(r %in% "YES"))
test <- as.data.frame(has_lowest_divergent_cluster)
test <- tibble::rownames_to_column(test)

OTU5_tiplabels <- tree$tip.label
OTU5_tiplabels <- as.data.frame(OTU5_tiplabels)
OTU5_tiplabels$number <- 1:1513
OTU5_tiplabels <- tail(OTU5_tiplabels, 1407)

colnames(OTU5_tiplabels) <- c("genome","number")
colnames(test) <- c("genome","has_lowest_divergent_cluster")
  
dat2 <- merge(OTU5_tiplabels, test)

p23_d7 <- subset(dat2, genome=="p23.D7")









################################################################
############### plot clusters to phylogeny #####################
################################################################
p <- ggtree(tree, layout = "circular",branch.length='none')
p
p3 <- gheatmap(p, all_df_matrix, colnames = FALSE, legend_title = "Cluster")
p3
p3 + geom_strip('p23.D7','p1.H11', barsize=2, color='purple',label='OTU5', offset=80, offset.text=1)
