# load libraries
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
tbl <- read.table("./lysogenic_fna/final_results2",header=FALSE,row.names=1)
colnames(tbl) <- rownames(tbl)

# plot k-means using 2 clusters (from elbow method)
tbl2 <- as.matrix(tbl)
#p2 <- kmeans(tbl2, 2)
#fviz_cluster(p2, data = tbl2, geom="point")

# plot k-means using 4 clusters (from silhouette method)
p4 <- kmeans(tbl2, 4,nstart=20)
fviz_cluster(p4, data = tbl2, geom="point")

# plot k-means using 10 clusters (from gap statistic)
#p1 <- kmeans(tbl2, 10)
#fviz_cluster(p1, data = tbl2, geom="point")

# load in the pseudomonas OTU5 phylogeny 
treex2 <- read.tree(file="./collapsed_otu5.nwk")

# cleaning the cluster data
# make p1 into a data frame to see if we can map the values to a phylogeny
clusters <- data.frame(p4$cluster)
fragmented_cluster_info <- data.frame(p4$cluster)
clusters$genome_ID <- row.names(clusters)
#clusters[] <- lapply(clusters, gsub, pattern='late', replacement='')
clusters$genome_ID <- gsub("^.{0,4}", "", clusters$genome_ID)
clusters[] <- lapply(clusters, gsub, pattern='.fna', replacement='')
rownames(clusters) <- NULL
clusters <- clusters[,c(2,1)]

# Add a column to clusters: Phage #
clusters['phages'] <- NA
df <- clusters %>% group_by(genome_ID) %>% dplyr::mutate(phages = row_number())

# analysis of phage data
# find largest # of phages in column phages
max(df$phages) # 8 is the maximum number of phages
ggplot(df, aes(x=phages)) + geom_histogram() + ggtitle("Phage distribution") 

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

# cleaning the tree
# remove rows in treex2 that aren't in df1 since that will have all genomes we need
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

#write.csv(all_df_matrix, "./lysogenic_ffn/lysogenic_circos_table.csv")

p <- ggtree(tree, layout='circular', branch.length='none')
p
p3 <- gheatmap(p, all_df_matrix, colnames = FALSE, legend_title = "Cluster")
p3
pal <- c("#ada77c","#77AADD", "#117777", "#771122", "#771155")
p3 + scale_fill_manual(values = pal)
#p3 + scale_fill_viridis(discrete = TRUE) + geom_strip('p9.C4','p4.C6', fontsize = 7,barsize = 1, label = "Other\nPseudomonads", offset = 90, offset.text=5) +
  #geom_strip('p27.D6','p5.C1', fontsize = 7, barsize=1,label='OTU5', offset=90, offset.text=40) + theme(legend.position = "left") + guides(fill=guide_legend(title= "Phage Genetic Cluster"))


# Analyzing each cluster separately
# Preparing the data for analysis
# subset each cluster to its own dataframe
df4 <- subset(fragmented_cluster_info, p4.cluster==4)
df3 <- subset(fragmented_cluster_info, p4.cluster==3)
df2 <- subset(fragmented_cluster_info, p4.cluster==2)
df1 <- subset(fragmented_cluster_info, p4.cluster==1)

#write.csv(df4, "./lysogenic_fna//cluster4_list.csv")
#write.csv(df3, "./lysogenic_fna/cluster3_list.csv")
#write.csv(df2, "./lysogenic_fna/cluster2_list.csv")
#write.csv(df1, "./lysogenic_fna/cluster1_list.csv")


# subset each mash distance matrix to its own cluster
df1.1 <- tbl[rownames(tbl) %in% rownames(df1), rownames(tbl) %in% rownames(df1)] 
df2.1 <- tbl[rownames(tbl) %in% rownames(df2), rownames(tbl) %in% rownames(df2)] 
df3.1 <- tbl[rownames(tbl) %in% rownames(df3), rownames(tbl) %in% rownames(df3)] 
df4.1 <- tbl[rownames(tbl) %in% rownames(df4), rownames(tbl) %in% rownames(df4)] 
# convert to data matrix 
df1.1 <- as.matrix(df1.1)
df2.1 <- as.matrix(df2.1)
df3.1 <- as.matrix(df3.1)
df4.1 <- as.matrix(df4.1)

# heatmaps
# visualize the mash distance within a cluster
image(df1.1)
image(df2.1)
image(df3.1)
image(df4.1)

