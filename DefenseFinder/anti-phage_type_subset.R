library(ggplot2)
library(ggpubr)
library(ape)
library(ggtree)
library(dplyr)
library(pals)

# read in data
df <- read.csv("./all_defense_finder_systems.csv")

# look at phylogenetic relationships of anti-phage systems
tree <- read.tree("./ps_1524_uncollapsed_5_2018.nwk")

tree_df <- df[ ,c(3,10)]
tree_df['system_type_count'] <- NA
tree_df2 <- tree_df %>% group_by(bacterial_strain) %>% dplyr::mutate(system_type_count = row_number())

# subset to only most common anti-phage systems (type still)
tree_df2_common <- subset(tree_df2, type == "RM" | type == "Abi2")
tree_df2_common['system_type_count'] <- NA
tree_df2 <- tree_df2_common %>% group_by(bacterial_strain) %>% dplyr::mutate(system_type_count = row_number())
max(tree_df2$system_type_count) #10

# subset df based off system_type_count column
df1 <- subset(tree_df2, system_type_count == 1) 
names(df1)[1] <- "system1"
df2 <- subset(tree_df2, system_type_count == 2) 
names(df2)[1] <- "system2"
df3 <- subset(tree_df2, system_type_count == 3) 
names(df3)[1] <- "system3"
df4 <- subset(tree_df2, system_type_count == 4) 
names(df4)[1] <- "system4"
df5 <- subset(tree_df2, system_type_count == 5) 
names(df5)[1] <- "system5"
df6 <- subset(tree_df2, system_type_count == 6) 
names(df6)[1] <- "system6"
df7 <- subset(tree_df2, system_type_count == 7) 
names(df7)[1] <- "system7"
df8 <- subset(tree_df2, system_type_count == 8) 
names(df8)[1] <- "system8"
df9 <- subset(tree_df2, system_type_count == 9) 
names(df9)[1] <- "system9"
df10 <- subset(tree_df2, system_type_count == 10) 
names(df10)[1] <- "system10"

# make comprehensive system type dataframe with the corresponding genome as rownames
# make all 15 dataframes 1510 observations long but with NA in places that don't have values
temp <- merge(df1, df2, by="bacterial_strain", all=TRUE)
temp <- merge(temp,df3, by="bacterial_strain", all=TRUE)
temp <- merge(temp,df4, by="bacterial_strain", all=TRUE)
temp <- merge(temp,df5, by="bacterial_strain", all=TRUE)
temp <- merge(temp,df6, by="bacterial_strain", all=TRUE)
temp <- merge(temp,df7, by="bacterial_strain", all=TRUE)
temp <- merge(temp,df8, by="bacterial_strain", all=TRUE)
temp <- merge(temp,df9, by="bacterial_strain", all=TRUE)
all_df <- merge(temp,df10, by="bacterial_strain", all=TRUE)

all_df <- all_df[ ,c(1,2,4,6,8,10,12,14,16,18,20)]

# convert to matrix 
all_df_matrix <- as.matrix(all_df)
rownames(all_df_matrix) <- all_df_matrix[,1]
all_df_matrix <- all_df_matrix[ ,c(2:11)]


# plot DefenseFinder data
p <- ggtree(tree, layout='circular', branch.length='none')
p1 <- gheatmap(p, all_df_matrix, colnames = FALSE, legend_title = "Anti-phage System")
p1 + geom_strip('p9.C4','p4.C6', fontsize = 7,barsize = 1, label = "Other\nPseudomonads", offset = 90, offset.text=5) +
  geom_strip('p27.D6','p5.C1', fontsize = 7, barsize=1,label='OTU5', offset=90, offset.text=40) + 
  theme(legend.position = "left") + guides(fill=guide_legend(title= "Anti-phage System"))