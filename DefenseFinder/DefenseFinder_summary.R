library(ggplot2)
library(ggpubr)
library(ape)
library(ggtree)
library(dplyr)
library(pals)

# read in data
df <- read.csv("./all_defense_finder_systems.csv")

png("./DefenseFinder_figures/types_barplot.png", width = 1300, height = 1300)
p1 <- ggplot(df, aes(x=reorder(type,type,function(x)-length(x)))) +
  geom_bar(width =.5, fill = "#4477AA") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Anti-phage system: Type") + ylab("Total genomes encoding systems") #+
  theme(text = element_text(size = 30))
dev.off()

subset(df, type == "Cas") # strains with CRISPR-Cas are p1.E6, p13.F5, p21.B10
# all are CAS_Class1-Type-I

png("./DefenseFinder_figures/sub-types_barplot.png", width = 1300, height = 1300)
p2 <- ggplot(df, aes(x=reorder(subtype,subtype,function(x)-length(x)))) +
  geom_bar(width =.5, fill = "#4477AA") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Anti-phage system: Sub-type") + ylab("Total genomes encoding systems") #+
  theme(text = element_text(size = 20))
dev.off()

png("./DefenseFinder_figures/type_and_subtype.png", width = 1000, height = 800)
ggarrange(p1 + font("x.text", size = 13), p2 + font("x.text", size = 13), labels=c("A","B"),nrow = 2, ncol = 1)
dev.off()


# look at phylogenetic relationships of anti-phage systems
tree <- read.tree("./ps_1524_uncollapsed_5_2018.nwk")

tree_df <- df[ ,c(3,10)]
tree_df['system_type_count'] <- NA
tree_df2 <- tree_df %>% group_by(bacterial_strain) %>% dplyr::mutate(system_type_count = row_number())

# find largest # of defense systems in one genome
max(tree_df2$system_type_count) #12

# histogram of # defense systems in one genome
png("./DefenseFinder_figures/distribution_histogram.png", width = 500, height = 500)
ggplot(tree_df2, aes(x=system_type_count)) + geom_bar(width=.5,fill = "#4477AA") + xlab("Defense systems in one genome") +
  ylab("Count") + scale_x_continuous(breaks=c(1,3,5,7,9,11,13))
dev.off()

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
df11 <- subset(tree_df2, system_type_count == 11) 
names(df11)[1] <- "system11"
df12 <- subset(tree_df2, system_type_count == 12) 
names(df12)[1] <- "system12"

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
temp <- merge(temp,df10, by="bacterial_strain", all=TRUE)
temp <- merge(temp,df11, by="bacterial_strain", all=TRUE)
all_df <- merge(temp,df12, by="bacterial_strain", all=TRUE)

all_df <- all_df[ ,c(1,2,4,6,8,10,12,14,16,18,20,22,24)]

# convert to matrix 
all_df_matrix <- as.matrix(all_df)
rownames(all_df_matrix) <- all_df_matrix[,1]
all_df_matrix <- all_df_matrix[ ,c(2:13)]


# plot DefenseFinder data
#pal <- c("#ada77c", "#7c82ad", "#4477AA","#AA7744","#77AADD","#DDAA77","#117777","#771111","#44AAAA","#AA4444","#77CCCC","#CC7777","#117744","#771144","#44AA77","#AA4477","#114477","#774411","#88CCAA","#CC88AA","#777711","#111177","#AAAA44","#4444AA","#DDDD77","#7777DD","#774411","#114477","#AA7744","#4477AA","#DDAA77","#77AADD","#771122","#117766","#DD7788","darkgoldenrod1","#771155","#117733","#AA4488","#44AA66","#CC99BB","#99CCAA","#AA7744","#4477AA","#AA4455","#846897","#7B9768","#4A87B5","#B5784A")
p <- ggtree(tree, layout='circular', branch.length='none')
p1 <- gheatmap(p, all_df_matrix, colnames = FALSE, legend_title = "Anti-phage System")
p1 + geom_strip('p9.C4','p4.C6', fontsize = 7,barsize = 1, label = "Other\nPseudomonads", offset = 90, offset.text=5) +
  geom_strip('p27.D6','p21.B10', fontsize = 7, barsize=1,label='OTU5', offset=90, offset.text=40) + 
  theme(legend.position = "left") + guides(fill=guide_legend(title= "Anti-phage System"))
