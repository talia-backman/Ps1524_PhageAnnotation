library(ggplot2)
library(ggpubr)
library(ape)
library(ggtree)

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



p <- ggtree(tree, layout='circular', branch.length='none')
p
gheatmap(p, tree_df)
