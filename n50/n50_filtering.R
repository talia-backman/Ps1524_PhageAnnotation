dat <- read.csv("./n50/n50.csv", sep = " ")
dat2 <- subset(dat, N50_bp >= 150000)
dat3 <- subset(dat, N50_bp < 150000)

write.csv(dat2, "./n50/genomes_to_keep.csv")
write.csv(dat3, "./n50/genomes_to_remove.csv")
