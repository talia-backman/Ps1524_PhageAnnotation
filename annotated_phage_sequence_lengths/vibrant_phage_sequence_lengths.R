library(tidyverse)

# looking at sequence lengths of the phages
my_data <- read_csv("./combined/count_length.txt",col_names = FALSE)

# make a new df split by the space between each number in column 2
dat <- read_csv("./combined/count_length.txt",col_names = FALSE)
s <- strsplit(dat$X2, split = " ")
dat2 <- data.frame(X1 = rep(dat$X1, sapply(s, length)), X2 = unlist(s))

# subset to only genomes with 15,000-16,000 bp
phages_15kb <- subset(dat2, X2 >= 15000)
phages_15kb <- subset(phages_15kb, X2 <= 16000)
write_csv(phages_15kb, "./combined/15kb_phages_combined.csv")

# separate X2 into many columns
my_data$phage1 <- sapply(strsplit(as.character(my_data$X2),' '), "[",1)
my_data$phage2 <- sapply(strsplit(as.character(my_data$X2),' '), "[",2)
my_data$phage3 <- sapply(strsplit(as.character(my_data$X2),' '), "[",3)
my_data$phage4 <- sapply(strsplit(as.character(my_data$X2),' '), "[",4)
my_data$phage5 <- sapply(strsplit(as.character(my_data$X2),' '), "[",5)
my_data$phage6 <- sapply(strsplit(as.character(my_data$X2),' '), "[",6)
my_data$phage7 <- sapply(strsplit(as.character(my_data$X2),' '), "[",7)
my_data$phage8 <- sapply(strsplit(as.character(my_data$X2),' '), "[",8)
my_data$phage9 <- sapply(strsplit(as.character(my_data$X2),' '), "[",9)
my_data$phage10 <- sapply(strsplit(as.character(my_data$X2),' '), "[",10)
my_data$phage11 <- sapply(strsplit(as.character(my_data$X2),' '), "[",11)
my_data$phage12 <- sapply(strsplit(as.character(my_data$X2),' '), "[",12)
my_data$phage13 <- sapply(strsplit(as.character(my_data$X2),' '), "[",13)
my_data$phage14 <- sapply(strsplit(as.character(my_data$X2),' '), "[",14)
my_data$phage15 <- sapply(strsplit(as.character(my_data$X2),' '), "[",15)
my_data$phage16 <- sapply(strsplit(as.character(my_data$X2),' '), "[",16)
my_data$phage17 <- sapply(strsplit(as.character(my_data$X2),' '), "[",17)
my_data$phage18 <- sapply(strsplit(as.character(my_data$X2),' '), "[",18)

my_data <- my_data[c(1,3:20)]
names(my_data)[names(my_data) == "X1"] <- "genome"


max(my_data) # max is 431326
min(my_data$X1) # min is 1066
mean(my_data$X1) # mean is 27214.85
median(my_data$X1) # median is 18109

small_df <- subset(my_data, X1<34000)
big_df <- subset(my_data,X1>=34000)

hist(small_df)
hist(big_df)


hist(my_data$X1)



# subset data to only the 15,000 bp phage
df <- subset(my_data, phage1 >= 15000)
df <- subset(df, phage1 <= 16000)


test <- subset(my_data, phage2 >=15000)
test <- subset(test, phage2 <= 16000)



