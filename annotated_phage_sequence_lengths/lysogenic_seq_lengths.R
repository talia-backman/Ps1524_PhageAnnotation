library(tidyverse)
library(tidyr)

# looking at sequence lengths of the phages
my_data <- read_csv("./lysogenic/z_lysogenic_count_length.txt",col_names = FALSE)
head(my_data)
my_data <- my_data %>% column_to_rownames(var="X1")



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

my_data <- my_data[c(2:13)]
#names(my_data)[names(my_data) == "X1"] <- "genome"

df1 <- lapply(my_data,as.numeric)

# max value
colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMax(df1) # largest 'lysogenic' phage is 431,272

# min value
colMin <- function(data) sapply(data, min, na.rm = TRUE)
colMin(df1) # smallest 'lysogenic' phage is 3,933


# plot histogram of reads we see
dat <- read_csv("./lysogenic/z_lysogenic_count_length.txt",col_names = FALSE)
dat2 <- data.frame(X2 = unlist(strsplit(as.character(dat$X2), " ")))
dat3 <- lapply(dat2,as.numeric)
dat3 <- as.data.frame(dat3)

hist(dat3$X2)

ggplot(dat3, aes(x=X2)) + geom_histogram() + xlim(3900,100000)

ggplot(dat3, aes(x=X2)) + geom_histogram() + xlim(100000,432000)
