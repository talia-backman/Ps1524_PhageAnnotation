library(ggplot2)
library(tidyr)

dat <- read.csv("./12_01_2021.csv",header=T,sep="\t")

# transform dataframe from wide to long format
df.long <- pivot_longer(dat, cols=2:97, names_to = "Sample", values_to = "OD600")

# subset df to only the samples/wells we used
df.long <- subset(df.long, Sample == "C11" | Sample == "F9" | Sample == "B9" | Sample == "H5" | Sample == "E2" | Sample == "B10" | Sample == "A6" | Sample == "F6" | Sample == "D6" | 
                    Sample == "E12" | Sample == "C9" | Sample == "G2" | Sample == "E1" | Sample == "G5" | Sample == "A5" | Sample == "C1" | Sample == "F2" | Sample == "H6" | 
                    Sample == "E6" | Sample == "C10" | Sample == "D5" | Sample == "B6" | Sample == "C6" | Sample == "C4")

df.long[df.long == "C11"] <- "p13.D10"
df.long[df.long == "F9"] <- "p13.D10"
df.long[df.long == "B9"] <- "p13.D10"
df.long[df.long == "H5"] <- "p20.H5"
df.long[df.long == "E2"] <- "p20.H5"
df.long[df.long == "B10"] <- "p20.H5"
df.long[df.long == "A6"] <- "p21.F9"
df.long[df.long == "F6"] <- "p21.F9"
df.long[df.long == "D6"] <- "p21.F9"
df.long[df.long == "E12"] <- "p23.D8"
df.long[df.long == "C9"] <- "p23.D8"
df.long[df.long == "G2"] <- "p23.D8"
df.long[df.long == "E1"] <- "p22.B2"
df.long[df.long == "G5"] <- "p22.B2"
df.long[df.long == "A5"] <- "p22.B2"
df.long[df.long == "C1"] <- "p23.D7"
df.long[df.long == "F2"] <- "p23.D7"
df.long[df.long == "H6"] <- "p23.D7"
df.long[df.long == "E6"] <- "p25.A12"
df.long[df.long == "C10"] <- "p25.A12"
df.long[df.long == "D5"] <- "p25.A12"
df.long[df.long == "B6"] <- "LB_ctrl"
df.long[df.long == "C6"] <- "LB_ctrl"
df.long[df.long == "C4"] <- "LB_ctrl"

# plot
ggplot(df.long, aes(x=Time..s., y=OD600,color=Sample)) + geom_point() + geom_smooth(se=F)

png("./figures/12_01_2021.png")
ggplot(df.long, aes(x=Time..s., y=OD600)) + geom_point() + geom_smooth(se=F) + facet_wrap(~Sample)
dev.off()

ggplot(df.long, aes(x=Time..s., y=OD600, color=Sample)) + geom_point() + geom_smooth() + facet_wrap(~Sample)
