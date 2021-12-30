library(ggplot2)
dat <- read.csv("./mitomycinC_trials2.csv", header = T, sep = ",")

ggplot(dat, aes(x=Time..s., y=OD600,color=Sample)) + geom_point() + geom_smooth(se=F)

ggplot(dat, aes(x=Time..s., y=OD600)) + geom_point() + geom_smooth(se=F) + facet_wrap(~Sample)


ggplot(dat, aes(x=Time..s., y=OD600, color=Sample)) + geom_point() + geom_smooth() + facet_wrap(~Sample)
