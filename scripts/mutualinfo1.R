## Calculation of mutual information from sampled data
## test 1

## libraries and colour-blind palette from http://www.sron.nl/~pault/
library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('png')
library('plot3D')
library('doParallel')
library('GA')
library('dplyr')

mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mycolours <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mycolours)
barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
dev.off()
mmtoin <- 0.0393701
plotsdir <- './'


## The script calculates the mutual information as in the genobayes.pdf notes, by calculating the frequencies for the combined variate and for the two variates separately in the data. Here we do this simply by counting. We choose a subset of genetic variations to check for calculation speed

## load all data
dpath = "./"
fs = dir(path = dpath,pattern = "dataset1")
d <- read.csv(paste0(dpath,fs[1]),sep=";")

## insomnia symptoms = col 6
## gene data = cols 8:102
ngenes <- 3
data <- d[,c(6,8:(7+ngenes))]
l <- length(data[1,])

## data frequencies for symptoms
fs <- tally(group_by_at(data,.vars=c(1)))
## data frequencies for genes
fg <- tally(group_by_at(data,.vars=c(2:l)))
## data frequencies for both
fx <- tally(group_by_at(data,.vars=c(1:l)))
