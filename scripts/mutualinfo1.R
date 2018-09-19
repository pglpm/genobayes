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
genes <- 1:10
ngenes <- length(genes)
data <- d[,c(6,7+genes)]
l <- length(data[1,])
n <- length(data[,1])
nn <- 5.3e6L
cs <- 8
cg <- 2^ngenes
cx <- cs*cg

## recurring quantities
nnc <- (nn+cx)
ncnn <- (n+cx)*nn
dn <- nn-n


## data frequencies for symptoms
fs <- data.matrix(tally(group_by_at(data,.vars=c(1))))
## data frequencies for genes
fg <- data.matrix(tally(group_by_at(data,.vars=c(genes+1))))
## data frequencies for both
fx <- data.matrix(tally(group_by_at(data,.vars=c(1,genes+1))))

notc <- cx-dim(fx)[1]


a <- log(ncnn) + sum(apply(fx,1,function(i){
    ## frequency of pair
    nx <- i[ngenes+2]
    ## marginal freqs
    ns <- fs[i[1],2]
    ng <- fg[which(apply(fg[,genes],1,function(r){all(r == i[genes+1])})),ngenes+1]

    numer <- nnc*nx+dn

    numer/ncnn*(log(numer)-log(nnc*ns+dn*cg)-log(nnc*ng+dn*cs))
})) +
    notc*dn/ncnn*log(dn)
