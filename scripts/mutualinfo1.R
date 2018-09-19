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
genes <- 1:15
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
cgminus <- cg-dim(fg)[1]


minfo <- log(ncnn) + sum(sapply(1:8,function(sig){
    apply(fg,1,function(gam){
        ## marginal freqs
        ns <- fs[sig,2]
        ng <- gam[ngenes+1]
        ## frequency of pair
        nx <- sum(apply(fx,1,function(i){all(i[c(1,genes+1)]==c(sig,gam[genes]))*i[ngenes+2]}))
        numer <- nnc*nx+dn

        numer/ncnn*(log(numer)-log(nnc*ns+dn*cg)-log(nnc*ng+dn*cs))})})) -
    cgminus*dn/ncnn*(log(cs) + sum(log(nnc*fs[,2]+dn*cg)))



sentropy <- -sum(sapply(1:8,function(sig){
    prob <- (nnc*fs[sig,2]+dn*cg)/ncnn
    prob*log(prob)}))

c(minfo,sentropy,minfo/sentropy)
