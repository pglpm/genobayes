## Calculation of mutual information from sampled data
## for 1 gene, 2 genes, etc, keeping at each step those with highest mutual info
## parameter A = 2

## libraries and colour-blind palette from http://www.sron.nl/~pault/
##runfunction <- function(aa=1e3L){
library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('png')
library('plot3D')
library('doParallel')
#library('GA')
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
savedir <- './'

xlogy <- function(x,y){if(x==0){0}else{x*log(y)}}
entropyf <- function(x){-sum(sapply(x,function(y){xlogy(y,y)}))}

dpath  <-  "./"
datafile <- 'dataset1_simple.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
d <- read.csv(paste0(dpath,nfile[1]))
n <- length(d[,1])

filename <- 'condprobs_1gene'
aa <- 0

n2 <- n+aa


allgenes <- 1:94
allsymptoms <- 1:3
cores <- 1

priorj <- matrix(1,2,2)/4
priorj1 <- priorj[2,]
priors <- matrix(1,2,3)/2

sfreq <- foreach(symp=allsymptoms, .combine=cbind) %do% {table(d[,symp])}
sprob <- (sfreq + aa*priors)/n2
sprob1 <- sprob[2,]
##sentropy <- apply(sprob,2,entropyf)

cl <- makeCluster(cores)
registerDoParallel(cl)

result <- foreach(gene=allgenes, .combine=cbind,.export=c('xlogy','n','d','priorj','allgenes','aa'), .packages=c('dplyr')) %do% {
    conds <- sapply(1:3,function(symp){
        jfreq <- table(d[,c(symp,3+gene)])
        jprob <- (jfreq + aa*priorj)/n2
        gprob <- apply(jprob,2,sum)
        jprob[2,]/gprob
    })
    -(sprob1-t(conds))/sprob1}
dim(result) <- c(3,2,length(allgenes))


inds <- sapply(allsymptoms,function(symp){arrayInd(which.max(abs(result[symp,,])), dim(result[symp,,]))})

indspos <- sapply(allsymptoms,function(symp){arrayInd(which.max((result[symp,,])), dim(result[symp,,]))})

stopCluster(cl)
