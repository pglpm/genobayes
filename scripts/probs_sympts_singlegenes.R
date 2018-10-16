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

dpath  <-  "./"
datafile <- 'dataset1_simple.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
d <- read.csv(paste0(dpath,nfile[1]))
n <- length(d[,1])

filename <- 'condprobs_1gene'
aa <- 0

n2 <- n+aa

priorjfreq <- matrix(1/4,2,2)

#fs <- foreach(symp=1:3, .combine=rbind) %do% {table[d[,symp]]/n}

allgenes <- 1:94
cores <- 30

cl <- makeCluster(cores)
registerDoParallel(cl)

result <- foreach(gene=allgenes, .combine=cbind,.export=c('xlogy','n','d','priorjfreq','allgenes','aa'), .packages=c('dplyr','doParallel')) %do% {

    parres <- foreach(symp=1:3, .combine=rbind) %do% {
        ## rows: symptom 0/1; columns: allele 0/1
        jfreq <- table(d[,c(symp,3+gene)])

        jprob <- (jfreq + aa*priorjfreq)/n2
        sprob <- apply(jprob,1,sum)
        gprob <- apply(jprob,2,sum)
        cprob <- t(t(jprob)/gprob)
        minfo <- sum(apply(cbind(c(jprob),c(jprob)/c(sprob %*% t(gprob))),1,function(x){xlogy(x[1],x[2])}))
        cbind(c(gene,minfo),sprob,cprob)
    }
    parres}

stopCluster(cl)

write.table(result,paste0(savedir,filename,'.csv'),sep=',',row.names=F,col.names=F)
