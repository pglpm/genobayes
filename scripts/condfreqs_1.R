## Calculation of expected value, std, and quantiles of conditional frequencies
## of symptoms given gene+allele
## Based on mackayetal1995. Uses a Dirichlet distribution for the conditional frequencies

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
#
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

## xlogy <- function(x,y){if(x==0){0}else{x*log(y)}}
## entropyf <- function(x){-sum(sapply(x,function(y){xlogy(y,y)}))}

dpath  <-  "./"
datafile <- 'dataset1_binarized.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
data <- read.csv(paste0(dpath,nfile[1]))
n <- length(data[,1])

filename <- 'condfreqs_1gene' # where to save the results
aa <- 0 # parameter A for Dirichlet

n2 <- n+aa

allgenes <- 1:94
allsymptoms <- 1:3
quantiles <- c(0.05,0.95)
cores <- 25

## row an column headers for results
gnames <- colnames(data)[allgenes+3]
rown <- c(rbind(paste0(gnames,'-AA'),paste0(gnames,'-Bx')))
coln <-  c('EV','STD','Q.05','Q.95')

## prior expected frequencies for each symptom (parameter alpha)
allsprior <- matrix(1,nrow=2,ncol=3)/2

cl <- makeCluster(cores)
registerDoParallel(cl)

for(i in allsymptoms){
    sdata <- data[,c(i,3+allgenes)]
    sprior <- allsprior[,i]
    res <- t(foreach(g=allgenes,
                     .combine=cbind, .export=c('sdata','sprior','aa')) %do% {
                       sapply(0:1,function(x){
                           f <- table(sdata[sdata[[1+g]]==x,c(1,1+g)])
                           a2 <- sum(f)+aa
                           f2 <- (f+aa*sprior)/a2
                           sstd <- sqrt(prod(f2)/(1+a2))
                           squant <- qbeta(quantiles,a2*f2[2],a2*f2[1])
                           c(f2[2],sstd,squant)
                       })
                     })
    rownames(res) <- rown
    colnames(res) <- coln

    write.csv(res,paste0(dpath,filename,'_s',c('A','B','C')[i],'.csv'))
    ##,sep=',',row.names=rown,col.names=coln,na='NA')
}

stopCluster(cl)
