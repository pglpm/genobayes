## Calculation of expected value, std, and quantiles of conditional frequencies
## of symptoms given gene+allele pairs
## Based on mackayetal1995. Uses a Dirichlet distribution for the conditional frequencies
## Calculates parameter Theta from the data, as in mackayetal1995

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
#library('gridExtra')
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
tobinary  <- function(number, noBits) {
    as.numeric(intToBits(number))[1:noBits]
#    binary_vector[-(1:(length(binary_vector) - noBits))]
}
frombinary  <- function(digits) {sum(2^(0:(length(digits)-1))*digits)}

dpath  <-  "./"
datafile <- 'dataset1_binarized.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
data <- read.csv(paste0(dpath,nfile[1]))[,-1]
n <- length(data[,1])

savedir <- './2gene-results/'
filename <- 'genes' # where to save the results
allgenes <- 1:94
allsymptoms <- 1:3
quantiles <- c(0.05,0.95)
cores <- 1

## row an column headers for results
alnames <- c('A','B')
gnames <- colnames(data)[allgenes+3]
ganames <- c(rbind(paste0(gnames,'-AA'),paste0(gnames,'-Bx')))
qnames <-  c('EV','STD','Q.05','Q.95',sapply(0:1,function(x){paste0('theta',x)}))


gridoutcomes <- t(apply(expand.grid(0:1,0:1),1,rev))
outcomes <- apply(gridoutcomes,1,frombinary)
binoutcomes <- lapply(outcomes,function(x){tobinary(x,2)})
noutcomes <- length(outcomes)
alcombnames <- sapply(binoutcomes,function(x){do.call(paste0,as.list(alnames[x+1]))})
ncombinations <- (length(allgenes)^2-length(allgenes))/2

## if(cores>1){
## cl <- makeCluster(cores)
## registerDoParallel(cl)
## }
for(i in allsymptoms){

    results <- list()
    spread <- matrix(NA,nrow=ncombinations,ncol=3)
    colnames(spread) <- c('spread','gene1','gene2')
    rownames(spread) <- NULL
    ii <- 0
    ymin <- +Inf
    ymax <- -Inf
    for(g1 in allgenes[-length(allgenes)]) {
        for(g2 in allgenes[-c(1:g1)]) {
            ii <- ii+1
            result <- as.matrix(read.csv(paste0(savedir,filename,g1,'-',g2,'_s',c('A','B','C')[i],'.csv'))[,-1])
            results[[ii]] <- result
            spread[ii,] <- c(max(result[1,])-min(result[1,]) ,g1,g2)
            ymin <- min(ymin,c(result[1,]-result[2,]))
            ymax <- max(ymax,c(result[1,]+result[2,]))
        }}
    print(paste('ymin',ymin,'ymax',ymax))
    write.csv(spread,paste0('spread_2genes_s',c('A','B','C')[i],'.csv'))

}
## if(cores>1){stopCluster(cl)}

## A min max
## 0.207918316907101, 0.32303495937424
##
## B min max
## 0.335406572681691, 0.495964560003159
##
## A min max
## 0.207918316907101, 0.32303495937424
