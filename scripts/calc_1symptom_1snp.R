## Calculation of expected value, std, thetas of conditional frequencies
## of symptom combos given SNP+allele
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

source('conditional_freqs_1.R')

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

savedir <- '1sympt_1snp/' # directory for saving results
filename <- 'freq_'
writeflag <- TRUE # whether to write the results for each case/snp combination in a file

cores <- 1 # for parallel processing

symptoms <- as.list(1:3)
namesymptoms <- c('A','B','C')
prefixsymptoms <- 'sym_' # for filename

#binarysymptomvariants <- c(0,1,2,4,3,5,6,7) # auxiliary quantity
symptomvariants <- 0:1
namesymptomvariants <- c('n','y')

snps <- as.list(3+(1:94))
namesnps <- 1:94 # colnames(data)[snps+3]
prefixsnps <- 'snp_' # for filename

snpvariants <- as.list(0:1)
namesnpvariants <- c('0','1')

namestatistics <- c(sapply(c('EV_','SD_','post.theta_','opt.theta_','max.spread_'),function(y){
    sapply(namesymptomvariants,function(x){paste0(y,x)})}))

## log-prior for thetas: see research notes
logpriortheta <- function(lt){dcauchy(lt,location=log(1000),scale=log(1000),log=TRUE)-lt}

## measure of spread, applied to the final matrix of quantities
spread <- function(q,numsymptomvariants,numsnpvariants){sapply(1:numsymptomvariants,function(co){max(abs(sapply(1:(numsnpvariants-1),function(x){sapply((x+1):numsnpvariants,function(y){(q[co,x]-q[co,y])/(q[co+numsymptomvariants,x]+q[co+numsymptomvariants,y])})})))})}


results <- condfreqstatistics(data,symptoms,symptomvariants,snps,snpvariants,namesymptoms,namesymptomvariants,namesnps,namesnpvariants,namestatistics,savedir,filename,logpriortheta,spread,writeflag=TRUE,cores=1)
