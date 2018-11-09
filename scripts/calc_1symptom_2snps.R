## Calculation of expected value, std, thetas, spread of conditional frequencies
## of each symptom given each pair of SNPs

#### PREAMBLE
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
#### END PREAMBLE

source('conditional_freqs_1.R') ## calls the function that calculates and writes the statistics

dpath  <-  "./"
datafile <- 'dataset1_binarized.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
data <- read.csv(paste0(dpath,nfile[1]))[,-1]
##n <- length(data[,1])

savedir <- '1sym_2snp/' # directory for saving results
filename <- 'freq-' # filename prefix
writeflag <- TRUE # whether to write the results for each case/snp combination in a file

cores <- 2 # for parallel processing

symptoms <- list(1,2,3) # symptoms A, B, C correspond to data indices 1, 2, 3
namesymptoms <- c('A','B','C')
prefixsymptoms <- 'sym_' # for filename

symptomvariants <- 0:1
namesymptomvariants <- c('n','y')

ngenes <- 94 # auxiliary quantity
snps <- unlist(lapply(1:(ngenes-1),function(x){lapply((x+1):ngenes,function(y){3+c(x,y)})}),
               recursive=FALSE) # list of all pairs of snps
namesnps <- sapply(snps,function(x){paste0(x[1]-3,'_',x[2]-3)})
prefixsnps <- 'snps_' # for filename

snpvariants <- list(c(0,0),c(0,1),c(1,0),c(1,1)) # list of allele-pair values
namesnpvariants <- c('00','01','10','11') # allele-pair names

## log-prior for thetas: see research notes
## 'lt' is the log of theta
logpriortheta <- function(lt){dcauchy(lt,location=log(1000),scale=log(1000),log=TRUE)-lt}

## measure of spread, applied to the final matrix of quantities
## it calculates abs((EV_freq1 - EV_freq2)/(SD_freq1 + SD_freq2))
## for all pairs and takes the maximum
spread <- function(q,numsymptomvariants,numsnpvariants){
    sapply(1:numsymptomvariants,function(co){
        max(abs(unlist(
            sapply(1:(numsnpvariants-1),function(x){
                sapply((x+1):numsnpvariants,function(y){
                    (q[co,x]-q[co,y])/(q[co+numsymptomvariants,x]+q[co+numsymptomvariants,y])
                })
            })
        )))
    })
}


results <- condfreqstatistics(data,symptoms,symptomvariants,snps,snpvariants,namesymptoms,namesymptomvariants,namesnps,namesnpvariants,savedir,filename,logpriortheta,spread,writeflag,cores)
