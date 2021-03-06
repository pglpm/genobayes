## Calculation of expected value, std, thetas, spread of conditional frequencies
## of each symptom given each SNP

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

savedir <- '1sym_1snp_cauchy/' # directory for saving results
filename <- 'freq-1_1_cauchy-' # filename prefix
writethreshold <- -Inf # write the results for each case/snp combination in a file when the spread is larger than this

cores <- 30 # for parallel processing

symptoms <- list(1,2,3) # symptoms A, B, C correspond to data indices 1, 2, 3
namesymptoms <- c('O','M','T')
prefixsymptoms <- 'sym_' # for filename

symptomvariants <- list(0,1)
namesymptomvariants <- c('absent','present')

snps <- as.list(3+(1:94)) # list of gene indices in data
namesnps <- sapply(snps,function(x){paste0(x[1]-3)})
## namesnps <- colnames(data)[(1:94)+3] # this was writing full SNP names
prefixsnps <- 'snp_' # for filename

snpvariants <- list(0, 1) # list of allele values
namesnpvariants <- c('A','B') # allele names

## log-prior for thetas: see research notes
## 'lt' is the log of theta

## first prior: the product of two Jeffreys priors for the two scale variables of the beta distribution = constant in log(variable), but regularized as a very broad Cauchy distribution. See research notes.
logpriortheta <- function(lt,t){sum(dcauchy(lt,location=log(1000),scale=log(1000),log=TRUE))-sum(lt)}

## second prior: constant in the frequency parameter and a very broad gamma density for the pseudocount parameter. See research notes.
#logpriortheta <- function(lt,t){dgamma(sum(t),shape=1,scale=1000,log=TRUE)-log(sum(t))}

## measure of spread, applied to the final matrix of quantities
## it calculates abs((EV_freq1 - EV_freq2)/sqrt(SD_freq1^2 + SD_freq2^2))
## for all pairs and takes the maximum.
## This corresponds to the mean/std of the distribution for the frequency difference
spread <- function(q,numsymptomvariants,numsnpvariants){
    sapply(1:numsymptomvariants,function(co){
        max(abs(unlist(
            sapply(1:(numsnpvariants-1),function(x){
                sapply((x+1):numsnpvariants,function(y){
                    (q[co,x]-q[co,y])/sqrt(q[co+numsymptomvariants,x]^2+q[co+numsymptomvariants,y]^2)
                })
            })
        )))
    })
}

results <- condfreqstatistics(data,symptoms,symptomvariants,snps,snpvariants,namesymptoms,namesymptomvariants,namesnps,namesnpvariants,savedir,filename,logpriortheta,spread,writethreshold,cores)
