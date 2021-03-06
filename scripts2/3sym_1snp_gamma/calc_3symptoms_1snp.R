## Calculation of expected value, std, thetas, spread of conditional frequencies
## of all symptom combos given each SNP

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

##auxiliary function to write the list of symptom combinations: 000, 100, ...,111
tobinary  <- function(number, noBits){as.numeric(intToBits(number))[1:noBits]}
#frombinary  <- function(digits) {sum(2^(0:(length(digits)-1))*digits)}

dpath  <-  "./"
datafile <- 'dataset1_binarized.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
data <- read.csv(paste0(dpath,nfile[1]))[,-1]
#n <- length(data[,1])

savedir <- '3sym_1snp_cauchy/' # directory for saving results
filename <- 'freq-3_1_cauchy' # filename prefix
writethreshold <- 2 # write the results in a file if spread>=writethreshold

cores <- 30 # for parallel processing

symptoms <- list(1:3) # all symptoms at once, they correspond to one index tuple: c(1,2,3)
namesymptoms <- '' # there's only one case
prefixsymptoms <- '' # for filename

binarysymptomvariants <- c(0,1,2,4,3,5,6,7) # auxiliary quantity
## below gives the list of all combinations: c(0,0,0), c(1,0,0), ... c(1,1,1)
symptomvariants <- lapply(binarysymptomvariants,function(x){tobinary(x,3)})
namesymptomvariants <- c('0','O','M','T','OM','OT','MT','OMT')

snps <- as.list(3+(1:94)) # list of gene indices in data
namesnps <- colnames(data)[(1:94)+3]
prefixsnps <- 'snp_' # for filename

snpvariants <- list(0, 1) # list of allele values
namesnpvariants <- c('A','B') # allele names


## log-prior for thetas: see research notes
## 'lt' is the log of theta
logpriortheta <- function(lt,t){sum(dcauchy(lt,location=log(1000),scale=log(1000),log=TRUE))-sum(lt)}
#logpriortheta <- function(lt,t){dgamma(sum(t),shape=1,scale=1000,log=TRUE)-7*log(sum(t))}

## measure of spread, applied to the final matrix of quantities
## it calculates abs((EV_freq1 - EV_freq2)/(SD_freq1 + SD_freq2))
## for all pairs and takes the maximum
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
