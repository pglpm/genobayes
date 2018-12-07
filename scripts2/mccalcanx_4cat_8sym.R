## Calculation of expected value, std, thetas, spread of conditional frequencies
## of each symptom given each SNP - extra data set

#### PREAMBLE
## libraries and colour-blind palette from http://www.sron.nl/~pault/
##runfunction <- function(aa=1e3L){
library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('png')
library('plot3D')
library('doParallel')
library('LaplacesDemon')
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

source('conditional_freqs_2.R') ## calls the function that calculates and writes the statistics

dpath  <-  "./"
datafile <- 'anxiety_dataset1_binarized.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
data <- read.csv(paste0(dpath,nfile[1]))[,-1]
##n <- length(data[,1])

savedir <- 'anxmc1e6_4cat_8sym_cons/' # directory for saving results
filename <- 'freq-4_8_cons-' # filename prefix

cores <- 30 # for parallel processing
mciterations <- 1e6 # number of Monte-Carlo samples

symptoms <- list(1) # anxiety category corresponds to data index 1
namesymptoms <- c('anxiety')
prefixsymptoms <- 'cat_' # for filename

symptomvariants <- as.list(0:3)
namesymptomvariants <- c(0:3)

snps <- list(c(2,3,4))
namesnps <- 'combo'
prefixsnps <- 'sym_' # for filename

snpvariants <- list(c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1),c(1,1,0),c(1,0,1),c(0,1,1),c(1,1,1)) # list of combo values
namesnpvariants <- c('N','O','M','T','OM','OT','MT','OMT') # allele-pair names

## combos of symptom values names (for conditional-frequency differences)
numsnpvariants <- length(snpvariants)
namesnpcombos <- unlist(sapply(1:(numsnpvariants-1),
                               function(y){sapply((y+1):numsnpvariants,
                                                  function(z){paste0(namesnpvariants[y],'-',namesnpvariants[z])}
                                                  )}
                               ))

## log-prior for thetas: see research notes
## 'lt' is the log of t

## first prior: constant in the frequency parameter and a very broad gamma density for the pseudocount parameter. See research notes.
logpriortheta <- function(lt,t){dgamma(sum(t),shape=1,scale=1000,log=TRUE)-3*log(sum(t))}

## second prior: constant. See research notes.
#logpriortheta <- FALSE

## third prior: constant in the frequency parameter and Jeffreys prior for the pseudocount variables of the beta distribution = constant in log(variable), but regularized as a very broad Cauchy distribution. See research notes.
#logpriortheta <- function(lt,t){dcauchy(log(sum(t)),location=log(1000),scale=log(1000),log=TRUE)-2*log(sum(t))}

## This function calculates the EVs and SDs of the marginals and all the differences of the conditional frequencies. For the latter also calculates the spreads ("significance"). It writes the results on two files if the max spread is larger than a given value

statsfunction <- function(f,samples,symptom,snp,numsymptomvariants,numsnpvariants,namesnpvariants=''){
    mval <- matrix(0:3,8,4,byrow=T)
    kam <- function(a,b){sum(abs(cumsum(a-b)))}

    diffsamples <- t(apply(samples,1,function(tt){
        fnew <- f+exp(tt)
        dd <- t(rdirichlet(numsnpvariants,t(fnew)))
        me <- colSums((0:3)*dd)
        va <- colSums(t(mval-me)^2*dd)
        sk <- colSums(t(mval-me)^3*dd)
        kdiffs <- sapply(2:8,function(z){kam(dd[,z],dd[,1])})
        
        c(me[-1]-me[1],
          va[-1]-va[1],
          sk[-1]-sk[1],
          kdiffs)
    }))

    colnames(diffsamples) <- c(paste0('diffmean_',namesnpvariants[-1]),
                               paste0('diffvar_',namesnpvariants[-1]),
                               paste0('diffskew_',namesnpvariants[-1]),
                               paste0('kantorovich_',namesnpvariants[-1]))

    writeflag <- TRUE
        write.csv(diffdata,paste0(savedir,filename,'spreads-',prefixsymptoms,namesymptoms[symptom],'-',prefixsnps,namesnps[snp],'.csv'))

    list(writeflag=writeflag)
    }

results <- condfreqstatistics(data,symptoms,symptomvariants,snps,snpvariants,namesymptoms,namesymptomvariants,namesnps,namesnpvariants,namesnpcombos,savedir,filename,logpriortheta,statsfunction,cores,mciterations)

## Call a function that samples without data, to plot the initial belief
## Not called in the case of constant prior
#
## if(is.function(logpriortheta)){
##     insamples <- priorsamples(data,symptoms,symptomvariants,snps,snpvariants,namesymptoms,namesymptomvariants,namesnps,namesnpvariants,namesnpcombos,savedir,filename,logpriortheta,statsfunction,cores,mciterations)
## }
