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

savedir <- 'testdirectory/' # directory for saving results
dir.create(savedir)
filename <- 'freq_'
writeflag <- TRUE # whether to write the results for each case/snp combination in a file

cores <- 1 # for parallel processing

## the script computes the probability for the conditional frequencies
## frequency(case | snp)
## where each case can assume values listed in 'combos'
## and each snp can assume values listed in 'alleles'
## such values can be tuples
cases <- list(1:3)
numcases <- length(cases)
namecases <- ''
prefixcases <- '' # for filename

binarycombos <- c(0,1,2,4,3,5,6,7) # auxiliary quantity
combos <- lapply(binarycombos,function(x){tobinary(x,3)})
numcombos <- length(combos)
namecombos <- c('0','A','B','C','AB','AC','BC','ABC')

snps <- as.list(3+(1:94))
numsnps <- length(snps)
namesnps <- 1:94 # colnames(data)[snps+3]
prefixsnps <- 'snp_' # for filename

alleles <- as.list(0:1)
numalleles <- length(alleles)
namealleles <- c('0','1')

namestatistics <- c(sapply(c('EV_','SD_','post.theta_','opt.theta_','max.spread_'),function(y){
    sapply(namecombos,function(x){paste0(y,x)})}))

## log-prior for thetas: see research notes
logpriortheta <- function(lt){dcauchy(lt,location=log(1000),scale=log(1000),log=TRUE)-lt}

## measure of spread, applied to the final matrix of quantities
spread <- function(q){sapply(1:numcombos,function(co){max(abs(sapply(1:(numalleles-1),function(x){sapply((x+1):numalleles,function(y){(q[co,x]-q[co,y])/(q[co+numcombos,x]+q[co+numcombos,y])})})))})}

## setup parallel processing
if(cores>1){
cl <- makeCluster(cores, outfile="")
registerDoParallel(cl)
}

results <- foreach(case=1:numcases,
                  .export=c('savedir','filename','writeflag','sdata','cases','prefixcases','namecases','combos','numcombos','namecombos','snps','prefixsnps','namesnps','alleles','numalleles','namealleles','namestatistics','logpriortheta')
                  ) %:%
    foreach(snp=1:numsnps,
                  .export=c('savedir','filename','writeflag','sdata','cases','prefixcases','namecases','combos','numcombos','namecombos','snps','prefixsnps','namesnps','alleles','numalleles','namealleles','namestatistics','logpriortheta')
            ) %do%
    {
        sdata <- data[,c(cases[[case]],snps[[snp]])] # load data for that case and snp

        ## calculate conditional frequencies
        ## one row per combo
        ## one column per allele
        f <- sapply(alleles,function(allele){
            sapply(combos,function(combo){
                sum(apply(sdata,1,function(z){all(z==c(combo,allele))}))
            })
        })

        ## minus log-probability for minimization
        ## we write thetas as exponentials of real numbers, to be positive
        ## see research notes
        logprob <- function(lt){
            t <- exp(lt)
            r2 <- f+t
            -(sum(lgamma(r2))
                - sum(lgamma(apply(r2,2,sum)))
                + numalleles * (lgamma(sum(t)) - sum(lgamma(t)))
                + sum(logpriortheta(lt)) )
        }
        ## gradient <- function(lt){
        ##     t <- exp(lt)
        ##     r2 <- f+t
        ##     -t*(apply(digamma(r2),1,sum)-
        ##       sum(digamma(apply(r2,2,sum))) +
        ##       noutcomes*(digamma(sum(t)) -
        ##          digamma(t)))
        ## }
        ## search parameter Theta with max evidence
        maxsearch <- nlm(f=logprob,p=log(apply(f,1,mean)),iterlim=1e6 )
        ## ## old method with optim() performed poorly
        ## maxsearch <- optim(par=log(apply(f,1,mean)/100),fn=logprob,control=list(maxit=1e6,reltol=1e-12),#,parscale=c(f[,1])),
        ##                 method="Nelder-Mead"
        ##                 #gr=gradient,method="BFGS"
        ##                 #method='L-BFGS-B',lower=c(1e-10,1e-10)
        ##                 )
        if(maxsearch$code>2){print(paste0('warn: ',maxsearch$code,' case=',case,' snp=',snp))}
        
        theta <- exp(maxsearch$estimate)
        fnew <- f+theta
        nnew <- apply(fnew,2,sum)
        quantities <- rbind(
            t(as.matrix(t(fnew)/nnew)), # EV
            ## STD
            t(sqrt(apply(fnew,1,function(x){x*(nnew-x)/(nnew^2*(1+nnew))}))),
            ## ## quantiles
            ## sapply(1:dim(fnew)[2],function(al){qbeta(quantiles,fnew[2,al],fnew[1,al])}),
            fnew, # posterior theta
            matrix(c(theta,rep(NA,numcombos*(numalleles-1))),nrow=numcombos) # theta from optimization, only in first column
        )
        quantities <- rbind(quantities,
                            matrix(c(spread(quantities),rep(NA,numcombos*(numalleles-1))),nrow=numcombos)) # spreads, only in first column

    rownames(quantities) <- namestatistics
    colnames(quantities) <- namealleles

        if(writeflag){
            write.csv(quantities,paste0(savedir,filename,prefixcases,namecases[case],prefixsnps,namesnps[snp],'.csv'))
        }
        
    quantities
}
if(cores>1){
stopCluster(cl)
}

