## Calculation of expected value, std, and quantiles of conditional frequencies
## of symptoms given gene+allele
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
library('GA')
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
data <- read.csv(paste0(dpath,nfile[1]))[,-1]
n <- length(data[,1])

filename <- 'allcondfreqs_1gene_thetamax' # where to save the results
allgenes <- 1:94
allsymptoms <- 1:3
quantiles <- c(0.05,0.95)
cores <- 1

## row an column headers for results
gnames <- colnames(data)[allgenes+3]
ganames <- c(rbind(paste0(gnames,'-AA'),paste0(gnames,'-Bx')))
qnames <-  c('EV','STD','Q.05','Q.95',sapply(0:1,function(x){paste0('theta',x)}))

## prior expected frequencies for each symptom (parameter alpha)

results <- list(A=NULL,B=NULL,C=NULL)

if(cores>1){
cl <- makeCluster(cores)
registerDoParallel(cl)
}
for(i in allsymptoms){
    sdata <- data[,c(i,3+allgenes)]
    result <- foreach(g=allgenes,
                     .combine=cbind, .export=c('sdata')) %do% {
                       f <- sapply(0:1,function(x){
                           table(sdata[sdata[[1+g]]==x,c(1,1+g)])
                       }) # one column per allele
                       ## functions for maximization
                       logprob <- function(t){
                           r2 <- f+t
                           -(sum(lgamma(r2))-
                             sum(lgamma(apply(r2,2,sum))) +
                             2*(lgamma(sum(t)) -
                                sum(lgamma(t))))
                       }
                       gradient <- function(t){
                           r2 <- f+t
                           -(apply(digamma(r2),1,sum)-
                             sum(digamma(apply(r2,2,sum))) +
                             2*(digamma(sum(t)) -
                                digamma(t)))
                       }
                       ## search parameter Theta with max evidence
                       maxsearch <- optim(par=c(1,1),fn=logprob,gr=gradient,control=list(maxit=1e6),
                                        method="Nelder-Mead"
                                        #method="BFGS"
                                        #method='L-BFGS-B',lower=c(1e-10,1e-10)
                                        )
                       fnew <- f+maxsearch$par
                       nnew <- apply(fnew,2,sum)
                       rbind(
                           t(as.matrix(t(fnew)[,-1]/nnew)), # EV
                           ## STD
                           sqrt(c(t(t(apply(fnew,2,prod))/((nnew^2)*(1+nnew))))),
                           ## quantiles
                           sapply(1:dim(fnew)[2],function(al){qbeta(quantiles,fnew[2,al],fnew[1,al])}),
                           matrix(rep(maxsearch$par,2),ncol=2) # Theta with max evidence
                       )
                       }
    rownames(result) <- qnames
    colnames(result) <- ganames

    write.csv(result,paste0(dpath,filename,'_s',c('A','B','C')[i],'.csv'))
}

if(cores>1){
stopCluster(cl)
}







stop('No error, just end of script')
    
    ## logprobold <- function(t0,t1){
    ##     sum(
    ##         apply(lgamma(result+c(t0,t1)),2,sum)-
    ##         lgamma(apply(result,2,sum)+t0+t1) +
    ##         lgamma(t0+t1) -
    ##         sum(lgamma(c(t0,t1)))
    ##     )
    ## }
## this is faster
gene <- 1
    logprob <- function(t){
        r2 <- result[,c(gene,gene+1)]+t
        -(sum(lgamma(r2))-
            sum(lgamma(apply(r2,2,sum))) +
            2*(lgamma(sum(t)) -
                            sum(lgamma(t))))
    }

    gradient <- function(t){
        r2 <- result[,c(gene,gene+1)]+t
        -(apply(digamma(r2),1,sum)-
            sum(digamma(apply(r2,2,sum))) +
            2*(digamma(sum(t)) -
                            digamma(t)))
    }
 maxres <-    ga(type='real-valued', fitness=logprob, lower=c(1e-10,1e-10),upper=c(1e12,1e12),popSize = 50, maxiter = 500, run = 100, parallel=FALSE,optim=FALSE,optimArgs=list(gr=gradient,control=list(maxit=2500)))

    
    maxpars <- optrobot@solution
    maxval <- -optrobot@fitnessValue



    
    optim(par=c(1,1),fn=logprob,gr=gradient,control=list(maxit=1e6),
                                        method="Nelder-Mead"
                                        #method="BFGS"
                                        #method='L-BFGS-B',lower=c(1e-10,1e-10)
          )

cl <- makeCluster(cores)
registerDoParallel(cl)

for(i in allsymptoms){
    sdata <- data[,c(i,3+allgenes)]
    sprior <- allsprior[,i]
    res <- t(foreach(g=allgenes,
                     .combine=cbind, .export=c('sdata','sprior','aa')) %dopar% {
                       sapply(0:1,function(x){
                           f <- table(sdata[sdata[[1+g]]==x,c(1,1+g)])
                           a2 <- sum(f)+aa
                           f2 <- (f+aa*sprior)/a2
                           sstd <- sqrt(prod(f2)/(1+a2))
                           squant <- qbeta(quantiles,a2*f2[2],a2*f2[1])
                           c(f2[2],sstd,squant,sum(f))
                       })
                     })
    rownames(res) <- rown
    colnames(res) <- coln

    write.csv(res,paste0(dpath,filename,'_s',c('A','B','C')[i],'.csv'))
    ##,sep=',',row.names=rown,col.names=coln,na='NA')
}

stopCluster(cl)
