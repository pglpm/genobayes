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

savedir <- './2gene-results_v2/'
filename <- 'genes' # where to save the results
allgenes <- 1:15
allsymptoms <- 1:3
quantiles <- c(0.05,0.95)
cores <- 25

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

if(cores>1){
cl <- makeCluster(cores)
registerDoParallel(cl)
}
for(i in allsymptoms){
    sdata <- data[,c(i,3+allgenes)]
    result <- foreach(g1=allgenes[-length(allgenes)], .export=c('sdata','binoutcomes','noutcomes')) %:%
        foreach(g2=allgenes[-c(1:(g1))],
                ##.combine=cbind,
                .export=c('sdata','binoutcomes','noutcomes')) %dopar% {
                    f <- sapply(binoutcomes,function(x){
                        table(factor(sdata[sdata[[1+g1]]==x[1]&sdata[[1+g2]]==x[2],1],levels=0:1))}) #one column per allele combination
                    ## functions for maximization
                    logprob <- function(lt){
                        t <- exp(lt)
                        r2 <- f+t
                        -(sum(lgamma(r2))-
                          sum(lgamma(apply(r2,2,sum))) +
                          noutcomes*(lgamma(sum(t)) -
                                     sum(lgamma(t))))
                    }
                    gradient <- function(lt){
                        t <- exp(lt)
                        r2 <- f+t
                        -t*(apply(digamma(r2),1,sum)-
                          sum(digamma(apply(r2,2,sum))) +
                          noutcomes*(digamma(sum(t)) -
                             digamma(t)))
                    }
                       ## search parameter Theta with max evidence
                       maxsearch <- optim(par=rep(0,2),fn=logprob,gr=gradient,control=list(maxit=1e6,reltol=1e-12),#,parscale=c(f[,1])),
                                        #method="Nelder-Mead"
                                        method="BFGS"
                                        #method='L-BFGS-B',lower=c(1e-10,1e-10)
                                        )
                    if(maxsearch$convergence>0){print(paste0('warn: ',maxsearch$convergence,' s',i,' g',g))}
                    theta <- exp(maxsearch$par)
                       fnew <- f+theta
                       nnew <- apply(fnew,2,sum)
                       quantities <- rbind(
                           t(as.matrix(t(fnew)[,-1]/nnew)), # EV
                           ## STD
                           sqrt(c(t(t(apply(fnew,2,prod))/((nnew^2)*(1+nnew))))),
                           ## quantiles
                           sapply(1:dim(fnew)[2],function(al){qbeta(quantiles,fnew[2,al],fnew[1,al])}),
                           matrix(rep(theta,noutcomes),nrow=2) # Theta with max evidence
                       )
                    rownames(quantities) <- qnames
                    colnames(quantities) <- alcombnames

                    write.csv(quantities,paste0(savedir,filename,g1,'-',g2,'_s',c('A','B','C')[i],'.csv'))
                    quantities
                }
    message('saving spread...')
    spread <- matrix(NA,nrow=ncombinations,ncol=3)
    colnames(spread) <- c('spread','gene1','gene2')
    rownames(spread) <- NULL
    ii <- 0
    ymin <- +Inf
    ymax <- -Inf
    for(g1 in allgenes[-length(allgenes)]) {
        for(g2 in allgenes[-c(1:g1)]) {
            ii <- ii+1
            resu <- as.matrix(read.csv(paste0(savedir,filename,g1,'-',g2,'_s',c('A','B','C')[i],'.csv'))[,-1])
            spread[ii,] <- c(max(resu[1,])-min(resu[1,]) ,g1,g2)
            ymin <- min(ymin,c(resu[1,]-resu[2,]))
            ymax <- max(ymax,c(resu[1,]+resu[2,]))
        }}
    print(paste('ymin',ymin,'ymax',ymax))
    write.csv(spread,paste0('spread_2genes_s',c('A','B','C')[i],'.csv'))

}

if(cores>1){
stopCluster(cl)
}
