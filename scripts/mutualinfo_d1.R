## Calculation of mutual information from sampled data
## test 1

## libraries and colour-blind palette from http://www.sron.nl/~pault/
library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('png')
library('plot3D')
library('doParallel')
library('GA')
library('dplyr')

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
plotsdir <- './'


## The script calculates the mutual information as in the genobayes.pdf notes, by calculating the frequencies for the combined variate and for the two variates separately in the data. Here we do this simply by counting. We choose a subset of genetic variations to check for calculation speed


## this function takes the ID of the genes we want to check, between 1 and 95, and returns the mutual info, its max value, and its normalized value
calculatemutualinfo <- function(whichgenes,cores=20,datafile='dataset1',messages=TRUE){

    if(any(whichgenes<1) | any(whichgenes>94)){print('Error: allowed genes within [1,94]')}
    ## load all data
    dpath  <-  "./"
    nfile  <-  dir(path = dpath,pattern = datafile)
    d <- read.csv(paste0(dpath,nfile[1]),sep=";")

    ## insomnia symptoms = col 6
    ## gene data = cols 8:102
    ## genes can be any subset of {1,...,95}
    genes <- c(whichgenes)
    ngenes <- length(genes)
    data <- d[#sample(1:dim(d)[1],4000)
       ,c(6,7+genes)]
    ## uncomment this to shuffle the symptom data
    ##    data[,1] <- sample(data[,1])
    l <- length(data[1,])
    n <- length(data[,1])
    #nn <- 5.3e6L
    cs <- length(unique(data[,1]))
    cg <- 2^ngenes
    cx <- cs*cg

    ## data frequencies for symptoms
    fs <- data.matrix(tally(group_by_at(data,.vars=c(1))))
    #print(dim(fs))
    ## data frequencies for genes
    fg <- data.matrix(tally(group_by_at(data,.vars=c(2:(1+ngenes)))))
    #print(dim(fg))
    ## data frequencies for both
    fx <- data.matrix(tally(group_by_at(data,.vars=c(1:(ngenes+1)))))
    #print(dim(fx))

    cgminus <- cg-dim(fg)[1]

    n2 <- n+2
    ## minfo <- log(ncnn) + sum(apply(fs,1,function(sig){
    ##     apply(fg,1,function(gam){
    ##         ## marginal freqs
    ##         ns <- sig[2]
    ##         ng <- gam[ngenes+1]
    ##         ## frequency of pair
    ##         nx <- sum(apply(fx,1,function(i){all(i[1:(ngenes+1)]==c(sig[1],gam[1:ngenes]))*i[ngenes+2]}))
    ##         numer <- nnc*nx+dn
            
    ##         numer/ncnn*(log(numer)-log(nnc*ns+dn*cg)-log(nnc*ng+dn*cs))})})) -
    ##     cgminus*dn/ncnn*sum(log(cs*(nnc*fs[,2]+dn*cg)))

        ## minfo2 <- log(ncnn) + sum(apply(fs,1,function(sig){
        ## apply(fg,1,function(gam){
        ##     ## marginal freqs
        ##     ns <- sig[2]
        ##     ng <- gam[ngenes+1]
        ##     ## frequency of pair
        ##     nx <- sum(apply(fx,1,function(i){all(i[1:(ngenes+1)]==c(sig[1],gam[1:ngenes]))*i[ngenes+2]}))
        ##     numer <- nnc*nx+dn
            
        ##     numer/ncnn*log(numer/((nnc*ns+dn*cg)*(nnc*ng+dn*cs)))})})) -
        ## cgminus*dn/ncnn*sum(log(cs*(nnc*fs[,2]+dn*cg)))

    cl <- makeCluster(cores)
    registerDoParallel(cl)
    minfo <- sum(apply(fs,1,function(sig){
if(messages==TRUE){print(sig[1])}
        foreach(i=1:(cg-cgminus), .combine=cbind,
                .export=c('fg','fx','sig','ngenes','cg','cs','cx','n','n2')
                ) %dopar% {
            gam <- fg[i,]
            ## marginal freqs
            ns <- sig[2]
            ng <- gam[ngenes+1]
            ## frequency of pair
            nx <- sum(apply(fx,1,function(i){all(i[1:(ngenes+1)]==c(sig[1],gam[1:ngenes]))*i[ngenes+2]}))

            probx <- (nx+2/cx)
            probs <- (ns+2/cs)
            probg <- (ng+2/cg)

            probx/n2*log(probx*n2/(probs*probg))}})) -
        cgminus*2/(n2*cx)*sum(log((fs[,2]+2/cs)*cx/cg/n2))
    stopCluster(cl)
    
    
    sentropy <- -sum(apply(fs,1,function(sig){
        prob <- (sig[2]+2/cs)/n2
        prob*log(prob)}))
if(messages==TRUE){
    print(paste0('max value = ',sentropy))
    print(paste0('mutual information = ',minfo))
    print(paste0('normalized value = ',minfo/sentropy))}
    c(minfo,sentropy)}

## with genes 6, 9, 21, 31, 32, 34, 46, 53, 57, 58, 59, 61, 64, 66, 68, 69, 72, 88, 89, 90 : info =  1.60235058302968 (92%)

## with 59, 68, 88, 89, 90 : info = 0.0704646052867898 (4%)

## with 6, 32, 57, 59, 64, 68, 72, 88, 89, 90 : info = 0.598136429089636 (34%)

## with 1:10 : info = 0.0305630433841963 (1.8%)

