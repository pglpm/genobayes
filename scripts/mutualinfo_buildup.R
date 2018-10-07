## Calculation of mutual information from sampled data
## for 1 gene, 2 genes, etc, keeping at each step those with highest mutual info

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

xlogy <- function(x,y){if(x==0){0}else{x*log(y)}}

dpath  <-  "./"
nfile  <-  dir(path = dpath,pattern = datafile)
d <- read.csv(paste0(dpath,nfile[1]),sep=";")

## this tuple contains the genes kept so far
keptgenes <- NA

## this tuple contains the genes to be checked in 
tocheck <- 1:94

changepar <- 2^7

for(iter=1:94){

    ## "if" to choose whether to parallelize outer or inner loop
    if(iter<changepar){# parallelize outer loop

        










    }







}



calculatemutualinfo <- function(whichgenes,a=2,cores=20,datafile='dataset1',messages=TRUE){

    if(any(whichgenes<1) | any(whichgenes>94)){print('Error: allowed genes within [1,94]')}
    ## load all data

    ## insomnia symptoms = col 6
    ## gene data = cols 8:102
    ## genes can be any subset of {1,...,95}
    genes <- c(whichgenes)
    ngenes <- length(genes)
    data <- d[#sample(1:dim(d)[1],4000)
       ,c(6,7+genes)]
    data[,2:dim(data)[2]] <- (data[,2:dim(data)[2]]>0)*1
    ## uncomment this to shuffle the symptom data
    ##    data[,1] <- sample(data[,1])
    l <- length(data[1,])
    n <- length(data[,1])
    nn <- 5.3e6L
    cs <- length(unique(data[,1]))
    cg <- 2^ngenes
    cx <- cs*cg

    if(a==0){aa  <- cx}else if(a==-2){aa <- sqrt(2*cx)}else if(a==-10){aa <- nn}else{aa <- a}

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
    cgplus <- dim(fg)[1]

    n2 <- n+aa
    dnna <- (nn-n)/(nn*n2)
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
        foreach(i=1:cgplus, .combine=cbind,
                .export=c('fg','fx','ngenes','cg','cs','cx','nn','n','n2','dnna','xlogy','aa')
                ) %dopar% {
            gam <- fg[i,]
            ## marginal freqs
            ns <- sig[2]
            ng <- gam[ngenes+1]
            ## frequency of pair
            nx <- sum(apply(fx,1,function(i){all(i[1:(ngenes+1)]==c(sig[1],gam[1:ngenes]))*i[ngenes+2]}))

            probx <- nx/nn+dnna*(nx+aa/cx)
            probs <- ns/nn+dnna*(ns+aa/cs)
            probg <- ng/nn+dnna*(ng+aa/cg)

            xlogy(probx,probx/(probs*probg))}})) -
        cgminus*2*dnna/cx*sum(log((fs[,2]/nn+dnna*(fs[,2]+aa/cs))*cx/cg))
    stopCluster(cl)
    
    
    sentropy <- -sum(apply(fs,1,function(sig){
        prob <- (sig[2]/nn+dnna*(sig[2]+aa/cs))
        xlogy(prob,prob)}))
if(messages==TRUE){
    print(paste0('max value = ',sentropy))
    print(paste0('mutual information = ',minfo))
    print(paste0('normalized value = ',minfo/sentropy))}
    c(minfo,sentropy)}

## all genes:

## [1] "max value = 1.73891986550377"
## [1] "mutual information = 1.73823079314495"
## [1] "normalized value = 0.999603735415018"
## There were 20 warnings (use warnings() to see them)

## "best genes" from pairwise analysis: {{32, 2}, {66, 2}, {6, 4}, {57, 1}, {75, 3}, {9, 1}, {58, 1}, {89, 1}, {53, 1}, {21, 1}, {56, 1}, {42, 1}, {68, 1}}
