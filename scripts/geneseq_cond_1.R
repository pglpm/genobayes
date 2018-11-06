## Calculation of beliefs for conditional frequencies
## for 1 gene, 2 genes, etc, keeping at each step those with highest spread

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

##xlogy <- function(x,y){if(x==0){0}else{x*log(y)}}

dpath  <-  "./"
datafile <- 'dataset1_binarized.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
data <- read.csv(paste0(dpath,nfile[1]))[,-1]
n <- length(data[,1])
symstrings <- c('A','B','C')

symptom <- 1
savedir <- paste0('cumulative_condfreq_s',symstrings[symptom],'/')

## max number of iterations
totalgenes <- 94
## this tuple contains the genes kept so far
keptgenes <- NULL
keptinfo <- NULL
allgenes <- 1:94
## this tuple contains the genes to be checked 
genegroup <- allgenes

cores <- 1

alnames <- c('A','B')
gnames <- colnames(data)[allgenes+3]

rname <- 'ac_'
warnstring <- 'warn1e3_'
cores <- 1
if(cores>1){cl <- makeCluster(cores)
registerDoParallel(cl)}

## iteration at which to switch parallel loop
changepar <- totalgenes - cores

for(ngenes in 1:totalgenes){
    message(paste0('iteration ',ngenes))
    ## "if" to choose whether to parallelize outer or inner loop
    if(ngenes<changepar){# parallelize outer loop

        allminfo <- foreach(i=genegroup, .combine=rbind, .export=c('xlogy','keptgenes','ngenes','n','cs','fs'),
                            .packages=c('dplyr','doParallel')) %dopar% {
            genes <- sort(c(keptgenes,i))
            data <- predata[,c(1,1+genes)]
            ## uncomment this to shuffle the symptom data
            ##    data[,1] <- sample(data[,1])
            cg <- 2^ngenes
            cx <- cs*cg
            
            ## data frequencies for genes
            fg <- data.matrix(tally(group_by_at(data,.vars=c(2:(1+ngenes)))))
            ## data frequencies for both
            fx <- data.matrix(tally(group_by_at(data,.vars=c(1:(ngenes+1)))))

            et <- dim(fx)[1]
            logpost <- function(x){lgamma(x)-lgamma(n+x)+sum(lgamma(fx[,ngenes+2]+x/cx))-et*lgamma(x/cx)-log(x)}
            ##aa <- optimize(f=logpost,interval=c(1,cx),maximum=T)$maximum
            ##if(aa==cx){write.table(c(ngenes,i),paste0(savedir,warnstring,ngenes,'_',i,'.csv'),sep=',',row.names=F,col.names=F)}
            aa <- cx ##print(aa)
            n2 <- n+aa
            dnna <- 1/n2 #(nn-n)/(nn*n2)
            
            cgminus <- cg-dim(fg)[1]
            cgplus <- dim(fg)[1]

            minfo <- sum(apply(fs,1,function(sig){
                ##  if(messages==TRUE){print(sig[1])}
                foreach(i=1:cgplus, .combine=cbind#,
                                        #  .export=c('fg','fx','ngenes','cg','cs','cx','n','n2','dnna','xlogy','aa')
                        ) %do% {
                                        gam <- fg[i,]
                                        ## marginal freqs
                                        ns <- sig[2]
                                        ng <- gam[ngenes+1]
                                        ## frequency of pair
                                        nx <- sum(apply(fx,1,function(i){all(i[1:(ngenes+1)]==c(sig[1],gam[1:ngenes]))*i[ngenes+2]}))
                                        
                                        probx <- dnna*(nx+aa/cx)
                                        probs <- dnna*(ns+aa/cs)
                                        probg <- dnna*(ng+aa/cg)

                                        xlogy(probx,probx/(probs*probg))}
            })) -
    cgminus*2*dnna/cx*sum(log((dnna*(fs[,2]+aa/cs))*cx/cg))
            c(minfo,aa)
        }
    } else {# parallelize inner loop
        allminfo <- foreach(i=genegroup, .combine=rbind #, .export=c('xlogy','keptgenes','ngenes','n','cs','fs','aa','n2','dnna'),
                           ) %do% {
            genes <- sort(c(keptgenes,i))
            data <- predata[,c(1,1+genes)]
            ## uncomment this to shuffle the symptom data
            ##    data[,1] <- sample(data[,1])
            cg <- 2^ngenes
            cx <- cs*cg

            ## data frequencies for genes
            fg <- data.matrix(tally(group_by_at(data,.vars=c(2:(1+ngenes)))))
            ## data frequencies for both
            fx <- data.matrix(tally(group_by_at(data,.vars=c(1:(ngenes+1)))))

            et <- dim(fx)[1]
            logpost <- function(x){lgamma(x)-lgamma(n+x)+sum(lgamma(fx[,ngenes+2]+x/cx))-et*lgamma(x/cx)-log(x)}
            ##aa <- optimize(f=logpost,interval=c(1,cx),maximum=T)$maximum
            ##if(aa==cx){write.table(c(ngenes,i),paste0(savedir,warnstring,ngenes,'_',i,'.csv'),sep=',',row.names=F,col.names=F)}
            aa <- cx
            n2 <- n+aa
            dnna <- 1/n2 #(nn-n)/(nn*n2)

            cgminus <- cg-dim(fg)[1]
            cgplus <- dim(fg)[1]

            minfo <- sum(apply(fs,1,function(sig){
                ##  if(messages==TRUE){print(sig[1])}
                foreach(i=1:cgplus, .combine=cbind , .export=c('fg','fx','ngenes','cg','cs','cx','n','n2','dnna','xlogy','aa')
                        ) %dopar% {
                                        gam <- fg[i,]
                                        ## marginal freqs
                                        ns <- sig[2]
                                        ng <- gam[ngenes+1]
                                        ## frequency of pair
                                        nx <- sum(apply(fx,1,function(i){all(i[1:(ngenes+1)]==c(sig[1],gam[1:ngenes]))*i[ngenes+2]}))
                                        
                                        probx <- dnna*(nx+aa/cx)
                                        probs <- dnna*(ns+aa/cs)
                                        probg <- dnna*(ng+aa/cg)

                                        xlogy(probx,probx/(probs*probg))}
            })) -
    cgminus*2*dnna/cx*sum(log((dnna*(fs[,2]+aa/cs))*cx/cg))
            c(minfo,aa)
        }
    }
    
        
    write.table(cbind(genegroup,allminfo),paste0(savedir,'iteration_',rname,ngenes,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
    ##print(allminfo)
    maxinfo <- max(allminfo[,1])
    keptinfo <- c(keptinfo,maxinfo)
    maxgene <- genegroup[which.max(allminfo[,1])]
    #maxaa <- allminfo[which.max(allminfo,1),2]
    keptgenes <- c(keptgenes,maxgene)
    genegroup <- allgenes[-keptgenes]
    write.table(cbind(keptgenes,keptinfo),paste0(savedir,'infoseq_',rname,ngenes,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
}
stopCluster(cl)

