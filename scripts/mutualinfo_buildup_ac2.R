## Calculation of mutual information from sampled data
## for 1 gene, 2 genes, etc, keeping at each step those with highest mutual info
## parameter A approx chosen via Bayes's theorem, optimized

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

xlogy <- function(x,y){if(x==0){0}else{x*log(y)}}

dpath  <-  "./"
datafile <- 'dataset1.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
d <- read.csv(paste0(dpath,nfile[1]),sep=";")
predata <- d[#sample(1:dim(d)[1],4000)
       ,c(6,7+(1:94))]
predata[,1+(1:94)] <- (predata[,1+(1:94)]>0)*1
n <- length(predata[,1])
cs <- length(unique(predata[,1]))
fs <- data.matrix(tally(group_by_at(predata,.vars=c(1))))
nn <- 5.3e6L

## iteration at which to switch parallel loop
changepar <- 7

totalgenes <- 1:94
## this tuple contains the genes kept so far
keptgenes <- NULL
keptinfo <- NULL
## this tuple contains the genes to be checked in 
genegroup <- totalgenes

cores <- 30
cl <- makeCluster(cores)
registerDoParallel(cl)

for(ngenes in totalgenes){
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
            aa <- 2^(0.7+0.8*(ngenes+3))
            n2 <- n+aa
            dnna <- 1/n2 #(nn-n)/(nn*n2)

            ## data frequencies for genes
            fg <- data.matrix(tally(group_by_at(data,.vars=c(2:(1+ngenes)))))
            ## data frequencies for both
            fx <- data.matrix(tally(group_by_at(data,.vars=c(1:(ngenes+1)))))

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
            minfo
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
            aa <- 2^(0.7+0.8*(ngenes+3))
            n2 <- n+aa
            dnna <- 1/n2 #(nn-n)/(nn*n2)
            ## data frequencies for genes
            fg <- data.matrix(tally(group_by_at(data,.vars=c(2:(1+ngenes)))))
            ## data frequencies for both
            fx <- data.matrix(tally(group_by_at(data,.vars=c(1:(ngenes+1)))))

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
            minfo
        }
    }
    
        
    write.table(cbind(genegroup,allminfo),paste0(savedir,'iteration_afromcopt_',ngenes,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
    maxinfo <- max(allminfo)
    keptinfo <- c(keptinfo,maxinfo)
    maxgene <- genegroup[which.max(allminfo)]
    keptgenes <- c(keptgenes,maxgene)
    genegroup <- totalgenes[-keptgenes]
    write.table(cbind(keptgenes,keptinfo),paste0(savedir,'infoseq_afromcopt_',ngenes,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
}
stopCluster(cl)

