## Divide insomia values into symptoms and binarize gene data

## Daniela's explanation:
## 1-controls
## 2-A (DIFFICULTIES IN FALLING ASLEEP)
## 3-B (MULTIPLE AWAKENINGS DURING THE NIGHT)
## 4-C (EARLY MORNING AWAKENINGS)
## 5-AB
## 6-BC
## 7-AC
## 8-ABC

symptoms <- function(s){sym <- c(0,0,0)
    aa <- c(2,5,7,8)
    bb <- c(3,5,6,8)
    cc <- c(4,6,7,8)
    if(any(s==aa)) sym[1] <- 1
    if(any(s==bb)) sym[2] <- 1
    if(any(s==cc)) sym[3] <- 1
    sym}

dpath  <-  "./"
datafile <- '_replica_dataset.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
d <- read.csv(paste0(dpath,nfile[1]),sep=";")
predata <- d[#sample(1:dim(d)[1],4000)
       ,c(1,5,6:(dim(d)[2]-1))]
n <- dim(predata)[1]
v <- dim(predata)[2]
predata[,3:v] <- (predata[,3:v]>0)*1

symp <- t(sapply(1:n,function(i){symptoms(predata[i,2])}))
data <- cbind(predata[,1],symp,predata[,3:v])
colnames(data) <- c('id','O','M','T',colnames(data)[-(1:4)])

data <- data[order(data[,1]),]
write.csv(data,paste0(dpath,'replica_datasetb_binarized.csv'))
## res <- sapply(1:n,function(i){all(c(symptoms(predata[i,1]),predata[i,2:v])==data[i,])})

## all(res==TRUE)
