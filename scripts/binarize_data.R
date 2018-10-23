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
datafile <- 'dataset1.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
d <- read.csv(paste0(dpath,nfile[1]),sep=";")
predata <- d[#sample(1:dim(d)[1],4000)
       ,c(6,7+(1:94))]
predata[,1+(1:94)] <- (predata[,1+(1:94)]>0)*1
n <- dim(predata)[1]
v <- dim(predata)[2]

symp <- t(sapply(1:n,function(i){symptoms(predata[i,1])}))
data <- cbind(symp,predata[,2:v])
colnames(data) <- c('A','B','C',colnames(data)[-(1:3)])

write.table(data,paste0(dpath,'dataset1_binarized.csv'),sep=',',row.names=F,na='NA')
res <- sapply(1:n,function(i){all(c(symptoms(predata[i,1]),predata[i,2:v])==data[i,])})

all(res==TRUE)
