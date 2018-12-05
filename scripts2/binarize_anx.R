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

#0-7 (normal), 8-10 (mild), 11-14 (moderate), 15-21 (severe).
categories <- function(s){
    if(is.na(s)){return(NA)}
    if(0<=s && s<=7){return(0)}
    else if(8<=s && s<=10){return(1)}
    else if(11<=s && s<=14){return(2)}
    else if(15<=s && s<=21){return(3)}
    else {return(NA)}
    }

dpath  <-  "./"
datafile <- '_dataset1.csv'
nfile  <-  dir(path = dpath,pattern = datafile)
d <- read.csv(paste0(dpath,nfile[1]),sep=";")
predata <- d[#sample(1:dim(d)[1],4000)
       ,c(1,2,6)]
n <- dim(predata)[1]
v <- dim(predata)[2]

symp <- t(sapply(1:n,function(i){symptoms(predata[i,3])}))
categ <- sapply(1:n,function(i){categories(predata[i,2])})
data <- na.omit(cbind(predata[,1],categ,symp))

colnames(data) <- c('id','anxiety','O','M','T')

data <- data[order(data[,1]),]
write.csv(data,paste0(dpath,'anxiety_dataset1_binarized.csv'),row.names=F)
## res <- sapply(1:n,function(i){all(c(symptoms(predata[i,1]),predata[i,2:v])==data[i,])})

## all(res==TRUE)
