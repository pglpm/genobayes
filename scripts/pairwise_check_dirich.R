source('mutualinfo_d1.R')

npairs <- 94*93/2

pb <- txtProgressBar(min = 1, max = 93,style=3)
res <- do.call(cbind,sapply(1:93,function(i){    setTxtProgressBar(pb, i)
    sapply((i+1):94,function(j){
        c(i,j,calculatemutualinfo(c(i,j),cores=10,messages=FALSE)[1])
    })}))
close(pb)


saveRDS(res,'results_pairs_d.rds')
write.table(res,paste0('results_pairs_d.csv'),sep=',',row.names=F,col.names=F,na='Infinity') 

