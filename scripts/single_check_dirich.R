source('mutualinfo_d1.R')

npairs <- 94*93/2

pb <- txtProgressBar(min = 1, max = 94,style=3)
res <- sapply(1:94,function(i){setTxtProgressBar(pb, i)
        calculatemutualinfo(c(i),cores=2,messages=FALSE)[1]})
close(pb)


saveRDS(res,'results_single_d.rds')
write.table(res,paste0('results_single_d.csv'),sep=',',row.names=F,col.names=F,na='Infinity') 
