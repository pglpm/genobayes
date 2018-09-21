source('mutualinfo2rand.R')

npairs <- 94*93/2

pb <- txtProgressBar(min = 0, max = 94,style=3)
res <- sapply(1:94,function(i){setTxtProgressBar(pb, i)
        calculatemutualinfo(c(i),messages=FALSE)[1]})
close(pb)


saveRDS(res,'results_single.rds')
write.table(res,paste0('results_single.csv'),sep=',',row.names=F,col.names=F,na='Infinity') 
