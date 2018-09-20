source('mutualinfo2ptest.R')
res <- calculatemutualinfo(1:94,cores=30)
write.table(res,paste0('res_mutualinfo94g.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
