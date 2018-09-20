source('mutualinfo2.R')

npairs <- 94*93/2

pb <- txtProgressBar(min = 0, max = 94)
res <- do.call(cbind,sapply(1:93,function(i){    setTxtProgressBar(pb, i)
    sapply((i+1):94,function(j){
        c(i,j,calculatemutualinfo(c(i,j),messages=FALSE)[1])
    })}))
close(pb)


saveRDS(res,'results_pairs.rds')
