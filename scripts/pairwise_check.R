source('mutualinfo2.R')

npairs <- 94*93/2
res <- matrix(NA,npairs,3)

pb <- txtProgressBar(min = 0, max = npairs)
k <- 0
sapply(1:6,function(i){sapply((i+1):7,function(j){
    k <- k+1
    setTxtProgressBar(pb, k)
    res[k,] <- c(i,j,calculatemutualinfo(c(i,j),messages=FALSE)[1])
})})
close(pb)

saveRDS(res,'results_pairs.rds')
