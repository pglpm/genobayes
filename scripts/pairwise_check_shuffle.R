source('mutualinfo2shuffle.R')

npairs <- 94*93/2

pb <- txtProgressBar(min = 0, max = 94,style=3)
res <- do.call(cbind,sapply(1:93,function(i){    setTxtProgressBar(pb, i)
    sapply((i+1):94,function(j){
        c(i,j,calculatemutualinfo(c(i,j),messages=FALSE)[1])
    })}))
close(pb)


saveRDS(res,'results_pairs2_shuffle.rds')


## c(6, 9, 17, 21, 30, 31, 32, 34, 38, 42, 46, 51, 53, 57, 58, 59, 61, 64, 66, 68, 69, 70, 72, 75, 78, 83, 88, 89, 90, 92, 93) give  0.00234867116380749

## c(6, 32, 57, 59, 64, 68, 72, 88, 89, 90) give  0.0791343853227863 (4% normalized)

## c(59, 68, 88, 89, 90) give
