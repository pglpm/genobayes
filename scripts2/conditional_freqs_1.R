## this script computes the probability for the conditional frequencies
## frequency(symptom | snp)
## where each symptom can assume values listed in 'symptomvariants'
## and each snp can assume values (alleles) listed in 'snpvariants'
## such values can be tuples
condfreqstatistics <- function(
                               data,# matrix of data
                               symptoms,# list of (groups of) indices in data
                               symptomvariants,# list of (combinations) of values for each symptom
                               snps,# list of (groups of) indices in data
                               snpvariants,# list of (combinations) of values for each symptom
                               namesymptoms,# names of each group of symptoms
                               namesymptomvariants,# names of each group of symptom values
                               namesnps,# names of each group of snps
                               namesnpvariants,# names of each group of snp values
                               savedir,# directory to save results
                               filename,# initial part of file name
                               logpriortheta,# prior log-probability of log(theta)
                               spread,# spread measure
                               writeflag=-Inf,# write results to files?
                               cores=1# cores for parallel processing, 1=no parallel
                               ){

    numsymptoms <- length(symptoms)
    numsymptomvariants <- length(symptomvariants)
    numsnps <- length(snps)
    numsnpvariants <- length(snpvariants)
    
    namestatistics <- c(sapply(c('EV_','SD_','post.theta_','opt.theta_','max.spread_'),function(y){
    sapply(namesymptomvariants,function(x){paste0(y,x)})}))

    if(writeflag==FALSE){writeflag <- +Inf}

    if(writeflag < +Inf){dir.create(savedir)}

## setup parallel processing
`%dox%` <- `%do%`
if(cores>1){
`%dox%` <- `%dopar%`
cl <- makeCluster(cores, outfile="")
registerDoParallel(cl)
}

result <- foreach(symptom=1:numsymptoms,
                  .export=c('savedir','filename','writeflag','data','symptoms','prefixsymptoms','namesymptoms','symptomvariants','numsymptomvariants','namesymptomvariants','snps','prefixsnps','namesnps','snpvariants','numsnpvariants','namesnpvariants','namestatistics','logpriortheta','spread')
                  ) %:%
    foreach(snp=1:numsnps,
                  .export=c('savedir','filename','writeflag','data','symptoms','prefixsymptoms','namesymptoms','symptomvariants','numsymptomvariants','namesymptomvariants','snps','prefixsnps','namesnps','snpvariants','numsnpvariants','namesnpvariants','namestatistics','logpriortheta','spread')
            ) %dox%
    {
        sdata <- data[,c(symptoms[[symptom]],snps[[snp]])] # load data for that symptom and snp

        ## calculate conditional frequencies
        ## one row per symptomvariants
        ## one column per snpvariant
        f <- sapply(snpvariants,function(snpvariant){
            sapply(symptomvariants,function(symptomvariant){
                sum(apply(sdata,1,function(z){all(z==c(symptomvariant,snpvariant))}))
            })
        })

        ## minus log-probability for minimization
        ## we write thetas as exponentials of real numbers, to be positive
        ## see research notes
        logprob <- function(lt){
            t <- exp(lt)
            r2 <- f+t
            -(sum(lgamma(r2))
                - sum(lgamma(apply(r2,2,sum)))
                + numsnpvariants * (lgamma(sum(t)) - sum(lgamma(t)))
                + logpriortheta(lt,t))
        }
        ## gradient <- function(lt){
        ##     t <- exp(lt)
        ##     r2 <- f+t
        ##     -t*(apply(digamma(r2),1,sum)-
        ##       sum(digamma(apply(r2,2,sum))) +
        ##       noutcomes*(digamma(sum(t)) -
        ##          digamma(t)))
        ## }
        ## search parameter Theta with max evidence
        mf <- apply(f,1,mean)
        maxsearch <- nlm(f=logprob,p=log(mf/sum(mf)),iterlim=1e6 )
        if(maxsearch$code>2){
            print(paste0('warning maximization: code=',maxsearch$code,' symptom=',symptom,' snp=',snp,' trying alterative...'))
            maxsearch <- nlm(f=logprob,p=rep(0,numsymptomvariants),iterlim=1e6 )
            }
        ## ## old method with optim() performed poorly
        ##                 method="Nelder-Mead"
        ##                 #gr=gradient,method="BFGS"
        ##                 #method='L-BFGS-B',lower=c(1e-10,1e-10)
        ##                 )
        if(maxsearch$code>2){print(paste0('failure maximization: code=',maxsearch$code,' symptom=',symptom,' snp=',snp))}
        
        theta <- exp(maxsearch$estimate)
        newtheta <- f+theta
        newA <- apply(newtheta,2,sum)
        quantities <- rbind(
            t(as.matrix(t(newtheta)/newA)), # EV
            ## STD
            t(sqrt(apply(newtheta,1,function(x){x*(newA-x)/(newA^2*(1+newA))}))),
            ## ## quantiles
            ## sapply(1:dim(newtheta)[2],function(al){qbeta(quantiles,newtheta[2,al],newtheta[1,al])}),
            newtheta, # posterior theta
            matrix(c(theta,rep(NA,numsymptomvariants*(numsnpvariants-1))),nrow=numsymptomvariants) # theta from optimization, only in first column
        )
        ## spreads
        spreads <- spread(quantities,numsymptomvariants,numsnpvariants)
        ## add spreads to matrix of quantities, only in first column
        quantities <- rbind(quantities,
                            matrix(c(spreads,rep(NA,numsymptomvariants*(numsnpvariants-1))),nrow=numsymptomvariants)
                            )

    rownames(quantities) <- namestatistics
    colnames(quantities) <- namesnpvariants

        if(max(spreads)>writeflag){
            write.csv(quantities,paste0(savedir,filename,prefixsymptoms,namesymptoms[symptom],'-',prefixsnps,namesnps[snp],'.csv'))
        }
    quantities
}
if(cores>1){
stopCluster(cl)
}
    print('writing results to file...')
    saveRDS(result,paste0(savedir,filename,'all.rds'))
    print('end')
    result}
