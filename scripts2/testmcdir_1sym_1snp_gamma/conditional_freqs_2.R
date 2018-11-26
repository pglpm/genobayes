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
                               namesnpcombos,# names of snp allele combinations
                               savedir,# directory to save results
                               filename,# initial part of file name
                               logpriortheta,# prior log-probability of log(theta)
                               statsfunction,# calculate relevant quantities & write to file
                               cores=1,# cores for parallel processing, 1=no parallel
                               mciterations=1e5,# number of Monte Carlo iterations
                               mcdiscard=500# number of initial Monte Carlo iterations to be discarded
                               ){

    numsymptoms <- length(symptoms)
    numsymptomvariants <- length(symptomvariants)
    numsnps <- length(snps)
    numsnpvariants <- length(snpvariants)

        dir.create(savedir)

    
## setup parallel processing
`%dox%` <- `%do%`
if(cores>1){
`%dox%` <- `%dopar%`
cl <- makeCluster(cores, outfile="")
registerDoParallel(cl)
}

result <- foreach(symptom=1:numsymptoms,
                  .export=c('savedir','filename','writethreshold','data','symptoms','prefixsymptoms','namesymptoms','symptomvariants','numsymptomvariants','namesymptomvariants','snps','prefixsnps','namesnps','snpvariants','numsnpvariants','namesnpvariants','namesnpcombos','logpriortheta','statsfunction','mciterations','mcdiscard'),
                  .packages=c('LaplacesDemon')
                  ) %:%
    foreach(snp=1:numsnps,
                  .export=c('savedir','filename','writethreshold','data','symptoms','prefixsymptoms','namesymptoms','symptomvariants','numsymptomvariants','namesymptomvariants','snps','prefixsnps','namesnps','snpvariants','numsnpvariants','namesnpvariants','namesnpcombos','logpriortheta','statsfunction','mciterations','mcdiscard'),
                  .packages=c('LaplacesDemon')
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

### Here is the set-up for the Markov-chain Monte Carlo sampling of the parameters
### Uses package LaplacesDemon

        PGF <- function(data){rep(0,numsymptomvariants)}

        mydata <- list(y=f, PGF=PGF,
               parm.names=paste0('u',(1:numsymptomvariants)-1),
               mon.names=c('')
               )

        logprob <- function(parm,data){
            t <- exp(parm)
            fnew <- data$y+t

            LL <- (sum(lgamma(fnew))
            - sum(lgamma(colSums(fnew)))
                + numsnpvariants * (lgamma(sum(t)) - sum(lgamma(t))))
            
            LP <- LL + logpriortheta(parm,t) + sum(parm)

            list(LP=LP, Dev=-2*LL, Monitor=1, yhat=1, parm=parm)
}

        ##Initial.Values <- c(0,0) #log(mf/sum(mf)) #GIV(prob, mydata, n=1000, PGF=T)
        
        Initial.Values <- GIV(logprob, mydata, n=1000, PGF=F)

        usamples <- LaplacesDemon(logprob, mydata, Initial.Values,
                        Covar=NULL,
                        Thinning=1,
                        Iterations=mciterations, Status=mciterations,
                        #Chains=nchains,CPUs=nchains,LogFile=paste0(filename,'_LDlog'), #Packages=c('Matrix'),#Type="MPI",
                        ##Algorithm="RDMH"#, Specs=list(B=list(1:d,d1:d2,d3:dnp))
                        ##Algorithm="Slice", Specs=list(B=NULL, Bounds=c(0,1), m=Inf, Type="Continuous", w=0.001)
                        Algorithm="AFSS", Specs=list(A=mcdiscard, B=NULL, m=100, n=0, w=1)
                        )

        statsfunction(f,usamples$Posterior2,numsymptomvariants,numsnpvariants)
}
if(cores>1){
stopCluster(cl)
}
    print('writing results to file...')
    saveRDS(result,paste0(savedir,'_',filename,'all.rds'))
    print('end')
    result}
