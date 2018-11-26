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
                  .export=c('savedir','filename','data','symptoms','prefixsymptoms','namesymptoms','symptomvariants','numsymptomvariants','namesymptomvariants','snps','prefixsnps','namesnps','snpvariants','numsnpvariants','namesnpvariants','namesnpcombos','logpriortheta','statsfunction','mciterations','mcdiscard'),
                  .packages=c('LaplacesDemon')
                  ) %:%
    foreach(snp=1:numsnps,
                  .export=c('savedir','filename','data','symptoms','prefixsymptoms','namesymptoms','symptomvariants','numsymptomvariants','namesymptomvariants','snps','prefixsnps','namesnps','snpvariants','numsnpvariants','namesnpvariants','namesnpcombos','logpriortheta','statsfunction','mciterations','mcdiscard'),
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

        if(min(c(f))<2){
            print(paste0('warning, sym=',symptom,', snp=',snp,': ',min(c(f)),' sample'))
        }

        if(is.function(logpriortheta)){
### Case of general hierarchic prior

### Here is the set-up for the Markov-chain Monte Carlo sampling of the parameters
### Uses package LaplacesDemon

        PGF <- function(data){rep(0,numsymptomvariants)}

        mydata <- list(y=f, PGF=PGF,
               parm.names=paste0('u',(1:numsymptomvariants)-1),
               mon.names=c('')
               )

        logprob <- function(parm,data){
            parm <- interval(parm,-700,700)
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
        samples <- usamples$Posterior2
        if(length(samples)<mciterations*numsymptomvariants/2){
            print(paste0('warning, sym=',symptom,', snp=',snp,': low acceptance rate ', length(samples)/(mciterations*numsymptomvariants),'; using non-stationary samples'))
            samples <- usamples$Posterior1
        }
        }else{### case of constant prior. Mimicked by zero-valued samples
            samples <- matrix(0,nrow=2,ncol=numsymptomvariants)
            fnew <- f+1
            usamples <- list(
                Summary1=rep(NA,7),
                Summary2=rep(NA,7),
                LML= sum(lgamma(fnew)) - sum(lgamma(colSums(fnew)))
                )
        }
        
        endstats <- statsfunction(f,samples,numsymptomvariants,numsnpvariants)
        if(endstats$writeflag){
        write.csv(rbind(usamples$Summary1,c(usamples$LML,rep(NA,6)),usamples$Summary2),paste0(savedir,filename,'mcsummary-',prefixsymptoms,namesymptoms[symptom],'-',prefixsnps,namesnps[snp],'.csv'))
        }
        endstats
}
if(cores>1){
stopCluster(cl)
}
    print('writing results to file...')
    saveRDS(result,paste0(savedir,'_',filename,'all.rds'))
    print('end')
    result}





## This function produces samples to plot the initial distribution
priorsamples <- function(
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



### Here is the set-up for the Markov-chain Monte Carlo sampling of the parameters
### Uses package LaplacesDemon

        PGF <- function(data){rep(0,numsymptomvariants)}

        mydata <- list(y=1, PGF=PGF,
               parm.names=paste0('u',(1:numsymptomvariants)-1),
               mon.names=c('')
               )

        logprob <- function(parm,data){
            parm <- interval(parm,-700,700)
            t <- exp(parm)

            LP <- logpriortheta(parm,t) + sum(parm)

            list(LP=LP, Dev=-2*LP, Monitor=1, yhat=1, parm=parm)
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
        
        samples <- usamples$Posterior2
        if(length(samples)<mciterations*numsymptomvariants/2){
            print(paste0('warning, sym=',symptom,', snp=',snp,': low acceptance rate ', length(samples)/(mciterations*numsymptomvariants),'; using non-stationary samples'))
            samples <- usamples$Posterior1
        }

        print('writing samples to file...')
        write.csv(samples,paste0(savedir,filename,'priorsamples.csv'))
        saveRDS(samples,paste0(savedir,'_',filename,'priorsamples.rds'))

    print('producing frequency samples and writing them...')
    fsamples <- t(apply(samples,1,function(tt){
        et <- exp(tt)
        c(rbeta(1,et[1],et[2]),rbeta(1,et[1],et[2]))
        }))
        write.csv(fsamples,paste0(savedir,filename,'priorfsamples.csv'))
        saveRDS(fsamples,paste0(savedir,'_',filename,'priorfsamples.rds'))

    print('end')
}

