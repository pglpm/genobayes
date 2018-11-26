## this script computes the probability for the conditional frequencies
## frequency(symptom | snp)
## where each symptom can assume values listed in 'symptomvariants'
## and each snp can assume values (alleles) listed in 'snpvariants'
## such values can be tuples
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
        saveRDS(result,paste0(savedir,'_',filename,'priorsamples.rds'))
        print('end')
}
