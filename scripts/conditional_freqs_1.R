## this script computes the probability for the conditional frequencies
## frequency(symptom | snp)
## where each symptom can assume values listed in 'symptomvariants'
## and each snp can assume values (alleles) listed in 'snpvariants'
## such values can be tuples
condfreqstatistics <- function(data,symptoms,symptomvariants,snps,snpvariants,namesymptoms,namesymptomvariants,namesnps,namesnpvariants,namestatistics,savedir,filename,logpriortheta,spread,writeflag=TRUE,cores=1){

numsymptoms <- length(symptoms)
numsymptomvariants <- length(symptomvariants)
numsnps <- length(snps)
numsnpvariants <- length(snpvariants)

if(writeflag){dir.create(savedir)}
    
## setup parallel processing
`%dox%` <- `%do%`
if(cores>1){
`%dox%` <- `%dopar%`
cl <- makeCluster(cores, outfile="")
registerDoParallel(cl)
}

result <- foreach(symptom=1:numsymptoms,
                  .export=c('savedir','filename','writeflag','sdata','symptoms','prefixsymptoms','namesymptoms','symptomvariants','numsymptomvariants','namesymptomvariants','snps','prefixsnps','namesnps','snpvariants','numsnpvariants','namesnpvariants','namestatistics','logpriortheta')
                  ) %:%
    foreach(snp=1:numsnps,
                  .export=c('savedir','filename','writeflag','sdata','symptoms','prefixsymptoms','namesymptoms','symptomvariants','numsymptomvariants','namesymptomvariants','snps','prefixsnps','namesnps','snpvariants','numsnpvariants','namesnpvariants','namestatistics','logpriortheta')
            ) %dox%
    {
        sdata <- data[,c(symptoms[[symptom]],snps[[snp]])] # load data for that symptom and snp

        ## calculate conditional frequencies
        ## one row per symptomvariants
        ## one column per snpvariant
        f <- sapply(snpvariants,function(snpvariant){
            sapply(symptomvariants,function(symptomvariants){
                sum(apply(sdata,1,function(z){all(z==c(symptomvariants,snpvariant))}))
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
                + sum(logpriortheta(lt)) )
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
        maxsearch <- nlm(f=logprob,p=rep(0,numsymptomvariants),iterlim=1e6 )
##        maxsearch <- nlm(f=logprob,p=log(apply(f,1,mean)/100),iterlim=1e6 )
        ## ## old method with optim() performed poorly
        ## maxsearch <- optim(par=log(apply(f,1,mean)/100),fn=logprob,control=list(maxit=1e6,reltol=1e-12),#,parscale=c(f[,1])),
        ##                 method="Nelder-Mead"
        ##                 #gr=gradient,method="BFGS"
        ##                 #method='L-BFGS-B',lower=c(1e-10,1e-10)
        ##                 )
        if(maxsearch$code>2){print(paste0('warn: ',maxsearch$code,' symptom=',symptom,' snp=',snp))}
        
        theta <- exp(maxsearch$estimate)
        fnew <- f+theta
        nnew <- apply(fnew,2,sum)
        quantities <- rbind(
            t(as.matrix(t(fnew)/nnew)), # EV
            ## STD
            t(sqrt(apply(fnew,1,function(x){x*(nnew-x)/(nnew^2*(1+nnew))}))),
            ## ## quantiles
            ## sapply(1:dim(fnew)[2],function(al){qbeta(quantiles,fnew[2,al],fnew[1,al])}),
            fnew, # posterior theta
            matrix(c(theta,rep(NA,numsymptomvariants*(numsnpvariants-1))),nrow=numsymptomvariants) # theta from optimization, only in first column
        )
        quantities <- rbind(quantities,
                            matrix(c(spread(quantities,numsymptomvariants,numsnpvariants),rep(NA,numsymptomvariants*(numsnpvariants-1))),nrow=numsymptomvariants)) # spreads, only in first column

    rownames(quantities) <- namestatistics
    colnames(quantities) <- namesnpvariants

        if(writeflag){
            write.csv(quantities,paste0(savedir,filename,prefixsymptoms,namesymptoms[symptom],prefixsnps,namesnps[snp],'.csv'))
        }        
    quantities
}
if(cores>1){
stopCluster(cl)
}
result}
