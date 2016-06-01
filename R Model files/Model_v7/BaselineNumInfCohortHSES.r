    # Most recent year number infected by birth cohort for HS and ES exit streams
    
    BaselineNumInfCohort=matrix(0,nCht,6)
    colnames(BaselineNumInfCohort)=c("HS detected [CI]","ES Detected [CI]","HS Infected [CI]","ES Infected [CI]","HS Missed [CI]","ES Missed [CI]")
    rownames(BaselineNumInfCohort)=c(chtStart:chtEnd)     

    eu1HS=formatC((baselineD$eHS)[,nTest],digits=2)
      eu1HS[(nCht-ageCapB[4]+1):nCht]=0
    eu1ES=formatC((baselineD$eES)[,nTest],digits=2)               #baseline test positives: ES (e.g. >72 months HS, etc...)
      eu1ES[(nCht-ageCapB[2]+1):nCht]=0
    eu1HSInf=formatC((truePrevD$eHSInf)[,nTest],digits=2)        #baseline test positives: HS (e.g. >72 months HS, etc...)
    eu1ESInf=formatC((truePrevD$eESInf)[,nTest],digits=2)        #baseline test positives: ES (e.g. >72 months HS, etc...)
    eu1HSMiss=round(as.numeric(eu1HSInf)-as.numeric(eu1HS) )     #number missed HS
          eu1HSMiss[eu1HSMiss<0]=0    #set negative values to zero - not missed any
    eu1ESMiss=round(as.numeric(eu1ESInf)-as.numeric(eu1ES) )#number missed ES
          eu1ESMiss[eu1ESMiss<0]=0   #set negative values to zero - not missed any
   # Confidence interval outputs
    eu1HSlci=formatC((qgamma(poL,as.numeric(eu1HS),rate=0,1)) ,digits=2)        #lower CI - baseline
    eu1HSuci=formatC((qgamma(poU,as.numeric(eu1HS),rate=0,1)),digits=2)         #upper CI - basleine
        eu1HSuci[as.numeric(eu1HSuci)<as.numeric(eu1HS)]=NA
    eu1ESlci=formatC((qgamma(poL,as.numeric(eu1ES),rate=0,1)),digits=2)         #lower CI - baseline
    eu1ESuci=formatC((qgamma(poU,as.numeric(eu1ES),rate=0,1)),digits=2)         #upper CI - basleine
        eu1ESuci[as.numeric(eu1ESuci)<as.numeric(eu1ES)]=NA

    eu1HSInflci=formatC((qgamma(poL,truePrevD$eHSInf,rate=0,1))[,nTest] ,digits=2)        #lower CI - baseline
    eu1HSInfuci=formatC((qgamma(poU,truePrevD$eHSInf,rate=0,1))[,nTest] ,digits=2)        #upper CI - basleine
                       eu1HSInfuci[as.numeric(eu1HSInfuci)<as.numeric(eu1HSInf)]=NA

    eu1ESInflci=formatC((qgamma(poL,truePrevD$eESInf,rate=0,1)) [,nTest] ,digits=2)       #lower CI - baseline
    eu1ESInfuci=formatC((qgamma(poU,truePrevD$eESInf,rate=0,1)) [,nTest] ,digits=2)       #upper CI - basleine
                             eu1ESInfuci[as.numeric(eu1ESInfuci)<as.numeric(eu1ESInf)]=NA
    eu1HSlciMiss=round(as.numeric(eu1HSInflci)-as.numeric(eu1HSlci)  )
          eu1HSlciMiss[eu1HSlciMiss<0]=0
    eu1HSuciMiss=round(as.numeric(eu1HSInfuci)-as.numeric(eu1HSuci) )
           eu1HSuciMiss[eu1HSuciMiss<0]=0
     eu1ESlciMiss=round(as.numeric(eu1ESInflci)-as.numeric(eu1ESlci)  )
        eu1ESlciMiss[eu1ESlciMiss<0]=0
    eu1ESuciMiss=round(as.numeric(eu1ESInfuci)-as.numeric(eu1ESuci) )
           eu1ESuciMiss[eu1ESuciMiss<0]=0
           
    BaselineNumInfCohort[,1]=paste(eu1HS,"[", eu1HSlci,",",eu1HSuci,"]")                 #baseline test positives: HS (e.g. >72 months HS, etc...)
    BaselineNumInfCohort[,2]=paste(eu1ES,"[", eu1ESlci,",",eu1ESuci,"]")                 #baseline test positives: ES (e.g. >72 months HS, etc...)
    BaselineNumInfCohort[,3]=paste(eu1HSInf,"[", eu1HSInflci,",",eu1HSInfuci,"]")        #baseline test positives: HS (e.g. >72 months HS, etc...)
    BaselineNumInfCohort[,4]=paste(eu1ESInf,"[", eu1ESInflci,",",eu1ESInfuci,"]")        #baseline test positives: ES (e.g. >72 months HS, etc...)
    BaselineNumInfCohort[,5]=paste(eu1HSMiss,"[", eu1HSlciMiss,",",eu1HSuciMiss,"]")     #number missed HS
    BaselineNumInfCohort[,6]=paste(eu1ESMiss,"[", eu1ESlciMiss,",",eu1ESuciMiss,"]")     #number missed ES
     
    BaselineNumInfCohort2<-BaselineNumInfCohort[8:dim(BaselineNumInfCohort)[1],]         # only output from 1992 to .csv file.
    
    write.csv(BaselineNumInfCohort2,paste(outDir,"/",MSName,"BaselineNumInfCohort.csv",sep=""))


