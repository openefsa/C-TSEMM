    #All years data for number of infected animals by birth cohort for HS and ES exit streams
    
    eu1HS=formatC((HSbyChtPos1203),digits=2)
      ## INPUT colnames and rownames....
      eu1HS[(nCht-ageCapB[4]+1):nCht]=0
    eu1ES=formatC((ESbyChtPos1203),digits=2)                  #baseline test positives: ES (e.g. >72 months HS, etc...)
      eu1ES[(nCht-ageCapB[2]+1):nCht]=0
    eu1HSInf=formatC((truePrevD$eHSInf),digits=2)             #baseline test positives: HS (e.g. >72 months HS, etc...)
    eu1ESInf=formatC((truePrevD$eESInf),digits=2)             #baseline test positives: ES (e.g. >72 months HS, etc...)
    eu1HSMiss=round(as.numeric(eu1HSInf)-as.numeric(eu1HS) )  #number missed HS
             eu1HSMiss[eu1HSMiss<0]=0
    eu1ESMiss=round(as.numeric(eu1ESInf)-as.numeric(eu1ES) )  #number missed ES
             eu1ESMiss[eu1ESMiss<0]=0                    
  
    eu1HSInflci=formatC((qgamma(poL,truePrevD$eHSInf,rate=0,1)) ,digits=2)          #lower CI - baseline
    eu1HSInfuci=formatC((qgamma(poU,truePrevD$eHSInf,rate=0,1)) ,digits=2)          #upper CI - basleine
               eu1HSInfuci[as.numeric(eu1HSInfuci)<as.numeric(eu1HSInf)]=NA

    eu1ESInflci=formatC((qgamma(poL,truePrevD$eESInf,rate=0,1)) ,digits=2)          #lower CI - baseline
    eu1ESInfuci=formatC((qgamma(poU,truePrevD$eESInf,rate=0,1)) ,digits=2)          #upper CI - basleine
               eu1ESInfuci[as.numeric(eu1ESInfuci)<as.numeric(eu1ESInf)]=NA
    eu1HSlciMiss=round(as.numeric(eu1HSInflci)-as.numeric(eu1HS)  )
                eu1HSlciMiss[eu1HSlciMiss<0]=0
    eu1HSuciMiss=round(as.numeric(eu1HSInfuci)-as.numeric(eu1HS) )
                eu1HSuciMiss[eu1HSuciMiss<0]=0
     eu1ESlciMiss=round(as.numeric(eu1ESInflci)-as.numeric(eu1ES)  )
                 eu1ESlciMiss[eu1ESlciMiss<0]=0
    eu1ESuciMiss=round(as.numeric(eu1ESInfuci)-as.numeric(eu1ES))
                eu1ESuciMiss[eu1ESuciMiss<0]=0
   #-------------------------------------------------------------------------------------
   ###Outputs for previous years
   
    HistNumInfChtHSPos=reshape(eu1HS,nCht,nTest)       #baseline test positives: HS (e.g. >72 months HS, etc...)
    colnames(HistNumInfChtHSPos)=c(testStart:testEnd)                 
    rownames(HistNumInfChtHSPos)=c(chtStart:chtEnd)
    
    HistNumInfChtESPos=reshape(eu1ES,nCht,nTest)          #baseline test positives: ES (e.g. >72 months HS, etc...)
    colnames(HistNumInfChtESPos)=c(testStart:testEnd)                 
    rownames(HistNumInfChtESPos)=c(chtStart:chtEnd)
    
    HistNumInfChtHS=reshape(paste(eu1HSInf,"[", eu1HSInflci,",",eu1HSInfuci,"]"),nCht,nTest)          #baseline test positives: HS (e.g. >72 months HS, etc...)
    colnames(HistNumInfChtHS)=c(testStart:testEnd)                 
    rownames(HistNumInfChtHS)=c(chtStart:chtEnd)
    
    HistNumInfChtES=reshape(paste(eu1ESInf,"[", eu1ESInflci,",",eu1ESInfuci,"]"),nCht,nTest)          #baseline test positives: ES (e.g. >72 months HS, etc...)
    colnames(HistNumInfChtES)=c(testStart:testEnd)                 
    rownames(HistNumInfChtES)=c(chtStart:chtEnd)
    
    HistNumInfChtHSMiss=reshape(paste(eu1HSMiss,"[", eu1HSlciMiss,",",eu1HSuciMiss,"]"),nCht,nTest)      #number missed HS
    colnames(HistNumInfChtHSMiss)=c(testStart:testEnd)                
    rownames(HistNumInfChtHSMiss)=c(chtStart:chtEnd)
    
    HistNumInfChtESMiss=reshape(paste(eu1ESMiss,"[", eu1ESlciMiss,",",eu1ESuciMiss,"]"),nCht,nTest)     #number missed ES
    colnames(HistNumInfChtESMiss)=c(testStart:testEnd)                 
    rownames(HistNumInfChtESMiss)=c(chtStart:chtEnd)
     
     write.csv(HistNumInfChtHSPos,paste(outDir,"/",MSName,"HistNumInfChtHSPos.csv",sep=""))
     write.csv(HistNumInfChtESPos,paste(outDir,"/",MSName,"HistNumInfChtESPos.csv",sep=""))
     write.csv(HistNumInfChtHS,paste(outDir,"/",MSName,"HistNumInfChtHS.csv",sep=""))
     write.csv(HistNumInfChtES,paste(outDir,"/",MSName,"HistNumInfChtES.csv",sep=""))
     write.csv(HistNumInfChtHSMiss,paste(outDir,"/",MSName,"HistNumInfChtHSMiss.csv",sep=""))
     write.csv(HistNumInfChtESMiss,paste(outDir,"/",MSName,"HistNumInfChtESMiss.csv",sep=""))
     
  