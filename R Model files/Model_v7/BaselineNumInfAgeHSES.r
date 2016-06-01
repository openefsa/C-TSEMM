    # Most recent year number infected by age for HS and ES exit streams
    
    BaselineNumInfAge=matrix(0,nAgeGrp2,18)
    colnames(BaselineNumInfAge)=c("HS detected","HS detected 2.5th", "HS detected 97.5th","ES detected","ES detected 2.5th", "ES detected 97.5th","HS infected","HS infected 2.5th", "HS infected 97.5th","ES infected","ES infected 2.5th", "ES infected 97.5th","HS missed","HS missed 2.5th", "HS missed 97.5th","ES missed","ES missed 2.5th", "ES missed 97.5th")
    rownames(BaselineNumInfAge)=c(rownames(TestStartIn),"p1","p2","p3","p4","p5","p6")

    eu1HS=formatC(chtConv1203(baselineD$eHS[,nTest]),digits=2)

    eu1ES=formatC(chtConv1203(baselineD$eES[,nTest]),digits=2)        #baseline test positives: ES (e.g. >72 months HS, etc...)

    eu1HSInf=formatC(chtConv1203(truePrevD$eHSInf[,nTest]),digits=2)        #baseline infecteds positives: HS (e.g. >72 months HS, etc...)
    eu1ESInf=formatC(chtConv1203(truePrevD$eESInf[,nTest]),digits=2)        #baseline infecteds positives: ES (e.g. >72 months HS, etc...)
    eu1HSMiss=round(as.numeric(eu1HSInf)-as.numeric(eu1HS))               #number missed HS
           eu1HSMiss[eu1HSMiss<0]=0
    eu1ESMiss=round(as.numeric(eu1ESInf)-as.numeric(eu1ES))              #number missed ES
           eu1ESMiss[ eu1ESMiss<0]=0
   # Confidence interval outputs
    eu1HSlci=formatC((qgamma(poL,as.numeric(eu1HS),rate=0,1)) ,digits=2)        #lower CI - baseline
    eu1HSuci=formatC((qgamma(poU,as.numeric(eu1HS),rate=0,1)),digits=2)         #upper CI - basleine
        eu1HSuci[as.numeric(eu1HSuci)<as.numeric(eu1HS)]=NA
    eu1ESlci=formatC((qgamma(poL,as.numeric(eu1ES),rate=0,1)),digits=2)        #lower CI - baseline
    eu1ESuci=formatC((qgamma(poU,as.numeric(eu1ES),rate=0,1)),digits=2)         #upper CI - basleine
        eu1ESuci[as.numeric(eu1ESuci)<as.numeric(eu1ES)]=NA

    eu1HSInflci=formatC(chtConv1203(qgamma(poL,truePrevD$eHSInf,rate=0,1)[,nTest]) ,digits=2)        #lower CI - baseline
    eu1HSInfuci=formatC(chtConv1203(qgamma(poU,truePrevD$eHSInf,rate=0,1)[,nTest]) ,digits=2)       #upper CI - basleine
                       eu1HSInfuci[as.numeric(eu1HSInfuci)<as.numeric(eu1HSInf)]=NA

    eu1ESInflci=formatC(chtConv1203(qgamma(poL,truePrevD$eESInf,rate=0,1) [,nTest]) ,digits=2)       #lower CI - baseline
    eu1ESInfuci=formatC(chtConv1203(qgamma(poU,truePrevD$eESInf,rate=0,1) [,nTest]) ,digits=2)      #upper CI - basleine
                             eu1ESInfuci[as.numeric(eu1ESInfuci)<as.numeric(eu1ESInf)]=NA
    eu1HSlciMiss=round(as.numeric(eu1HSInflci)-as.numeric(eu1HSlci))
    eu1HSlciMiss[eu1HSlciMiss<0]=0
    eu1HSuciMiss=round(as.numeric(eu1HSInfuci)-as.numeric(eu1HSuci) )
           eu1HSuciMiss[eu1HSuciMiss<0]=0
     eu1ESlciMiss=round(as.numeric(eu1ESInflci)-as.numeric(eu1ESlci))
          eu1ESlciMiss[eu1ESlciMiss<0]=0
    eu1ESuciMiss=round(as.numeric(eu1ESInfuci)-as.numeric(eu1ESuci) )
             eu1ESuciMiss[ eu1ESuciMiss<0]=0
             
    BaselineNumInfAge[,1]=paste(eu1HS)       #baseline test positives: HS (e.g. >72 months HS, etc...)
       BaselineNumInfAge[,2]=paste(eu1HSlci) 
       BaselineNumInfAge[,3]=paste(eu1HSuci) 
    BaselineNumInfAge[,4]=paste(eu1ES)        #baseline test positives: ES (e.g. >72 months HS, etc...)
      BaselineNumInfAge[,5]=paste(eu1ESlci)
      BaselineNumInfAge[,6]=paste(eu1ESuci)
    BaselineNumInfAge[,7]=paste(eu1HSInf)        #baseline test positives: HS (e.g. >72 months HS, etc...)
      BaselineNumInfAge[,8]=paste(eu1HSInflci)
      BaselineNumInfAge[,9]=paste(eu1HSInfuci)
    BaselineNumInfAge[,10]=paste(eu1ESInf)        #baseline test positives: ES (e.g. >72 months HS, etc...)
      BaselineNumInfAge[,11]=paste(eu1ESInflci)
      BaselineNumInfAge[,12]=paste(eu1ESInfuci)
    BaselineNumInfAge[,13]=paste(eu1HSMiss)    #number missed HS
      BaselineNumInfAge[,14]=paste(eu1HSlciMiss)
      BaselineNumInfAge[,15]=paste(eu1HSuciMiss)
    BaselineNumInfAge[,16]=paste(eu1ESMiss)   #number missed ES
        BaselineNumInfAge[,17]=paste(eu1ESlciMiss)
        BaselineNumInfAge[,18]=paste(eu1ESuciMiss)
    
    write.csv(BaselineNumInfAge,paste(outDir,"/",MSName,"BaselineNumInfAge.csv",sep=""))

