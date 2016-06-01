###############################################################################
#### All years number infected by age for HS and ES exit streams   ############
###############################################################################       
### ----------------------- SET up output matrices ----------------------------              

eu1HSInf=matrix(0,nAgeGrp2,nTest)
eu1ESInf=matrix(0,nAgeGrp2,nTest)
eu1HSInflci=matrix(0,nAgeGrp2,nTest)
eu1HSInfuci=matrix(0,nAgeGrp2,nTest)
eu1ESInflci=matrix(0,nAgeGrp2,nTest)
eu1ESInfuci=matrix(0,nAgeGrp2,nTest)

#---Convert Cohort to Age based output -----------------------------------------  
for (i in 1:nTest) {                            # need for loop as will only work on one column at a time
        
    eu1HSInf[,i]=chtConv1203(truePrevD$eHSInf[1:(nCht-nTest+i),i])        #baseline test positives: HS (e.g. >72 months HS, etc...)
    eu1ESInf[,i]=chtConv1203(truePrevD$eESInf[1:(nCht-nTest+i),i])        #baseline test positives: ES (e.g. >72 months HS, etc...)
    
    eu1HSInflci[,i]=chtConv1203(qgamma(poL,truePrevD$eHSInf[1:(nCht-nTest+i),i],rate=0,1))       #lower CI - baseline
    eu1HSInfuci[,i]=chtConv1203(qgamma(poU,truePrevD$eHSInf[1:(nCht-nTest+i),i],rate=0,1))        #upper CI - basleine
    eu1ESInflci[,i]=chtConv1203(qgamma(poL,truePrevD$eESInf[1:(nCht-nTest+i),i],rate=0,1))       #lower CI - baseline
    eu1ESInfuci[,i]=chtConv1203(qgamma(poU,truePrevD$eESInf[1:(nCht-nTest+i),i],rate=0,1))       #upper CI - basleine
 }
    
eu1HSInfuci[as.numeric(eu1HSInfuci)<as.numeric(eu1HSInf)]=NA          # set weird outputs to NA (-ve values etc...)
eu1ESInfuci[as.numeric(eu1ESInfuci)<as.numeric(eu1ESInf)]=NA         # set weird ouputs to NA  (-ve values etc...)
  
eu1HSMiss=round(as.numeric(eu1HSInf[1:nAgeGrp,])-(HSbyAgePos))   #number missed HS
          eu1HSMiss[eu1HSMiss<0]=0
eu1ESMiss=round(as.numeric(eu1ESInf[1:nAgeGrp,])-(ESbyAgePos))   #number missed ES   
          eu1ESMiss[eu1ESMiss<0]=0    
eu1HSlciMiss=round(as.numeric(eu1HSInflci[1:nAgeGrp,])-(HSbyAgePos))
          eu1HSlciMiss[eu1HSlciMiss<0]=0
eu1HSuciMiss=round(as.numeric(eu1HSInfuci[1:nAgeGrp,])-(HSbyAgePos))
          eu1HSuciMiss[eu1HSuciMiss<0]=0
eu1ESlciMiss=round(as.numeric(eu1ESInflci[1:nAgeGrp,])-(ESbyAgePos))
          eu1ESlciMiss[eu1ESlciMiss<0]=0
eu1ESuciMiss=round(as.numeric(eu1ESInfuci[1:nAgeGrp,])-(ESbyAgePos))
          eu1ESuciMiss[eu1ESuciMiss<0]=0

# ------------------- Outputs for previous years ---------------------------------------------------
    HistNumInfAgeHS=reshape(paste(formatC(eu1HSInf,digits=2),"[", formatC(eu1HSInflci,digits=2),",",formatC(eu1HSInfuci,digits=2),"]"),nAgeGrp2,nTest)   
    colnames(HistNumInfAgeHS)=c(testStart:testEnd)            
    rownames(HistNumInfAgeHS)=c(rownames(TestStartIn),"p1","p2","p3","p4","p5","p6")
    
    HistNumInfAgeES=reshape(paste(formatC(eu1ESInf,digits=2),"[", formatC(eu1ESInflci,digits=2),",",formatC(eu1ESInfuci,digits=2),"]"),nAgeGrp2,nTest)          #baseline test positives: ES (e.g. >72 months HS, etc...)
    colnames(HistNumInfAgeES)=c(testStart:testEnd)               
    rownames(HistNumInfAgeES)=c(rownames(TestStartIn),"p1","p2","p3","p4","p5","p6")
   
    #put value, lci & uci in same matrix 
        HistNumInfAgeHS<-cbind(eu1HSInf[,1],eu1HSInflci[,1],eu1HSInfuci[,1])
       for (i in 2:nTest)  {HistNumInfAgeHS<-cbind(HistNumInfAgeHS,eu1HSInf[,i],eu1HSInflci[,i],eu1HSInfuci[,i])  }
       colnames(HistNumInfAgeHS)=paste(c(rep(testStart:testEnd,1,each=3)),  c(rep(c("Base","lci","uci"),nTest)))             
    rownames(HistNumInfAgeHS)=c(rownames(TestStartIn),"p1","p2","p3","p4","p5","p6")
    
      HistNumInfAgeES<-cbind(eu1ESInf[,1],eu1ESInflci[,1],eu1ESInfuci[,1])
       for (i in 2:nTest)  {HistNumInfAgeES<-cbind(HistNumInfAgeES,eu1ESInf[,i],eu1ESInflci[,i],eu1ESInfuci[,i])  }
       colnames(HistNumInfAgeES)=paste(c(rep(testStart:testEnd,1,each=3)),  c(rep(c("Base","lci","uci"),nTest)))             
    rownames(HistNumInfAgeES)=c(rownames(TestStartIn),"p1","p2","p3","p4","p5","p6")
   ## miss   
           HistNumInfAgeHSMiss<-cbind(eu1HSMiss[,1],eu1HSlciMiss[,1],eu1HSuciMiss[,1])
       for (i in 2:nTest)  {HistNumInfAgeHSMiss<-cbind(HistNumInfAgeHSMiss,eu1HSMiss[,i],eu1HSlciMiss[,i],eu1HSuciMiss[,i])  }
       colnames(HistNumInfAgeHSMiss)=paste(c(rep(testStart:testEnd,1,each=3)),  c(rep(c("BaseMiss","lci","uci"),nTest)))             
    rownames(HistNumInfAgeHSMiss)=c(rownames(TestStartIn))
    
      HistNumInfAgeESMiss<-cbind(eu1ESMiss[,1],eu1ESlciMiss[,1],eu1ESuciMiss[,1])
       for (i in 2:nTest)  {HistNumInfAgeESMiss<-cbind(HistNumInfAgeESMiss,eu1ESMiss[,i],eu1ESlciMiss[,i],eu1ESuciMiss[,i])  }
       colnames(HistNumInfAgeESMiss)=paste(c(rep(testStart:testEnd,1,each=3)),  c(rep(c("Base","lci","uci"),nTest)))             
    rownames(HistNumInfAgeESMiss)=c(rownames(TestStartIn))
     
     write.csv(HSbyAgePos,paste(outDir,"/",MSName,"HSbyAgePos.csv",sep=""))
     write.csv(ESbyAgePos,paste(outDir,"/",MSName,"ESbyAgePos.csv",sep=""))
     write.csv(HistNumInfAgeHS,paste(outDir,"/",MSName,"HistNumInfAgeHS.csv",sep=""))
     write.csv(HistNumInfAgeES,paste(outDir,"/",MSName,"HistNumInfAgeES.csv",sep=""))
     write.csv(HistNumInfAgeHSMiss,paste(outDir,"/",MSName,"HistNumInfAgeHSMiss.csv",sep=""))
     write.csv(HistNumInfAgeESMiss,paste(outDir,"/",MSName,"HistNumInfAgeESMiss.csv",sep=""))
                                                                