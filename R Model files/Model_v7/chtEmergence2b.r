 ## EFSA BSE: new Emergence
chtEmergence2b<-function(aa,aaT,bb,cc,dd,trendUp0,paras,split0) {

#--------------- Set up matrices -------------------------------------------------------------------
NumInfAgeEmer<-matrix(0,nAgeGrp2,fP+1)    #fP=maximum number of years to predict forward 
rownames(NumInfAgeEmer)=c(rownames(TestStartIn),"p1","p2","p3","p4","p5","p6")
 colnames(NumInfAgeEmer)=c("Current",1:(fP))
NumDetAgeEmer<-NumInfAgeEmer 
NumInfAgeEmerCSFS<-NumInfAgeEmer 
NumDetAgeEmerCSFS<-NumInfAgeEmer 
NumInfAgeEmerHSES<-NumInfAgeEmer
NumDetAgeEmerHSES<-NumInfAgeEmer
NumInfAgeEmerHS<-NumInfAgeEmer 
NumDetAgeEmerHS<-NumInfAgeEmer
NumInfAgeEmerES<-NumInfAgeEmer 
NumDetAgeEmerES<-NumInfAgeEmer
NumInfAgeEmerFS<-NumInfAgeEmer 
NumDetAgeEmerFS<-NumInfAgeEmer 
NumInfAgeEmerCS<-NumInfAgeEmer 
NumDetAgeEmerCS<-NumInfAgeEmer
totDetEmer=  matrix(0,1,fP+1)
totInfEmer=matrix(0,1,fP+1)
# ----------------------  inital calculations ---------------------------------   
## Split up to CS & FS - based on prop test positive
pCSAge<-c(split0$CStp,rep(split0$CStp[length(split0$CStp)],6))              #6 extra values here (p1 - p6): assign them all the last age value of CStp (204 months)  
pFSAge<-c(split0$FStp,rep(split0$FStp[length(split0$FStp)],6))  
pHSAge<-c(split0$HStp,rep(split0$HStp[length(split0$HStp)],6))     
pESAge<-c(split0$EStp,rep(split0$EStp[length(split0$EStp)],6))  

pCSAgeInf<-c(splitTrue$CStp,rep(splitTrue$CStp[length(splitTrue$CStp)],6))  #6 extra values here (p1 - p6): assign them all the last age value of CStp (204 months) 
pFSAgeInf<-c(splitTrue$FStp,rep(splitTrue$FStp[length(splitTrue$FStp)],6)) 
pHSAgeInf<-c(splitTrue$HStp,rep(splitTrue$HStp[length(splitTrue$HStp)],6))    
pESAgeInf<-c(splitTrue$EStp,rep(splitTrue$HStp[length(splitTrue$HStp)],6)) 
        
NumDetAgeEmer[,1]=chtConv1203(aa$eTot[,nTest])                                  # num detected by age in current year
NumInfAgeEmer[,1]=chtConv1203(aaT$eTotInf[,nTest])                              # num infected by age in current year  - based on true prev
NumDetAgeEmerCSFS[,1]<-chtConv1203(aa$eCSFS[,nTest])           ;  NumDetAgeEmerHSES[,1]<-chtConv1203(aa$eHSES[,nTest])
NumInfAgeEmerCSFS[,1]<-chtConv1203(truePrevD$eCSFSInf[,nTest]) ;  NumInfAgeEmerHSES[,1]<-chtConv1203(truePrevD$eHSESInf[,nTest])
      NumDetAgeEmerCS[,1]=NumDetAgeEmerCSFS[,1] * pCSAge;
      NumDetAgeEmerFS[,1]=NumDetAgeEmerCSFS[,1] * pFSAge;
      NumDetAgeEmerHS[,1]=NumDetAgeEmerHSES[,1] * pHSAge;
      NumDetAgeEmerES[,1]=NumDetAgeEmerHSES[,1] * pESAge;
      
      NumInfAgeEmerCS[,1]=NumInfAgeEmerCSFS[,1] * pCSAgeInf;
      NumInfAgeEmerFS[,1]=NumInfAgeEmerCSFS[,1] * pFSAgeInf;
      NumInfAgeEmerHS[,1]=NumInfAgeEmerHSES[,1] * pHSAgeInf;
      NumInfAgeEmerES[,1]=NumInfAgeEmerHSES[,1] * pESAgeInf;
#-------------------- Get trend and K again -----------------------------------------
#paras<-MLEout$paras
jj<-1:(nCht)
ifelse(Wbl==1, {
   A=paras[1]; B=paras[2];C=paras[3]    #weibull trend
   trend=C* ((1:nCht) / A)^(B-1)* exp(-((1:nCht) / A)^B); 
   k=paras[4];K=exp(k)/(1+exp(k));  
},{
  trend<-paras[1]*exp(paras[2]*(jj-1))    #exponential trend
  k=paras[3];K=exp(k)/(1+exp(k));  
})
#extend assuming an x% increase in trend for cohorts born each year after the current year (nCht)
trend[(nCht+1):(nCht+fP)]=trend[nCht]*(1+incRate)^(1:fP)
trend0=trend
#--------------------------------------------------------------------
#Number dead
HSESNd0<-matrix(0,nCht+fP,nTest+fP); CSFSNd0<-matrix(0,nCht+fP,nTest+fP)     
       HSESNd0[1:nCht,1:nTest]= tIn$HSESNumDead ; CSFSNd0[1:nCht,1:nTest]= tIn$CSFSNumDead
#Number Tested
HSESNt0<-matrix(0,nCht+fP,nTest+fP); CSFSNt0<-matrix(0,nCht+fP,nTest+fP)
       HSESNt0[1:nCht,1:nTest]= dd$HSESNumTest ; CSFSNt0[1:nCht,1:nTest]= dd$CSFSNumTest 
#PC/PD/PS          
PC0<-matrix(0,nCht+fP,nTest+fP); PS0<-matrix(0,nCht+fP,nTest+fP) ;  PD0<-matrix(0,nCht+fP,nTest+fP)
       PC0[1:nCht,1:nTest]= PC ; PS0[1:nCht,1:nTest]= PS ; PD0[1:nCht,1:nTest]= PD
#Number detected - expected, use trend       
NumDetEmerHSES<-matrix(0,nCht+fP,nTest+fP); NumDetEmerCSFS<-matrix(0,nCht+fP,nTest+fP)
       NumDetEmerHSES[1:nCht,1:nTest]<-aa$eHSES ; NumDetEmerCSFS[1:nCht,1:nTest]<-aa$eCSFS
#Number Infected - expected, use trend      
NumInfEmerHSES<-matrix(0,nCht+fP,nTest+fP); NumInfEmerCSFS<-matrix(0,nCht+fP,nTest+fP)
       NumInfEmerHSES[1:nCht,1:nTest]<-aa$eHSESInf ; NumInfEmerCSFS[1:nCht,1:nTest]<-aa$eCSFSInf
#----------------- Extend Number of animals/PC/PS/PD matrices forward fP years ------------------------------     
for(i in 1:fP) {
    HSESNd0[,nTest+i]=c(0,HSESNd0[1:(nCht+fP-1),nTest+i-1] )
    CSFSNd0[,nTest+i]=c(0,CSFSNd0[1:(nCht+fP-1),nTest+i-1] )
    HSESNt0[,nTest+i]=c(0,HSESNt0[1:(nCht+fP-1),nTest+i-1] )
    CSFSNt0[,nTest+i]=c(0,CSFSNt0[1:(nCht+fP-1),nTest+i-1] )    
    PS0[1:(nCht+i),nTest+i]=PsurvAge[(ageMat[1,nTest]+i+1):1]        # new PS values for  year nTest+i - animals > 42 years will get NA value: we dont care about these, assume there are no animals this old.
    PC0[1:(nCht+i),nTest+i]= PclinAge[(ageMat[1,nTest]+i+1):1]       # new PC values for  year nTest+i - animals > 42 years will get NA value: we dont care about these, assume there are no animals this old.
    PD0[1:(nCht+i),nTest+i]=  PdetectedAge[(ageMat[1,nTest]+i+1):1]  # new PD values for  year nTest+i - animals > 42 years will get NA value: we dont care about these, assume there are no animals this old.
       }    #end
# data formatting      
PC0[is.na(PC0)]=0  ; PS0[is.na(PS0)]=0 ; PD0[is.na(PD0)]=0  #NA's = 0 probability    
#----------------- Extend Number detected/infected forward fP years ------------------------------    
for(i in 1:fP) { 
      trend1=trend[1:(nCht+i)] #get trend up to future year i 
#update test risk rates
        rr1<-   trend1*(1-K)*PC0[1:(nCht+i),nTest+i]*MaxSensitivity
        rrH1<-PD0[1:(nCht+i),nTest+i]*(1-PS0[1:(nCht+i),nTest+i])*trend1*K
#update infected risk rates 
        rrInf1<-trend1*(1-K)*PC0[1:(nCht+i),nTest+i]                   #new CSFS risk rate for year nTest +i
        rrHInf1<-trend1*K*(1-PS0[1:(nCht+i),nTest+i])                 #new HSES risk rate for year nTest +i  
#update number detected
        NumDetEmerHSES[1:(nCht+i),(nTest+i)] = rrH1*HSESNt0[1:(nCht+i),(nTest+i)]
        NumDetEmerCSFS[1:(nCht+i),(nTest+i)] = rr1*CSFSNt0[1:(nCht+i),(nTest+i)]      
#update number infected 
        NumInfEmerHSES[1:(nCht+i),(nTest+i)] = rrHInf1*HSESNd0[1:(nCht+i),(nTest+i)]
        NumInfEmerCSFS[1:(nCht+i),(nTest+i)] = rrInf1*CSFSNd0[1:(nCht+i),(nTest+i)]
}     #end                                       
totInfEmer[1]=sum( (NumInfAgeEmer[,1]))  
totDetEmer[1]=sum( (NumDetAgeEmer[,1]))
NumDetTot<-NumDetEmerHSES+NumDetEmerCSFS     #total number detected
nSumTot<-colSums(NumDetTot)                  #sum total per testing period
#-------------- Convert to by Age --------------------------------------------------
for (i in 2:(fP+1))  {
 if ( (i/50)==(i %/% 50)) print (i)
  NumDetAgeEmerCSFS[,i]<-chtConv1203(NumDetEmerCSFS[1:(nCht+i-1),nTest+i-1]) ;  NumDetAgeEmerHSES[,i]<-chtConv1203(NumDetEmerHSES[1:(nCht+i-1),nTest+i-1])
  NumInfAgeEmerCSFS[,i]<-chtConv1203(NumInfEmerCSFS[1:(nCht+i-1),nTest+i-1]) ;  NumInfAgeEmerHSES[,i]<-chtConv1203(NumInfEmerHSES[1:(nCht+i-1),nTest+i-1])
        NumDetAgeEmerCS[,i]=NumDetAgeEmerCSFS[,i] * pCSAge;
        NumDetAgeEmerFS[,i]=NumDetAgeEmerCSFS[,i] * pFSAge;
        NumDetAgeEmerHS[,i]=NumDetAgeEmerHSES[,i] * pHSAge;
        NumDetAgeEmerES[,i]=NumDetAgeEmerHSES[,i] * pESAge;
        
        NumInfAgeEmerCS[,i]=NumInfAgeEmerCSFS[,i] * pCSAgeInf;
        NumInfAgeEmerFS[,i]=NumInfAgeEmerCSFS[,i] * pFSAgeInf;
        NumInfAgeEmerHS[,i]=NumInfAgeEmerHSES[,i] * pHSAgeInf;
        NumInfAgeEmerES[,i]=NumInfAgeEmerHSES[,i] * pESAgeInf;
  
  NumDetAgeEmer[,i]<-NumDetAgeEmerCSFS[,i]+NumDetAgeEmerHSES[,i]
  NumInfAgeEmer[,i]<-NumInfAgeEmerCSFS[,i]+NumInfAgeEmerHSES[,i]
  totInfEmer[i]=sum(NumInfAgeEmer[,i],na.rm=TRUE)  #will get NAs after ~16 rows, remove them
 totDetEmer[i]=sum(NumDetAgeEmer[,i],na.rm=TRUE)
  }
#*******************************************************************************
trendUp=eval(trendUp0)
nSkip<-max(1,find(diff(nSumTot[nTest:(nTest+fP)])>=0)[1] )  #find when expected total stops decreasing - minus off up to current year  - dont let emergence detection happen until after this point
if(is.na(nSkip)) nSkip<-fP
AgeDetEmer=AgeDetEmer2   

#----------------- Expected number detected - deterministic: dont use------------------------------  
#this method is just using the baseline values - it will give a much higher output than the simulations and is not as good.  
#DONT COMPARE OUTPUT OF THIS WITH CIs OF SIMULATION METHOD
YrDetD<-nSkip-1
emerCV<-round(trendUp)-1
while (emerCV<round(trendUp) & YrDetD<fP) {
         emerCV<- sum(round(NumDetAgeEmer[AgeDetEmer,YrDetD+1]))      #+1 because first column is current year
        YrDetD<-YrDetD+1
} # end while emerCV<=round(trendUp)
#----------------- Simulate number detected ------------------------------  
nSim<-10000                                       #number of simulated runs to do for number detected
YrDetSim<-matrix(nSkip-1,1,nSim)                    #number of years since increasing trend started 
for (i in 1:nSim)   {
    emerCV<-round(trendUp)-1                       #set inital value of emerCV to something that wont exit while loop
    simDetAgeEmer<-reshape(rpois(matrix(0,nAgeGrp2,fP+1),NumDetAgeEmer),nAgeGrp2,fP+1)
  while (emerCV<round(trendUp) & YrDetSim[i]<fP) {
         emerCV<- sum(round(simDetAgeEmer[AgeDetEmer,YrDetSim[i]+1]))      #+1 because first column is current year
         YrDetSim[i]<-YrDetSim[i]+1
  } # end while emerCV<=round(trendUp)
} #end for (i in 1:nSim) 
YrDetMed<- ceil(quantile(YrDetSim,0.5))    #Median (50th pctile) of simulated detection years -round up to nearest year
YrDetAvg<- ceil(mean(YrDetSim))            #Mean number of simulated detection years -round up to nearest year

YrDetLci<- ceil(quantile(YrDetSim,0.025)) #lower confidence interval = 2.5th percentile of simulated detection years -round up to nearest year
YrDetUci<- ceil(quantile(YrDetSim,0.975)) #upper confidence interval = 97.5th percentile of simulated detection years - round up to nearest year
#*****************************************************************************   
#--------------------- Output   ------------------------------------
i<-YrDetAvg

NumMissedHSES=(NumInfAgeEmerHSES-NumDetAgeEmerHSES )

totByYr<-matrix(0,3,i+1)
colnames(totByYr)=c("Total","Current",1:(i-1))
rownames(totByYr)=c("Detected","Infected","Missed")
totByYr[1,2:(i+1)]=totDetEmer[1:i] ; totByYr[1,1]=sum(totByYr[1,3:(i+1)])   #sum from 3rd place as dont want the current year
totByYr[2,2:(i+1)]=totInfEmer[1:i] ; totByYr[2,1]=sum(totByYr[2,3:(i+1)])
totByYr[3,2:(i+1)]=colSums(NumMissedHSES[,1:i]) ; totByYr[3,1]=sum(totByYr[3,3:(i+1)])
# Number Missed Output
  EmerAABase<-matrix(0,nAgeGrp2,6)
  rownames(EmerAABase)=c(rownames(TestStartIn),"p1","p2","p3","p4","p5","p6")
  colnames(EmerAABase)=c("nMissDetYr","nMissDetYrLci","nMissDetYrUci","nMissSum","nMissSumLci","nMissSumUci")
    EmerAABase[,1]<-NumMissedHSES[,YrDetAvg+1]
    EmerAABase[,2]<-NumMissedHSES[,YrDetLci+1]
    EmerAABase[,3]<-NumMissedHSES[,YrDetUci+1]
    EmerAABase[,4]<-rowSums(NumMissedHSES[,2:YrDetAvg])           #row 1=current
    EmerAABase[,5]<-rowSums(NumMissedHSES[,2:YrDetLci] )
    EmerAABase[,6]<-rowSums(NumMissedHSES[,2:YrDetUci])
    ## ISSUE:  The lower confidence interval of the year could be at a time when there is an increasing trend in younger animals but still a decreasing trend in older animals
# In these cases the upperCI of number missed for those older animals will actually be higher than the lower CI 
# This is not useful for the model so we replace with NaNs

EmerAABase[(EmerAABase[,3]<EmerAABase[,2] |  EmerAABase[,3]<EmerAABase[,1]),3]=NA
EmerAABase[(EmerAABase[,2]>EmerAABase[,3] |  EmerAABase[,2]>EmerAABase[,1]),2]=NA 
EmerAABase[(EmerAABase[,6]<EmerAABase[,5] |  EmerAABase[,6]<EmerAABase[,4]),6]=NA
EmerAABase[(EmerAABase[,5]>EmerAABase[,6] |  EmerAABase[,5]>EmerAABase[,4]),5]=NA 

  write.csv(EmerAABase,paste(outDir,"/Emer2013/",MS,"nMissEmer",cc,".csv",sep=""))
# write to .csv files 
if (!file.exists(file.path(outDir,"/Emer2013/",cc))) dir.create( file.path(outDir,"/Emer2013/",cc))      
write.csv((NumInfAgeEmer),paste(outDir,"/Emer2013/",cc,"/",MS,"nInfEmer(",cc,trendMeth,").csv",sep=""))
write.csv((NumDetAgeEmer),paste(outDir,"/Emer2013/",cc,"/",MS,"nDetEmer(",cc,trendMeth,").csv",sep=""))
write.csv((NumMissedHSES),paste(outDir,"/Emer2013/",cc,"/",MS,"nMissedHSES(",cc,trendMeth,").csv",sep=""))
write.csv((totByYr),paste(outDir,"/Emer2013/",cc,"/",MS,"totByYr(",cc,trendMeth,").csv",sep=""))

write.csv((NumInfAgeEmerHS),paste(outDir,"/Emer2013/",cc,"/",MS,"nInfEmerHS(",cc,trendMeth,").csv",sep=""))
write.csv((NumInfAgeEmerES),paste(outDir,"/Emer2013/",cc,"/",MS,"nInfEmerES(",cc,trendMeth,").csv",sep=""))
write.csv((NumInfAgeEmerFS),paste(outDir,"/Emer2013/",cc,"/",MS,"nInfEmerFS(",cc,trendMeth,").csv",sep=""))
write.csv((NumInfAgeEmerCS),paste(outDir,"/Emer2013/",cc,"/",MS,"nInfEmerCS(",cc,trendMeth,").csv",sep=""))
write.csv((NumInfAgeEmerCSFS),paste(outDir,"/Emer2013/",cc,"/",MS,"nInfEmerCSFS(",cc,trendMeth,").csv",sep=""))
write.csv((NumInfAgeEmerHSES),paste(outDir,"/Emer2013/",cc,"/",MS,"nInfEmerHSES(",cc,trendMeth,").csv",sep=""))

write.csv((NumDetAgeEmerHS),paste(outDir,"/Emer2013/",cc,"/",MS,"nDetEmerHS(",cc,trendMeth,").csv",sep=""))
write.csv((NumDetAgeEmerES),paste(outDir,"/Emer2013/",cc,"/",MS,"nDetEmerES(",cc,trendMeth,").csv",sep=""))
write.csv((NumDetAgeEmerFS),paste(outDir,"/Emer2013/",cc,"/",MS,"nDetEmerFS(",cc,trendMeth,").csv",sep=""))
write.csv((NumDetAgeEmerCS),paste(outDir,"/Emer2013/",cc,"/",MS,"nDetEmerCS(",cc,trendMeth,").csv",sep=""))
write.csv((NumDetAgeEmerHSES),paste(outDir,"/Emer2013/",cc,"/",MS,"nDetEmerHSES(",cc,trendMeth,").csv",sep=""))
write.csv((NumDetAgeEmerCSFS),paste(outDir,"/Emer2013/",cc,"/",MS,"nDetEmerCSFS(",cc,trendMeth,").csv",sep=""))

write.csv((NumInfAgeEmerHS-NumDetAgeEmerHS),paste(outDir,"/Emer2013/",cc,"/",MS,"nMissedHS(",cc,trendMeth,").csv",sep=""))
write.csv((NumInfAgeEmerES-NumDetAgeEmerES),paste(outDir,"/Emer2013/",cc,"/",MS,"nMissedES(",cc,trendMeth,").csv",sep=""))
write.csv((NumInfAgeEmerFS-NumDetAgeEmerFS),paste(outDir,"/Emer2013/",cc,"/",MS,"nMissedFS(",cc,trendMeth,").csv",sep=""))
write.csv((NumInfAgeEmerCS-NumDetAgeEmerCS),paste(outDir,"/Emer2013/",cc,"/",MS,"nMissedCS(",cc,trendMeth,").csv",sep=""))

list(nInf=NumInfAgeEmer,nDet=NumDetAgeEmer,nMiss=NumMissedHSES,tot=totByYr,yrDet=YrDetAvg,yrDetLci=YrDetLci,yrDetUci=YrDetUci)  
}  # end function
