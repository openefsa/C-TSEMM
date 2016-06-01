################# Input Matrices ###############################################
inputData<-function(ageCap,ageCapEnd,prTest)   {
  CSbyAge1<-ScenarioIn(CSbyAge,ageCap[1],ageCapEnd[1],prTest[1]);  ESbyAge1<-ScenarioIn(ESbyAge,ageCap[2],ageCapEnd[2],prTest[2])
  FSbyAge1<-ScenarioIn(FSbyAge,ageCap[3],ageCapEnd[3],prTest[3]);  HSbyAge1<-ScenarioIn(HSbyAge,ageCap[4],ageCapEnd[4],prTest[4])

  ESNumDeadbyAge1<-ScenarioIn(ESNumDeadbyAge,ageCap[2],ageCapEnd[2],prTest[2])       #note: CSNumDead=CSNumTest by definition
  FSNumDeadbyAge1<-ScenarioIn(FSNumDeadbyAge,ageCap[3],ageCapEnd[3],prTest[3]) ;  HSNumDeadbyAge1<-ScenarioIn(HSNumDeadbyAge,ageCap[4],ageCapEnd[4],prTest[4]) 

  CSbyAgePos1<-ScenarioIn(CSbyAgePos,ageCap[1],ageCapEnd[1],prTest[1]) ;  ESbyAgePos1<-ScenarioIn(ESbyAgePos,ageCap[2],ageCapEnd[2],prTest[2])
  FSbyAgePos1<-ScenarioIn(FSbyAgePos,ageCap[3],ageCapEnd[3],prTest[3]) ;  HSbyAgePos1<-ScenarioIn(HSbyAgePos,ageCap[4],ageCapEnd[4],prTest[4])
 
 ######### Data converted to by cohort rather than by age
  CSnt=AgeConv(CSbyAge1);  ESnt=AgeConv(ESbyAge1) ;  FSnt=AgeConv(FSbyAge1) ;  HSnt=AgeConv(HSbyAge1)
  ESnd=AgeConv( ESNumDeadbyAge1) ;  FSnd=AgeConv( FSNumDeadbyAge1) ;  HSnd=AgeConv( HSNumDeadbyAge1)
  CStp=AgeConv(CSbyAgePos1);  EStp=AgeConv(ESbyAgePos1) ;  FStp=AgeConv(FSbyAgePos1)    ;  HStp=AgeConv(HSbyAgePos1)
 
  CSFStp=CStp+FStp; HSEStp=HStp+EStp  ; CSFSnt=CSnt+FSnt; HSESnt=HSnt+ESnt ; CSFSnd=CSnt+FSnd; HSESnd=HSnd+ESnd
   list(CSNumTest=CSnt,ESNumTest=ESnt,FSNumTest=FSnt,HSNumTest=HSnt,ESNumDead=ESnd,FSNumDead=FSnd,HSNumDead=HSnd,CSFSNumTest=CSFSnt,HSESNumTest=HSESnt,CSFSNumDead=CSFSnd,HSESNumDead=HSESnd,CSTestPos=CStp,ESTestPos=EStp,FSTestPos=FStp,HSTestPos=HStp,CSFSTestPos=CSFStp,HSESTestPos=HSEStp)
  }
  ########### Poisson Confidence intervals ######################################
poisInUCI<-function(bb,cc)   {
     ### Calculate upper and lower cofidence intervals of observed data      
      CSUci<- qgamma(cc,bb$CSTestPos+1,1);   CSUci[CSTestPos==0]=0  
      FSUci<- qgamma(cc,bb$CSTestPos+1,1);    FSUci[FSTestPos==0]=0  
      HSUci<- qgamma(cc,bb$HSTestPos+1,1);   HSUci[HSTestPos==0]=0  
      ESUci<- qgamma(cc,bb$CSTestPos+1,1);    ESUci[ESTestPos==0]=0  
      CSFSUci<- qgamma(cc,bb$CSFSTestPos+1,1);  CSFSUci[CSFSTestPos==0]=0  
      HSESUci<- qgamma(cc,bb$HSESTestPos+1,1); HSESUci[HSESTestPos==0]=0  
      ESCSFSUci<- qgamma(cc,bb$ESCSFSTestPos+1,1); ESCSFSUci[ESCSFSTestPos==0]=0 
      list(CSNumTest=bb$CSNumTest,ESNumTest=bb$ESNumTest,FSNumTest=bb$FSNumTest,HSNumTest=bb$HSNumTest,ESNumDead=bb$ESNumDead,FSNumDead=bb$FSNumDead,HSNumDead=bb$HSNumDead,CSFSNumTest=bb$CSFSNumTest,HSESNumTest=bb$HSESNumTest,CSFSNumDead=bb$CSFSNumDead,HSESNumDead=bb$HSESNumDead,CSTestPos=CSUci,ESTestPos=ESUci,FSTestPos=FSUci,HSTestPos=HSUci,CSFSTestPos=CSFSUci,HSESTestPos=HSESUci)
  }
poisInLCI<-function(bb,cc)   {
      CSLci<- qgamma(cc,bb$CSTestPos+1,1);    CSLci[CSTestPos==0]=0   #otherwise lowerCI much greater than 0 when there are 0 test positives
      FSLci<- qgamma(cc,bb$FSTestPos+1,1);    FSLci[FSTestPos==0]=0   #
      HSLci<- qgamma(cc,bb$HSTestPos+1,1);     HSLci[HSTestPos==0]=0  #otherwise lowerCI>0
      ESLci<- qgamma(cc,bb$ESTestPos+1,1);    ESLci[ESTestPos==0]=0   #
      CSFSLci<- qgamma(cc,bb$CSFSTestPos+1,1);    CSFSLci[CSFSTestPos==0]=0   #otherwise lowerCI>0   
      HSESLci<- qgamma(cc,bb$HSESTestPos+1,1);     HSESLci[HSESTestPos==0]=0  #otherwise lowerCI>0
      ESCSFSLci<- qgamma(cc,bb$ESCSFSTestPos+1,1);    ESCSFSLci[ESCSFSTestPos==0]=0   #otherwise lowerCI>0  
      list(CSNumTest=bb$CSNumTest,ESNumTest=bb$ESNumTest,FSNumTest=bb$FSNumTest,HSNumTest=bb$HSNumTest,ESNumDead=bb$ESNumDead,FSNumDead=bb$FSNumDead,HSNumDead=bb$HSNumDead,CSFSNumTest=bb$CSFSNumTest,HSESNumTest=bb$HSESNumTest,CSFSNumDead=bb$CSFSNumDead,HSESNumDead=bb$HSESNumDead,CSTestPos=CSLci,ESTestPos=ESLci,FSTestPos=FSLci,HSTestPos=HSLci,CSFSTestPos=CSFSLci,HSESTestPos=HSESLci)
  }
### Risk Rate calculations    ######################################################
RiskRate<-function(paras,Wbl) {
## calculate trends   -------------------------------------------------------- 
    jj<-1:(nCht)
   ifelse(Wbl==1, {A=paras[1]; B=paras[2];C=paras[3]
      trend=C* ((1:nCht) / A)^(B-1)* exp(-((1:nCht) / A)^B); 
       k=paras[4];K=exp(k)/(1+exp(k));  
       },{
    trend<-paras[1]*exp(paras[2]*(jj-1))
    k=paras[3];K=exp(k)/(1+exp(k));  
    })
#-----------Test positive risk rate ---------------------------------------------------                                    
    #ESTIMATE CSFS risk rate: 
    RiskRate=PC *trend*(1-K) *MaxSensitivity
        RiskRate[RiskRate<0]=0 ;     #Cant realistically have risk<0 ; values likely due to rounding error 
    #ESTIMATE HSES risk rate: 
    RiskRateHS <- PD*(1-PS)*trend*K
        RiskRateHS[RiskRateHS<0]=0;   #Cant realistically have risk<0 ; values likely due to rounding error
#-------------TRUE/Infected RISK RATE (i.e. not accounting for probability of detection) --------------------
      RiskRateTrue=PC *trend*(1-K) 
          RiskRateTrue[RiskRateTrue<0]=0 ;     #Cant realistically have risk<0 ; values likely due to rounding error
      RiskRateTrueHS= PC*(1-PS)*trend*K
          RiskRateTrueHS[RiskRateTrueHS<0]=0;   #Cant realistically have risk<0 ; values likely due to rounding error 
#---------------Stnading Population prevalence ----------------------------------------------------------------
    numTrend<-numel(trend)
    numStand<-  dim(standPop)[1]
    numStandEx<-   numTrend-dim(standPop)[1]
    standPop2<-matrix(0,numTrend,1)
    standPop2[1:numStand]=standPop[1:numStand,]
   
standPrev=PD*(1-PS)*trend *standPop2[numTrend:1,]
#-------------- Infection prevalence Matrix ----------------------------
#matrix of infection prevalence by birth cohort: PREVALENCE IS THE SAME REGARDLESS OF AGE OF ANIMAL

infPrevMat=trend*(1-K)*PC
  
 usePS=1;    
ifelse(usePS==1, { infPrevH=trend*(1-PS)*(K); infPrevH[infPrevH>1]=1  ;  infPrevHMat=infPrevH},
                 {  infPrevH=trend*(K); infPrevH[infPrevH>1]=1  ;  infPrevHMat <-repmat(matrix(infPrevH, length(infPrevH), 1),1,nTest)  }
                 )

#-----------------------------Output--------------------------------------------------------------------------------- 
      list(rr=RiskRate,rrH=RiskRateHS, rrTrue=RiskRateTrue,rrTrueH=RiskRateTrueHS,standP=standPrev,rrInf=infPrevMat,rrInfH=infPrevHMat)
  }
#=====================================================================================================================
#====================================================================================================================
prTestEU25<-function(prTest,ageCap,ageCapEnd,prMeth) {
   if(prMeth==1) {    #mean of all years 
  HStp<-rowSums(HSPosProxy*prTest[4])/rowSums(HSPosProxy*prTest[4]+ESPosProxy*prTest[2])
  CStp<-rowSums(CSPosProxy*prTest[1])/rowSums(CSPosProxy*prTest[1]+FSPosProxy*prTest[3])
  }
  if(prMeth==2) { #mean of first 3 years
         HStp<-rowSums(HSPosProxy[,1:3]*prTest[4])/rowSums(HSPosProxy[,1:3]*prTest[4]+ESPosProxy[,1:3]*prTest[2])
         CStp<-rowSums(CSPosProxy[,1:3]*prTest[1])/rowSums(CSPosProxy[,1:3]*prTest[1]+FSPosProxy[,1:3]*prTest[3])
         }
 if(prMeth==3) {  #mean of last 4 years
         HStp<-rowSums(HSPosProxy[,7:10]*prTest[4])/rowSums(HSPosProxy[,7:10]*prTest[4]+ESPosProxy[,7:10]*prTest[2])
         CStp<-rowSums(CSPosProxy[,7:10]*prTest[1])/rowSums(CSPosProxy[,7:10]*prTest[1]+FSPosProxy[,7:10]*prTest[3])
         }
   HStp[is.na(HStp)]=HStp[find(!is.na(HStp))[1]]       #replace NaNs with value from youngest age group with a non-nan value
   CStp[is.na(CStp)]=CStp[find(!is.na(CStp))[1]]       #replace NaNs with value from youngest age group with a non-nan value
   EStp<-1-HStp
  FStp<-1-CStp
  
  HStp[0:ageCap[4]]=0                                                   # set number tested=0 for all animals under minimum test age
  if(!is.na(ageCapEnd[4])) HStp[ageCapEnd[4]:nrow(HStp)]=0        # set number tested=0 for all animals over maximum test age
 
  EStp[0:ageCap[2]]=0                                                   # set number tested=0 for all animals under minimum test age
  if(!is.na(ageCapEnd[2])) EStp[ageCapEnd[3]:nrow(EStp)]=0        # set number tested=0 for all animals over maximum test age
   
  CStp[0:ageCap[1]]=0                                                   # set number tested=0 for all animals under minimum test age
  if(!is.na(ageCapEnd[4])) CStp[ageCapEnd[4]:nrow(CStp)]=0        # set number tested=0 for all animals over maximum test age
 
  FStp[0:ageCap[3]]=0                                                   # set number tested=0 for all animals under minimum test age
  if(!is.na(ageCapEnd[2])) FStp[ageCapEnd[3]:nrow(FStp)]=0        # set number tested=0 for all animals over maximum test age

  list(HStp=HStp,CStp=CStp,EStp=EStp,FStp=FStp)
}

EstTrends<-function( aa,bb) { 
#------------------------ ESTIMATE TEST Positives ----------------------------------------#
  ###### CSFS estimated number test positive
  EstNumPosCSFS=aa$rr *bb$CSFSNumTest
  ## Split up to CS & FS - based on prop test positive
   pCSTestPos= bb$CSTestPos/(bb$CSFSTestPos)  #est proportion of CSFS test positives that are  CS => based on observed data
   pCSTestPos[is.na(pCSTestPos)]=( bb$CSNumTest/bb$CSFSNumTest)[is.na( pCSTestPos)]
   pCSTestPos[is.na(pCSTestPos)]=0
  
  pFSTestPos=1-pCSTestPos        #est proportion of CSFS that are FS => i.e. 1-CS
  EstNumPosCS=EstNumPosCSFS * pCSTestPos;
  EstNumPosFS=EstNumPosCSFS * pFSTestPos;
  ## HSES estimated number test positive
  EstNumPosHSES=aa$rrH * (bb$HSESNumTest)
  # split up to HS and ES
   pHSTestPos= bb$HSTestPos/(bb$HSESTestPos)  #est proportion of HSES test positives that are  HS => based on observed data
   pHSTestPos[is.na(pHSTestPos)]=( bb$HSNumTest/bb$HSESNumTest)[is.na( pHSTestPos)]
   pHSTestPos[is.na(pHSTestPos)]=0
pESTestPos=1-pHSTestPos        #est proportion of HSES that are ES => i.e. 1-HS
  EstNumPosHS=EstNumPosHSES * pHSTestPos;
  EstNumPosES=EstNumPosHSES * pESTestPos;
  
  EstNumPosESCSFS = EstNumPosES+EstNumPosCS+EstNumPosFS
  EstTotalNumPos= EstNumPosHSES+ EstNumPosCSFS
 
  #------------------------ ESTIMATE True infected -----------------------------# 
  EstNumPosTrueCSFS=aa$rrTrue *bb$CSFSNumTest  
    EstNumPosTrueCS=EstNumPosTrueCSFS * pCSTestPos;
    EstNumPosTrueFS=EstNumPosTrueCSFS * pFSTestPos;

  EstNumPosTrueHSES=aa$rrTrueH * (bb$HSESNumTest) 
    EstNumPosTrueHS=EstNumPosTrueHSES * pHSTestPos;
    EstNumPosTrueES=EstNumPosTrueHSES * pESTestPos;  
  EstNumPosTrue= EstNumPosTrueHSES+ EstNumPosTrueCSFS
#-----------------------------Output--------------------------------------------------------------------------------- 
  list(eCS=EstNumPosCS, eFS=EstNumPosFS, eHS=EstNumPosHS, eES=EstNumPosES, eCSFS=EstNumPosCSFS, eHSES=EstNumPosHSES,eESCSFS=EstNumPosESCSFS, eTot=EstTotalNumPos,eCSTrue=EstNumPosTrueCS,eFSTrue=EstNumPosTrueFS, eHSTrue=EstNumPosTrueHS, eESTrue=EstNumPosTrueES, eTotTrue=EstNumPosTrue)
}    #end EstTrends
################################################################################
################################################################################
EstTrendsDead<-function( aa,bb) { 
## Same method as EstTrends but using Number Dead inputs rather than Number tested
  ##--------------- ESTIMATE TEST Positives ---------------------------------###
  ###### CSFS estimated number test positive
  EstNumPosCSFS=aa$rr *bb$CSFSNumDead
  ## Split up to CS & FS - based on prop test positive
  pCSTestPos= bb$CSTestPos/(bb$CSFSTestPos)  #est proportion of CSFS test positives that are  CS => based on observed data
    pCSTestPos[is.na(pCSTestPos)]= pCSTestPosEU25[is.na(pCSTestPos)]
      pCSTestPos[is.na(pCSTestPos)]=0
  pFSTestPos=1-pCSTestPos        #est proportion of CSFS that are FS => i.e. 1-CS
  EstNumPosCS=EstNumPosCSFS * pCSTestPos;
  EstNumPosFS=EstNumPosCSFS * pFSTestPos;
  ## HSES estimated number test positive
  EstNumPosHSES=aa$rrH * (bb$HSESNumDead)
  # split up to HS and ES
  pHSTestPos= bb$HSTestPos/(bb$HSESTestPos)  #est proportion of HSES test positives that are  HS => based on observed data
    
    pHSTestPos[is.na(pHSTestPos)]= pHSTestPosEU25[is.na(pHSTestPos)]
      pHSTestPos[is.na(pHSTestPos)]=0
  pESTestPos=1-pHSTestPos        #est proportion of HSES that are ES => i.e. 1-HS
  EstNumPosHS=EstNumPosHSES * pHSTestPos;
  EstNumPosES=EstNumPosHSES * pESTestPos;
  
  EstNumPosESCSFS = EstNumPosES+EstNumPosCS+EstNumPosFS
  
  EstTotalNumPos= EstNumPosHSES+ EstNumPosCSFS
     
  #------------------------- ESTIMATE True clinically infected ----------------------------## 
  EstNumPosTrueCSFS=aa$rrTrue *bb$CSFSNumDead 
    EstNumPosTrueCS=EstNumPosTrueCSFS * pCSTestPos;
    EstNumPosTrueFS=EstNumPosTrueCSFS * pFSTestPos
  # Estimate True HSES 
  EstNumPosTrueHSES=aa$rrTrueH * (bb$HSESNumDead) 
    EstNumPosTrueHS=EstNumPosTrueHSES * pHSTestPos;
    EstNumPosTrueES=EstNumPosTrueHSES * pESTestPos;  
  EstNumPosTrue= EstNumPosTrueHSES+ EstNumPosTrueCSFS
#---------------------- Estimate true infected (assume independent of age) -----------------------------------------------
 EstNumPosInfCSFS=aa$rrInf *bb$CSFSNumDead 
    EstNumPosInfCS=EstNumPosInfCSFS * pCSTestPos;
    EstNumPosInfFS=EstNumPosInfCSFS * pFSTestPos
  # Estimate True HSES 
  EstNumPosInfHSES=aa$rrInfH * (bb$HSESNumDead) 
    EstNumPosInfHS=EstNumPosInfHSES * pHSTestPos;
    EstNumPosInfES=EstNumPosInfHSES * pESTestPos;  
  EstNumPosInf= EstNumPosInfHSES+ EstNumPosInfCSFS
  
#-----------------------------Output--------------------------------------------------------------------------------- 
  list(eCS=EstNumPosCS, eFS=EstNumPosFS, eHS=EstNumPosHS, eES=EstNumPosES, eCSFS=EstNumPosCSFS, eHSES=EstNumPosHSES,eESCSFS=EstNumPosESCSFS, eTot=EstTotalNumPos,eCSTrue=EstNumPosTrueCS,eFSTrue=EstNumPosTrueFS, eHSTrue=EstNumPosTrueHS, eESTrue=EstNumPosTrueES, eTotTrue=EstNumPosTrue,eHSESInf=EstNumPosInfHSES,eCSFSInf=EstNumPosInfCSFS,eHSInf=EstNumPosInfHS,eESInf=EstNumPosInfES,eCSInf=EstNumPosInfCS,eFSInf=EstNumPosInfFS,eTotInf=EstNumPosInf)
}    #end EstTrends
################################################################################
################################################################################
Emergence<-function(trendInc,trendIncTrue,trendUp) {
    trendInc1=trendInc[1]; i=1
    ##simulate until the number of cases in the current year (trendInc1) exceeds the threshold value (trendUp) or the number of years exceeds 'fp' 
    #(outputs only set to go up to fp - if you want to run for longer then change the value of fp in the BSE_parameters.R file)
    while (trendInc1<trendUp & i<fP) {   
           i=i+1
           trendInc[i]=trendInc[i-1]+incRate*trendInc[i-1]
           trendIncTrue[i]=trendIncTrue[i-1]+incRate*trendIncTrue[i-1]
           trendInc1=trendInc[i]
    }
     tFail<-i 
     nPos<-sum(trendInc[1:tFail])   
     nPosTrue<-sum(trendIncTrue[1:tFail])
     if(is.na(tFail)) tFail<-paste("detection has not occurred within",fP, "years")
     list(nPos=nPos,nPosTrue=nPosTrue,tFail=tFail)
 }