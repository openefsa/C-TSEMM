################################################################################
## Number to Test to achieve design prevalence, dP (default 1 in 100,000)  #####
################################################################################

designPrev<-function(aa,bb,bb2,cc,ageCap,prTest,dPType,standEst) {
  baseTP= dim(bb$eHSTrue)[2]                    #column of matrices which contain the basline testing period
  TotNumTest=aIn$HSESNumTest+aIn$CSFSNumTest    # total number of animals tested
  tP=(colSums(cc$eHSTrue+cc$eESTrue+cc$eCSTrue+cc$eFSTrue)/colSums(TotNumTest)) [baseTP] #true prevalence
  ifelse(standEst==1,piTot1<-colSums(baseR$standP)[baseTP]/colSums(standPop), piTot1<- (colSums(cc$eTot)/colSums(TotNumTest))[baseTP] ) # detection prevalence
         if (colSums(TotNumTest)[baseTP]==0) piTot=0;
  #### model predicted stream prevalences
  piCS1=(colSums(bb$eCS)/colSums(aa$CSNumTest) )[baseTP]; 
        if (colSums(aa$CSNumTest)[baseTP]==0) piCS1=0;
  piES1=(colSums(bb$eES)/colSums(aa$ESNumTest))[baseTP] ;
        if (colSums(aa$ESNumTest)[baseTP]==0) piES1=0;
  piFS1=(colSums(bb$eFS)/colSums(aa$FSNumTest))[baseTP]  ;
        if (colSums(aa$FSNumTest)[baseTP]==0) piFS1=0;
  piHS1=(colSums(bb$eHS)/colSums(aa$HSNumTest))[baseTP]  ;       
        if (colSums(aa$HSNumTest)[baseTP]==0) piHS1=0;      
  ## set dC=0 if tP or piTot1=0 or model crashes    
  ifelse(dPType=="inf", {ifelse(tP==0,dC<-0,dC<-dP/tP)},{ifelse(piTot1==0,dC<-0,dC<-dP/piTot1)})   #choose design conversion factor for either infection or detection total prevalence           
  #### scale model predicted prevalence to design prevalence
  piES= (piES1)*dC
  piFS=(piFS1)*dC
  piHS=(piHS1)*dC
  piTot=piTot1*dC
  #For clinical suspects need to scale number tested (and test positive) rather than prevalence
  #because if we are reducing/increasing the number of BSE infected animals then this impacts the number tested
  #(if an animal is no longer infected then it wont have clinical signs and so won't be tested)
  #Not doing this could result in a CS prevalence >1 if we are scaling up
  piCS=(piCS1);  CSNumTestDp2=aa$CSNumTest*dC;  
  
   #number tested by stream
   CSnt=colSums(CSNumTestDp2)[baseTP]
   ESnt= colSums(aa$ESNumTest)[baseTP]
   FSnt=colSums(aa$FSNumTest)[baseTP]
   HSnt=colSums(aa$HSNumTest)[baseTP]
   Totnt= colSums(TotNumTest)[baseTP]
   #number tested if testing in other streams
  ndpHS=    (log(1-tau)- log(1-piCS)*(CSnt)-log(1-piES)*(ESnt)-log(1-piFS)*( FSnt))/log(1-piHS)
  if(piHS==0) ndpHS=0;  
  ndpCS=    (log(1-tau)- log(1-piHS)*( HSnt)-log(1-piES)*(ESnt)-log(1-piFS)*( FSnt))/log(1-piCS)
  if(piCS==0) ndpCS=0;  
  ndpES=    (log(1-tau)- log(1-piCS)*(CSnt)-log(1-piHS)*( HSnt)-log(1-piFS)*( FSnt))/log(1-piES)
  if(piES==0) ndpES=0;  
  ndpFS=    (log(1-tau)- log(1-piCS)*(CSnt)-log(1-piES)*(ESnt)-log(1-piHS)*( HSnt))/log(1-piFS)
  if(piFS==0) ndpFS=0;  
  ndpTot=  log(1-tau)/log(1-piTot)
   #number tested if not testing in other streams
  ndpHS0=    (log(1-tau))/log(1-piHS)
  ndpCS0=    (log(1-tau))/log(1-piCS)
  ndpES0=    (log(1-tau))/log(1-piES)
  ndpFS0=    (log(1-tau))/log(1-piFS)
   #if dC set to 0 then set outputs to Inf
  if (dC==0) {ndpHS=Inf;  ndpES=Inf;ndpFS=Inf;ndpCS=Inf;ndpHS0=Inf;ndpES0=Inf;ndpFS0=Inf;ndpCS0=Inf;}
   ##Outputs   
    nToTest=matrix(0,7,5)
  colnames(nToTest)=c("HS","ES","FS","CS","Overall")
  nToTest[1,]= c(paste(ageCap[4],",",prTest[4],",",bb2,",",dPType), paste(ageCap[2],",",prTest[2],",",bb2,",",dPType),paste(ageCap[3],",",prTest[3],",",bb2,",",dPType),paste(ageCap[1],",",prTest[1],",",bb2,",",dPType),"NA")
  nToTest[2,]=round(c(1/piHS1,1/piES1,1/piFS1,1/piCS1,1/piTot1),digits=0)
  nToTest[3,]=round(c(1/piHS,1/piES,1/piFS,1/piCS,1/piTot),digits=0)
  nToTest[4,]=round(c(HSnt,ESnt,FSnt,CSnt,Totnt),digits=0)
  nToTest[5,]=round(c(max(0,ndpHS),max(0,ndpES),max(0,ndpFS),max(0,ndpCS),max(0,ndpTot)),digits=0)
  nToTest[6,]=round(c(max(0,ndpHS/HSnt),max(0,ndpES/ESnt),max(0,ndpFS/FSnt),max(0,ndpCS/CSnt),max(0,ndpTot/Totnt)),digits=2)
  nToTest[7,]=round(c(ndpHS0,ndpES0,ndpFS0,ndpCS0,ndpTot),digits=0)
 
  rownames( nToTest)=c("Monitoring Scenario", "Estimated observed prevalence in exit stream (expressed as: 1 in X)","Estimated scaled prevalence in exit stream (expressed as: 1 in X)", "Actual number tested","Number tested to detect design prevalence if other streams are tested as prescribed by user","Proportion of actual number tested","Number to test if only this stream tested")#,"detectable prevalence")
if (UserInput[1,8]==TRUE)  write.csv(nToTest[,1:4],paste(outDir,"/",MSName,"NumberToTest(",bb2,"_",dPType,").csv",sep=""))       #not outputting overall column now.
  list(nToTest=nToTest)
}     #end DesignPrev

################################################################################
##### 'Design Prevalence' of Current system ####################################
################################################################################

nDetect<-function(aa,bb,bb2,cc,dPType,standEst)   {
ifelse(colSums(cc$eTot)==0,{scenarioDp<-0;dpSolver(0,1,aa,bb,cc,dPType,standEst)} ,{
  tmp<-nlminb(dpSolve0,dpSolver,gradient=NULL,hessian=NULL,0,aa,bb,cc,dPType,standEst,control=list(eval.max=1000),upper=upDpLim)        
  ifelse(tmp$par==0,scenarioDp<-0,scenarioDp<-1/tmp$par) 
  dpSolve1=dpSolve0;  #initial value for parameter to change initial starting value if it fails.
  #if estimated Dp is <= upper limit then assume hasn't converged: +1 to account for rounding errors. 
  #Change initial start value of solver routine and try again.
  # If still no good value when inital value gets to 1e-10 then stop. 
  while(scenarioDp<=(1/upDpLim+1) && dpSolve1>1e-9) { 
         dpSolve1= 10^(log10(dpSolve1)-1)   #i.e. reduce initial starting value by 1 log10 (e.g. 1e-5 -> 1e-6)
         tmp<-nlminb(dpSolve1,dpSolver,gradient=NULL,hessian=NULL,0,aa,bb,cc,dPType,standEst,control=list(eval.max=1000),upper=upDpLim)
         ifelse(tmp$par==0,scenarioDp<-0,scenarioDp<-1/tmp$par) 
         }     
   dpSolver(tmp$par,1,aa,bb,cc,dPType,standEst)
   } )
   fixedScenario[q,1]=MS          #member state id
   fixedScenario[q,2:4]=numMissedOut[1:3,nTest]

   dpOutS=read.csv(paste(wd,"/dpOutS.csv",sep=""),header=TRUE,row.names=1) 
   ifelse( abs( (sum((dpOutS)[2,]) - sum((dpOutS)[1,]))/sum((dpOutS)[1,])) <=0.1, fixedScenario[q,5]<-formatC(scenarioDp,digits=0,big.mark=",",format="f"),fixedScenario[q,5]<-"N/A")
   fixedScenario[q,6]= formatC((dpOutS)[1,1],digits=0,big.mark=",",format="f")
   fixedScenario[q,7]= formatC((dpOutS)[2,1],digits=0,big.mark=",",format="f")
   colnames(fixedScenario)=c("MS","Baseline test positive", "Monitoring scenario test positive", "Number missed by monitoring scenario", "baseline 'design prevalence' (expressed as 1 in X)","observed number tested","number tested to acheive dp")
 if (UserInput[1,13]==TRUE)  write.csv((fixedScenario[q,5:7]),paste(outDir,"/",MSName,"BaselineDp(",bb2,"_",dPType,").csv",sep="")) 
   list(nDetOut=fixedScenario)
 }
################################################################################
## Function to find value of dP at which ndpH~=observed number tested for all streams
## very similar setup to designPrev function 
################################################################################
dpSolver<-function(dP,writeOut,aa,bb,cc,dPType,standEst) {
  TotNumTest=aIn$HSESNumTest+aIn$CSFSNumTest
  baseTP= dim(bb$eHSTrue)[2]  
  tP=(colSums(cc$eHSTrue+cc$eESTrue+cc$eCSTrue+cc$eFSTrue)/colSums(TotNumTest)) [baseTP] #true prevalence
   ifelse(standEst==1,piTot1<-colSums(baseR$standP)[baseTP]/colSums(standPop), piTot1<- (colSums(cc$eTot)/colSums(TotNumTest))[baseTP] ) # detection prevalence
  
        if (colSums(TotNumTest)[baseTP]==0) piTot=0;
        
   #### model predicted stream prevalence
  piCS1=(colSums(bb$eCS)/colSums(aa$CSNumTest) )[baseTP]; 
        if (colSums(aa$CSNumTest)[baseTP]==0) piCS1=0;
  piES1=(colSums(bb$eES)/colSums(aa$ESNumTest))[baseTP] ;
        if (colSums(aa$ESNumTest)[baseTP]==0) piES1=0;
  piFS1=(colSums(bb$eFS)/colSums(aa$FSNumTest))[baseTP]  ; 
        if (colSums(aa$FSNumTest)[baseTP]==0) piFS1=0;
  piHS1=(colSums(bb$eHS)/colSums(aa$HSNumTest))[baseTP]  ;
        if (colSums(aa$HSNumTest)[baseTP]==0) piHS1=0;
    
  ## set dC=0 if tP or piTot1=0 or model crashes    
    ifelse(dPType=="inf", {ifelse(tP==0,dC<-0,dC<-dP/tP)},{ifelse(piTot1==0,dC<-0,dC<-dP/piTot1)})   #choose design conversion factor for either infection or detection total prevalence           
 #### scale model predicted prevalence to design prevalence
  piES= (piES1)*dC
  piFS=(piFS1)*dC
  piHS=(piHS1)*dC
  piTot=piTot1*dC
   piCS=(piCS1);  CSNumTestDp2=aa$CSNumTest*dC;   
   if(is.na(piHS)) piHS=0   ;       if(is.na(piES)) piES=0   ;   if(is.na(piFS)) piFS=0   ;   if(is.na(piCS)) piCS=0   ;   if(is.na(piTot)) piTot=0   ;
  #pi's can't exceed tau or will break, so set equal to tau if this happens (shouldn't do usually but prevents crash if it does)    
  if(piHS>tau) piHS=tau; if(piES>tau) piES=tau; if(piFS>tau) piFS=tau;     if(piCS>tau) piCS=tau

   CSnt=colSums(CSNumTestDp2)[baseTP];  
  ESnt= colSums(aa$ESNumTest)[baseTP]
   FSnt=colSums(aa$FSNumTest)[baseTP]
   HSnt=colSums(aa$HSNumTest)[baseTP]
   Totnt= colSums(TotNumTest)[baseTP]
  ndpHS=    (log(1-tau)- log(1-piCS)*(CSnt)-log(1-piES)*(ESnt)-log(1-piFS)*( FSnt))/log(1-piHS)
            if(piHS==0) ndpHS=0;   #set=0 if piHS=0 or get -ve logs.  Also not good if ndpHS<0 so penalise heavily - PL wasn't converging without this...
  ndpCS=    (log(1-tau)- log(1-piHS)*( HSnt)-log(1-piES)*(ESnt)-log(1-piFS)*( FSnt))/log(1-piCS)
            if(piCS==0) ndpCS=0;  
  ndpES=    (log(1-tau)- log(1-piCS)*(CSnt)-log(1-piHS)*( HSnt)-log(1-piFS)*( FSnt))/log(1-piES)
            if(piES==0) ndpES=0;   
  ndpFS=    (log(1-tau)- log(1-piCS)*(CSnt)-log(1-piES)*(ESnt)-log(1-piHS)*( HSnt))/log(1-piFS)
            if(piFS==0) ndpFS=0;    
  ndpTot=  log(1-tau)/log(1-piTot)
            if(piTot==0) ndpTot=0;
  if (writeOut==1) { 
        #if dC set to 0 then set outputs to Inf
        if (dC==0) {ndpHS=Inf;  ndpES=Inf;ndpFS=Inf;ndpCS=Inf;ndpHS0=Inf;ndpES0=Inf;ndpFS0=Inf;ndpCS0=Inf;}                                           
        dpOutS=matrix(0,2,4)
        dpOutS[1,]=c(colSums(aa$HSNumTest)[baseTP],colSums(aa$ESNumTest)[baseTP],colSums(aa$FSNumTest)[baseTP],colSums(CSNumTestDp2)[baseTP] )
        dpOutS[2,]=c(ndpHS,ndpES,ndpFS,ndpCS)
        rownames(dpOutS)=c("observed number tested (CS scaled by dP)","number tested to achieve dp")
        colnames(dpOutS)=c("HS","ES","FS","CS")
        write.csv(dpOutS,paste(wd,"/dpOutS.csv",sep="")) 
     }
   # ansDp is function we are trying to minimise
  ansDp=((colSums(aa$HSNumTest)[baseTP]-ndpHS)^2 + (colSums(CSNumTestDp2)[baseTP]-ndpCS)^2 +(colSums(aa$FSNumTest)[baseTP]-ndpFS)^2+(colSums(aa$ESNumTest)[baseTP]-ndpES)^2)
  }
    
                                                                                          