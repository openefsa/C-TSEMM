##########################################################################################
##### Inital Code to get input data into the right format to read into model #####
#########################################################################################

###Derive Age of animal from birth Cohort c at testing period t
jtmp=repmat(1:nTest,nCht,1)
ctmp=t(repmat(1:nCht,nTest,1))
ageMat=(jtmp+testStart)-(ctmp+chtStart)
ageMatNA=ageMat
ageMatNA[ageMatNA<0]=NA  #used in converting Age data to Cohort Data
ageMatNA[ageMatNA>maxAgeYr+1]=NA
ageMat[ageMat<0]=0  

#####  Set up FULL data (for true/infected prevlaence)#########################
CSNumTest=AgeConv(CSbyAge) ;ESNumTest=AgeConv(ESbyAge) ;FSNumTest=AgeConv(FSbyAge)   ;HSNumTest=AgeConv(HSbyAge)
ESNumDeadbyAge[is.na(ESNumDeadbyAge)]=0
 FSNumDeadbyAge[is.na(FSNumDeadbyAge)]=0
 HSNumDeadbyAge[is.na(HSNumDeadbyAge)]=0
ESNumDead=AgeConv(ESNumDeadbyAge);  FSNumDead=AgeConv(FSNumDeadbyAge);HSNumDead=AgeConv(HSNumDeadbyAge);
CSTestPos=AgeConv(CSbyAgePos)  ;ESTestPos=AgeConv(ESbyAgePos) ;FSTestPos=AgeConv(FSbyAgePos)   ;HSTestPos=AgeConv(HSbyAgePos)
################################################################################

####### Sums
CSFSNumTest=CSNumTest+FSNumTest
HSESNumTest=HSNumTest+ESNumTest

CSFSTestPos=CSTestPos+FSTestPos 
HSESTestPos=HSTestPos+ESTestPos 

ESCSFSTestPos=CSTestPos+FSTestPos+ESTestPos 

CSFSTestNeg=CSFSNumTest-CSFSTestPos; CSFSTestNeg[CSFSTestNeg<0]=0   #max to avoid getting negative numbers if inconsistency in data sets
HSESTestNeg=HSESNumTest-HSESTestPos; HSESTestNeg[HSESTestNeg<0]=0   #max to avoid getting negative numbers if inconsistency in data sets

TotTestPos=CSFSTestPos+HSESTestPos
TotNumTest=CSFSNumTest+HSESNumTest

##### calculate the required numerical integration  - these are used throughout and much quicker to do it once here than repeatedly every time it is needed
##put a 0 first for zero probability at age 0
   PclinAge=   c(0, sapply(1:ageMAX,function(ageA) quad(AgeOnsetDist, ageA-0.5,  (ageA+0.5),tol=1e-6,,Poly,knot,Onset ) )  )     # probability clinical onset at age a
   PdetectedAge= c(0, sapply(1:ageMAX,function(ageA) ifelse( (ageA+ageAdj)<=0,0,quad (AOTLC, ageA+ageAdj, (ageA+ageAdj + TestDetectionLength),tol=1e-6,, ageA+ageAdj,Poly,knot,Onset ))));   #Probability detection at age a, given test sensitivity
   ## Bound to prevent value for 0 being identically zero.
   PdetectedAge= c( sapply(0:ageMAX,function(ageA) ifelse( (ageA+ageAdj)<=-10,0,quad (AOTLC, ageA+ageAdj, (ageA+ageAdj + TestDetectionLength),tol=1e-6,, ageA+ageAdj,Poly,knot,Onset ))));   #Probability detection at age a, given test sensitivity
    PclinAge=   c(sapply(0:ageMAX,function(ageA) quad(AgeOnsetDist, ageA-0.5,  (ageA+0.5),tol=1e-6,,Poly,knot,Onset ) )  )     # probability clinical onset at age a

   PsurvAge= c(0, sapply(1:ageMAX,function(ageA) ifelse(ageA<=0,0,quad(AgeOnsetDist, 0 , ageA,tol=1e-6,,Poly,knot,Onset)))) ; #probability survival to age a

   PclinAgeL=   c(0, sapply(1:ageMAX,function(ageA) quad(AgeOnsetDist, ageA-0.5,  (ageA+0.5),tol=1e-6,,Poly,knot,OnsetL ) )  )     # probability clinical onset at age a
   PdetectedAgeL= c(0, sapply(1:ageMAX,function(ageA) ifelse( (ageA+ageAdj)<=0,0,quad (AOTLC, ageA+ageAdj, (ageA+ageAdj + TestDetectionLength),tol=1e-6,, ageA+ageAdj,Poly,knot,OnsetL ))));   #Probability detection at age a, given test sensitivity
 ## Bound to prevent value for 0 being identically zero.
  PclinAgeL=   c(sapply(0:ageMAX,function(ageA) quad(AgeOnsetDist, ageA-0.5,  (ageA+0.5),tol=1e-6,,Poly,knot,OnsetL ) )  )
  PdetectedAgeL= c( sapply(0:ageMAX,function(ageA) ifelse( (ageA+ageAdj)<=-10,0,quad (AOTLC, ageA+ageAdj, (ageA+ageAdj + TestDetectionLength),tol=1e-6,, ageA+ageAdj,Poly,knot,OnsetL ))));   #Probability detection at age a, given test sensitivity

   PsurvAgeL= c(0, sapply(1:ageMAX,function(ageA) ifelse(ageA<=0,0,quad(AgeOnsetDist, 0 , ageA,tol=1e-6,,Poly,knot,OnsetL)))) ; #probability survival to age a

   PclinAgeU=   c(0, sapply(1:ageMAX,function(ageA) quad(AgeOnsetDist, ageA-0.5,  (ageA+0.5),tol=1e-6,,Poly,knot,OnsetU ) )  )     # probability clinical onset at age a
   PdetectedAgeU= c(0, sapply(1:ageMAX,function(ageA) ifelse( (ageA+ageAdj)<=0,0,quad (AOTLC, ageA+ageAdj, (ageA+ageAdj + TestDetectionLength),tol=1e-6,, ageA+ageAdj,Poly,knot,OnsetU ))));   #Probability detection at age a, given test sensitivity
    ## Bound to prevent value for 0 being identically zero.
      PclinAgeU=   c(sapply(0:ageMAX,function(ageA) quad(AgeOnsetDist, ageA-0.5,  (ageA+0.5),tol=1e-6,,Poly,knot,OnsetU ) )  )  
     PdetectedAgeU= c(sapply(0:ageMAX,function(ageA) ifelse( (ageA+ageAdj)<=-10,0,quad (AOTLC, ageA+ageAdj, (ageA+ageAdj + TestDetectionLength),tol=1e-6,, ageA+ageAdj,Poly,knot,OnsetU ))));   #Probability detection at age a, given test sensitivity
   
   PsurvAgeU= c(0, sapply(1:ageMAX,function(ageA) ifelse(ageA<=0,0,quad(AgeOnsetDist, 0 , ageA,tol=1e-6,,Poly,knot,OnsetU)))) ; #probability survival to age a

   PC=matrix(0,nCht,nTest) ;     PCL=matrix(0,nCht,nTest) ;   PCU=matrix(0,nCht,nTest) 
   PD= matrix(0,nCht,nTest) ; PDL= matrix(0,nCht,nTest); PDU= matrix(0,nCht,nTest)
   PS=matrix(0,nCht,nTest) ; PSU=matrix(0,nCht,nTest); PSL=matrix(0,nCht,nTest)
   for (cht in 1:nCht) { 
      PC[cht,]=PclinAge[ageMat[cht,]+1]   ;     PCL[cht,]=PclinAgeL[ageMat[cht,]+1];     PCU[cht,]=PclinAgeU[ageMat[cht,]+1]
      PD[cht,]=PdetectedAge[ageMat[cht,]+1] ;   PDL[cht,]=PdetectedAgeL[ageMat[cht,]+1]; PDU[cht,]=PdetectedAgeU[ageMat[cht,]+1]
      PS[cht,]=PsurvAge[ageMat[cht,]+1] ;       PSL[cht,]=PsurvAgeL[ageMat[cht,]+1];     PSU[cht,]=PsurvAgeU[ageMat[cht,]+1]
   }
   
   #######Set up Full regime data   ############################################
   ageCapT=matrix(0,4)
   ageCapTend=matrix(NA,4)
   prTestT<-matrix(1,4)
   
   #######Set Up Adult testing Data  ########################################################
   # age above which animals are tested
   ageCapA=matrix(1,4)
   ageCapStrA=c(">24",">24",">24",">24")
   ageCapAend=matrix(NA,4)
  prTestA=matrix(1,4)
  
   ######Set Up Baseline regime Data  ########################################################
   # age above which animals are tested
   ageCapB=matrix(0,4)
     ageCapB[1]=pmatch(sub(">","",as.character(UserInput[1,37])),rownames(CSbyAge))-1 
     ageCapB[2]=  pmatch(sub(">","",as.character(UserInput[1,35])),rownames(ESbyAge))-1 
     ageCapB[3]=  pmatch(sub(">","",as.character(UserInput[1,36])),rownames(FSbyAge))-1 
     ageCapB[4]=  pmatch(sub(">","",as.character(UserInput[1,34])),rownames(HSbyAge))-1 
   ageCapStrB=c(as.character(UserInput[1,37]),as.character(UserInput[1,35]),as.character(UserInput[1,36]),as.character(UserInput[1,34]))
   ageCapBend=matrix(0,4)
     ageCapBend[1]=pmatch(sub(">","",as.character(UserInput[1,45])),rownames(CSbyAge))-1 
     ageCapBend[2]=  pmatch(sub(">","",as.character(UserInput[1,43])),rownames(ESbyAge))-1 
     ageCapBend[3]=  pmatch(sub(">","",as.character(UserInput[1,44])),rownames(FSbyAge))-1 
     ageCapBend[4]=  pmatch(sub(">","",as.character(UserInput[1,42])),rownames(HSbyAge))-1 
  
  #proportion of animals tested above the ageCap 
  prTestB=matrix(0,4)
  prTestB[1]=UserInput[1,41] ;  prTestB[2]=UserInput[1,39] ;  prTestB[3]=UserInput[1,40] ;  prTestB[4]=UserInput[1,38] 
   
 ################################################################################
 #### Set up Scenario data ########################################################
  ageCapS=matrix(0,4)
  ageCapS[1]= pmatch(sub(">","",as.character(UserInput[1,19])),rownames(CSbyAge))-1    #max age at which animals are not tested - in terms of row number of matrix: 1=<24, 2=24-29 etc... need to -1 as pmatch identifies lower end of age band - these animals will be tested, all below will not. Using sub as need to remove the ">" sign from the User Log   
    if(is.na(ageCapS[1])) {ageCapS[1]=0 ;  }   #if get bad value then assume ageCap doesn't exist
  ageCapS[2]= pmatch(sub(">","",as.character(UserInput[1,17])),rownames(ESbyAge))-1     
    if(is.na(ageCapS[2])) {ageCapS[2]=0 ; }
  ageCapS[3]= pmatch(sub(">","",as.character(UserInput[1,18])),rownames(FSbyAge))-1 
    if(is.na(ageCapS[3])) {ageCapS[3]=0 ;   }    
  ageCapS[4]= pmatch(sub(">","",as.character(UserInput[1,16])),rownames(HSbyAge))-1   
      if(is.na(ageCapS[4])) {ageCapS[4]=0 ;  } 
  
  ageCapSend=matrix(0,4)
    ageCapSend[1]= pmatch(sub("<","",as.character(UserInput[1,12])),rownames(CSbyAge))-1    #max age at which animals are not tested - in terms of row number of matrix: 1=<24, 2=24-29 etc... need to -1 as pmatch identifies lower end of age band - these animals will be tested, all below will not. Using sub as need to remove the ">" sign from the User Log   
    if(is.na(ageCapSend[1])) {ageCapSend[1]=NA ;  }   #if get bad value then assume ageCap doesn't exist
  ageCapSend[2]= pmatch(sub("<","",as.character(UserInput[1,10])),rownames(ESbyAge))-1     
    if(is.na(ageCapSend[2])) {ageCapSend[2]=NA ; }
  ageCapSend[3]= pmatch(sub("<","",as.character(UserInput[1,11])),rownames(FSbyAge))-1 
    if(is.na(ageCapSend[3])) {ageCapSend[3]=NA ; }    
  ageCapSend[4]= pmatch(sub("<","",as.character(UserInput[1,9])),rownames(HSbyAge))-1   
      if(is.na(ageCapSend[4])) {ageCapSend[4]=NA ;  } 
 
  ageCapStrS=c( as.character(UserInput[1,19]),as.character(UserInput[1,17]),as.character(UserInput[1,18]),as.character(UserInput[1,16]))
  ## Proportion of animals tested
  prTestS=matrix(0,4)
  prTestS[1]=UserInput[1,23]  ;  prTestS[2]=UserInput[1,21] ;  prTestS[3]=UserInput[1,22] ;  prTestS[4]=UserInput[1,20] 
  
  