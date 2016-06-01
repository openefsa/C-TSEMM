################################################################################
## EFSA BSE Surveillance model: BSE_Parameters.r 
################################################################################

### LOAD R LIBRARIES ###########################################################
library(pracma)
library(splines)
#### Setup #####################################################################
runPoisInCI=0;     #set=1 to run poisson confidence intervals on A & B estimates - default off as takes v long time and not necessary.
useOnsetCI=FALSE   #use age of onset confidence intervals  - false uses poisson about input values

###########Baseline Testing Regieme
trAgeHS=  UserInput[34]
trAgeES=  UserInput[35]
trAgeFS= UserInput[36]            
trAgeCS= UserInput[37]

########################################################################################################
### READ IN THE DATA
### In order to read in the data correctly, the name of the file must match with the text file and
### the number of columns must be the same (i.e. ncol = X; X = number of columns).
#########################################################################################################
ifelse(runSeparate=="FALSE", {                                                ##If running model for all MS's then already defined MS so don't overwrite here.
          MSName=   as.character(UserInput[1,3]);                             #full MS's - include all MS's for use in output directory name
          MultiSt=unlist(strsplit(as.character(UserInput[1,3]), split="-"))   #split string to find mutiple MS's        
          nMS= length(MultiSt) },                                             #not doing multiple MS's so number MS's=1. 
          { MSName=MS; nMS=1;MultiSt="NA" } 
       )
#-------------------------------------------------------------------------------
for (i in 1:nMS)    {
       if(runSeparate=="FALSE") { MS=MultiSt[i]} 
    ##--------------------Poison confidence limits--------------------------------
    poU=0.975; poL=0.025    #Upper and lower values for input Poisson Confidence intervals
    poUu=0.999;             #upper limit for Emergence
    poUl=0.9;               #lower limit for Emergence                           
    ##------------------------- Estimating Trend  ----------------------------------
      estTrend=1;                                               #1 for country estimation of trend, 2 for EU trend and 3 for UserInput trend  
      paras=c(0.6,-2,2)                                         # inital parameter estimates for A,B & K
      parasWbl<-c(3.0069329 ,1.2288016,1.3177901,-0.3863629)    #inital estimates for Weibull distribution
      ifelse(runSeparate=="TRUE",EUTrendPos<- pmatch(MS,MSnoTrend[,1]), EUTrendPos<- max(pmatch(MultiSt,MSnoTrend[,1])))          #whether to use EU data or not
      EUNum="NA";           #initial value for this variable
      ifelse( is.na(EUTrendPos),{EUNum="NA"},{estTrend=2 ; EUNum=as.character(MSnoTrend[EUTrendPos,2])} )
      if(EUNum=="EU8")  EU8skip=as.numeric(read.csv(paste(wd,"/","EU8skip.csv",sep=""),row.names=1))  
      if (estTrend<1 || estTrend>2) stop("value of estTrend is invalid")
        ## ------------------------ Number Dead --------------------------------  
    setwd(paste(fileInput,"Number Dead",sep="" ) )
       ifelse(i==1, {ESNumDeadbyAge <-(read.csv(paste(MS,"- ES Number Dead.csv"),header=TRUE,row.names=1 )   )
                    FSNumDeadbyAge <-(read.csv(paste(MS,"- FS Number Dead.csv"),header=TRUE,row.names=1 )   )
                    HSNumDeadbyAge <-(read.csv(paste(MS,"- HS Number Dead.csv"),header=TRUE,row.names=1 )   )
                    }  ,  {   # multi MS run so add to previous input
                    ESNumDeadbyAge <- ESNumDeadbyAge+(read.csv(paste(MS,"- ES Number Dead.csv"),header=TRUE,row.names=1 ) ) 
                    FSNumDeadbyAge <- FSNumDeadbyAge+(read.csv(paste(MS,"- FS Number Dead.csv"),header=TRUE,row.names=1 ) )
                    HSNumDeadbyAge <- HSNumDeadbyAge+(read.csv(paste(MS,"- HS Number Dead.csv"),header=TRUE,row.names=1 ) )
                    }
             )  
    ## ------------------------ Standing Population --------------------------------  
    # for extra outputs.
    setwd(paste(fileInput,"Standing Population",sep="" ) )
       ifelse(i==1, {standPop <-(read.csv(paste(MS,".csv",sep=""),header=FALSE )   )
                    }  ,  {   # multi MS run so add to previous input
                   standPop <- standPop+(read.csv(paste(MS,".csv",sep=""),header=FALSE ) ) 
                    }
             )           
    ## -------------------- Number Tested ------------------------------------------      
    #### Read in data dependent on testing periods
    setwd(paste(fileInput,"Number Tested",sep="" ) ) 
      TestStartIn <-  (read.csv(paste(MS,"- CS Number Tested.csv") ,row.names=1,check.names=FALSE   )   )     #to get numeric column headers, first of which is the first testing period - R doesn't like numeric column headings so not doing this for other files  
      CSbyAgeTmp <-as.matrix(read.csv(paste(MS,"- CS Number Tested.csv"),row.names=1 )  ) 
      ESbyAgeTmp <-as.matrix(read.csv(paste(MS,"- ES Number Tested.csv"),row.names=1 )  ) 
      FSbyAgeTmp <-as.matrix(read.csv(paste(MS,"- FS Number Tested.csv"),row.names=1 )  )
      HSbyAgeTmp <-as.matrix(read.csv(paste(MS,"- HS Number Tested.csv"),row.names=1 )  ) 
        
       ###if multi MS run & i>1, add to previous input
      ifelse(i==1, CSbyAge <- CSbyAgeTmp,  CSbyAge <- CSbyAge+CSbyAgeTmp  )    
      ifelse(i==1, ESbyAge <- ESbyAgeTmp,  ESbyAge <- ESbyAge+ESbyAgeTmp  )
      ifelse(i==1, FSbyAge <- FSbyAgeTmp,  FSbyAge <- FSbyAge+FSbyAgeTmp  )
      ifelse(i==1, HSbyAge <- HSbyAgeTmp,  HSbyAge <- HSbyAge+HSbyAgeTmp  )
    ## ------------------------ Test Positives for each strain type -------------------------------------
    setwd(paste(fileInput,"Test Positive",sep="" ) )  
     ifelse(i==1, CSbyAgePosC <-(read.csv(paste(MS,"- CS C Test Positive.csv"),header=FALSE,row.names=1  ) ), 
                  CSbyAgePosC <-CSbyAgePosC+(read.csv(paste(MS,"- CS C Test Positive.csv"),header=FALSE,row.names=1  ) )  )
     ifelse(i==1, CSbyAgePosH <-(read.csv(paste(MS,"- CS H Test Positive.csv"),header=FALSE,row.names=1    ) ), 
                  CSbyAgePosH <-CSbyAgePosH+(read.csv(paste(MS,"- CS H Test Positive.csv"),header=FALSE,row.names=1    ) ) ) 
     ifelse(i==1, CSbyAgePosL <-(read.csv(paste(MS,"- CS L Test Positive.csv"),header=FALSE,row.names=1     ) ),  
                  CSbyAgePosL <-CSbyAgePosL+(read.csv(paste(MS,"- CS L Test Positive.csv"),header=FALSE,row.names=1    ) ) ) 
     ifelse(i==1, CSbyAgePosU <-(read.csv(paste(MS,"- CS U Test Positive.csv"),header=FALSE ,row.names=1   ) ),  
                  CSbyAgePosU <-CSbyAgePosU+(read.csv(paste(MS,"- CS U Test Positive.csv"),header=FALSE,row.names=1    ) ) ) 
      
     ifelse(i==1, ESbyAgePosC <-(read.csv(paste(MS,"- ES C Test Positive.csv"),header=FALSE  ,row.names=1  ) ),
                  ESbyAgePosC <-ESbyAgePosC+(read.csv(paste(MS,"- ES C Test Positive.csv"),header=FALSE  ,row.names=1  ) ))
     ifelse(i==1, ESbyAgePosH <-(read.csv(paste(MS,"- ES H Test Positive.csv"),header=FALSE ,row.names=1   ) ), 
                  ESbyAgePosH <-ESbyAgePosH+(read.csv(paste(MS,"- ES H Test Positive.csv"),header=FALSE  ,row.names=1  ) ))
     ifelse(i==1, ESbyAgePosL <-(read.csv(paste(MS,"- ES L Test Positive.csv") ,header=FALSE ,row.names=1   ) ), 
                  ESbyAgePosL <-ESbyAgePosL+(read.csv(paste(MS,"- ES L Test Positive.csv"),header=FALSE  ,row.names=1  ) ))
     ifelse(i==1, ESbyAgePosU <-(read.csv(paste(MS,"- ES U Test Positive.csv") ,header=FALSE ,row.names=1   ) ), 
                  ESbyAgePosU <-ESbyAgePosU+(read.csv(paste(MS,"- ES U Test Positive.csv"),header=FALSE  ,row.names=1  ) ))   
       
     ifelse(i==1, FSbyAgePosC <-(read.csv(paste(MS,"- FS C Test Positive.csv"),header=FALSE  ,row.names=1   ) ), 
                  FSbyAgePosC <-FSbyAgePosC+(read.csv(paste(MS,"- FS C Test Positive.csv"),header=FALSE  ,row.names=1   ) )) 
     ifelse(i==1, FSbyAgePosH <-(read.csv(paste(MS,"- FS H Test Positive.csv"),header=FALSE  ,row.names=1   ) ),
                  FSbyAgePosH <-FSbyAgePosH+(read.csv(paste(MS,"- FS H Test Positive.csv"),header=FALSE  ,row.names=1   ) )) 
     ifelse(i==1, FSbyAgePosL <-(read.csv(paste(MS,"- FS L Test Positive.csv") ,header=FALSE  ,row.names=1  ) ),
                  FSbyAgePosL <-FSbyAgePosL+(read.csv(paste(MS,"- FS L Test Positive.csv"),header=FALSE  ,row.names=1   ) )) 
     ifelse(i==1, FSbyAgePosU <-(read.csv(paste(MS,"- FS U Test Positive.csv") ,header=FALSE ,row.names=1    ) ),
                  FSbyAgePosU <-FSbyAgePosU+(read.csv(paste(MS,"- FS U Test Positive.csv"),header=FALSE  ,row.names=1   ) )) 
           
     ifelse(i==1, HSbyAgePosC <-(read.csv(paste(MS,"- HS C Test Positive.csv"),row.names=1    )   )              ,
                  HSbyAgePosC <-HSbyAgePosC+(read.csv(paste(MS,"- HS C Test Positive.csv"),row.names=1    )   )) 
     ifelse(i==1, HSbyAgePosH <-(read.csv(paste(MS,"- HS H Test Positive.csv")  ,header=FALSE ,row.names=1  )   ),
                  HSbyAgePosH <-HSbyAgePosH+(read.csv(paste(MS,"- HS H Test Positive.csv")  ,header=FALSE ,row.names=1  )   ) )
     ifelse(i==1, HSbyAgePosL <-(read.csv(paste(MS,"- HS L Test Positive.csv")  ,header=FALSE ,row.names=1  )   ),
                  HSbyAgePosL<-HSbyAgePosL+ (read.csv(paste(MS,"- HS L Test Positive.csv")  ,header=FALSE ,row.names=1  )   ))
     ifelse(i==1, HSbyAgePosU <-(read.csv(paste(MS,"- HS U Test Positive.csv")  ,header=FALSE ,row.names=1  )   ),
                  HSbyAgePosU <-HSbyAgePosU+(read.csv(paste(MS,"- HS U Test Positive.csv")  ,header=FALSE ,row.names=1  )   ) )
 }     #end for (i in 1:nMS)
#-------------------------------------------------------------------------------
standPop[1:2,]=0   #set animals <24 months=0 (i.e. from birth cohorts of last 2 years)

#Sum up to get baseline Test Positive Files for required strain type(s):   
 strainType=unlist(strsplit(as.character(UserInput[1,31]),"\""))[1]
 if (strainType=="CLASSICAL AND UNKNOWN") st="CU";
 if (strainType=="TYPE H")   st="H"
 if (strainType=="TYPE L")   st="L"
 if (strainType=="TYPE H AND L") st="HL"
 if (strainType=="C H L AND UNKNOWN") st="CUHL"
    
 if(estTrend==1) {
     if (strainType=="CLASSICAL AND UNKNOWN") {  
        CSbyAgePos<- CSbyAgePosC+CSbyAgePosU ; ESbyAgePos<- ESbyAgePosC+ESbyAgePosU
        FSbyAgePos<- FSbyAgePosC+FSbyAgePosU ; HSbyAgePos<- HSbyAgePosC+HSbyAgePosU
               }  #end if
      if (strainType=="TYPE H") { 
        CSbyAgePos<- CSbyAgePosH ; ESbyAgePos<- ESbyAgePosH
        FSbyAgePos<- FSbyAgePosH ; HSbyAgePos<- HSbyAgePosH
               }    #end if
      if (strainType=="TYPE L") { 
        CSbyAgePos<- CSbyAgePosL ; ESbyAgePos<- ESbyAgePosL
        FSbyAgePos<- FSbyAgePosL ; HSbyAgePos<- HSbyAgePosL
               }    #end if
     if (strainType=="TYPE H AND L") { 
        CSbyAgePos<- CSbyAgePosH+CSbyAgePosL ; ESbyAgePos<- ESbyAgePosH+ESbyAgePosL
        FSbyAgePos<- FSbyAgePosH+FSbyAgePosL ; HSbyAgePos<- HSbyAgePosH+HSbyAgePosL
               }     #end if
     if (strainType=="C H L AND UNKNOWN") { 
        CSbyAgePos<- CSbyAgePosC+CSbyAgePosU+ CSbyAgePosH+CSbyAgePosL ; 
        ESbyAgePos<- ESbyAgePosC+ESbyAgePosU+ ESbyAgePosH+ESbyAgePosL
        FSbyAgePos<- FSbyAgePosC+FSbyAgePosU+ FSbyAgePosH+FSbyAgePosL ; 
        HSbyAgePos<- HSbyAgePosC+HSbyAgePosU+ HSbyAgePosH+HSbyAgePosL
               }      #end if
     }  # end if (estTrend==1)
 ####If using EU trend data then just input zero matrices - some input files have no data and can crash if try to do estTrend==1 loop           
 if (estTrend==2) {
    CSbyAgePos<-  matrix(0,dim(CSbyAge)[1],dim(CSbyAge)[2]); rownames(CSbyAgePos)<-rownames(CSbyAge)  ; colnames(CSbyAgePos)<-colnames(CSbyAge)
    ESbyAgePos<-  matrix(0,dim(CSbyAge)[1],dim(CSbyAge)[2]); rownames(CSbyAgePos)<-rownames(CSbyAge)  ; colnames(CSbyAgePos)<-colnames(CSbyAge)
    FSbyAgePos<-  matrix(0,dim(CSbyAge)[1],dim(CSbyAge)[2]); rownames(CSbyAgePos)<-rownames(CSbyAge)  ; colnames(CSbyAgePos)<-colnames(CSbyAge)
    HSbyAgePos<-  matrix(0,dim(CSbyAge)[1],dim(CSbyAge)[2]); rownames(CSbyAgePos)<-rownames(CSbyAge)  ; colnames(CSbyAgePos)<-colnames(CSbyAge)
 }  
##-------------------------- Age Of Onset-------------------------------------##
ageAdj=-0.5             #Value to adjust Age by when integrating: to account for time step being 12 months, start/end at halfway through previous/next year when evaluating integral
maxAge=50               #Maximum age of animal, above which assume AgeOnsetDist=0   - need to define as must normalise AgeOnset for numerical integtation so need defined end point.
ageMAX=40               #how many years to calculate the numerical integration steps for - pick a number to exceed the maximum number of years you expect to forward predict - but not so large that it takes forever to calculate the integrals
TestDetectionLength=12; #number of years used in the AOTLC integral (technically goes up to infinity, but as TestSensitivity ->0 all values ~0 after a certain point -which is assumed to be this parameter, essentially anything over 2 years is negligible
### Test Sensitivity: If UserInput[1,5]=TRUE => predefined estimates   ELSE estimates are user defined
### NOTE that if user defined inputs are ill defined model may produce strange results        
ifelse(  UserInput[1,5]==TRUE, {testA<-5.94; testB<- -3.4*12},{testA<-UserInput[1,6]; testB<-UserInput[1,7]})
 #### determine age of onset values - dependent on strain type and if UK is included      
ifelse(runAllMS=="FALSE",{isUK= pmatch("UK",MSList,nomatch=0)},{ isUK= pmatch("UK",MS,nomatch=0)})

if(UserInput[1,4]=="TRUE") {  
    if(strainType=="C H L AND UNKNOWN")  {Onset<-c(2.0062,0.2816); OnsetL<-c(1.9953,0.2742); OnsetU<-c(2.017,0.2896)}
    if ((MS=="EU25"||MS=="EU17"||MS=="EU27")&& strainType=="C H L AND UNKNOWN") {Onset<-c(2.0703,0.2904); OnsetL<-c(2.0621,0.2847); OnsetU<-c(2.0785,0.2963)}
    if(isUK==1 && strainType=="C H L AND UNKNOWN")    {Onset<-c(2.1439,0.2828); OnsetL<-c(2.1322,0.2731); OnsetU<-c(2.1536,0.2896)}    
    if (strainType=="CLASSICAL AND UNKNOWN")  { Onset<-c(1.9996,0.2760); OnsetL<-c(1.9889,0.2686); OnsetU<-c(2.0103,0.2837)   } 
    if ((MS=="EU25"||MS=="EU17"||MS=="EU27")&& strainType=="CLASSICAL AND UNKNOWN")  {Onset<-c(2.0662,0.2873); OnsetL<-c(2.058,0.2816); OnsetU<-c(2.0743,0.2931)}   #My classical estimates
    if(isUK==1 && strainType=="CLASSICAL AND UNKNOWN") { Onset<-c(2.1419,0.2811); OnsetL<-c(2.1303,0.2731); OnsetU<-c(2.1536,0.2896)   }  
    if(strainType=="TYPE H"||strainType=="TYPE L"||strainType=="TYPE H AND L")  {Onset<-c(2.561,0.2433); OnsetL<-c(2.4832,0.1933); OnsetU<-c(2.6388,0.3124)}   #My classical values  
        } #end if userInput[1,4]==True
if(UserInput[1,4]=="FALSE") { Onset<-as.numeric(c(UserInput[1,32],UserInput[1,33]));OnsetL<-c(1.747,0.000122); OnsetU<-c(1.791,0.332)}    
                        
################################################################################ 
## Ammend data: Check if test positive data exists and omit periods where there is not
###check if more missing data than previously       
   TotTestCheck<-colSums(CSbyAge+ESbyAge+FSbyAge+HSbyAge )
        colSkip=0;   #some countries don't test animals until 2003 or later - if these are used need to start analysis from this point.
        colSkip0=0
        while(TotTestCheck[ colSkip0+1]<=0.1*TotTestCheck[ colSkip0+2]) { colSkip0= colSkip0+1 }  
        if(colSkip0>colSkip) colSkip=colSkip0
        if (MS=="EU8")  write.csv(colSkip,paste(wd,"/","EU8skip.csv",sep=""))
          
##### Start & End dates    #####################################################
  testStart=as.numeric(colnames(TestStartIn)[1])+ colSkip    #first testing period - accounting for if need to skip any due to lack of data
     nAgeGrp=dim(TestStartIn)[1]                             #number of age groups in model
     ageHead<-rownames(TestStartIn)
 
  testEnd=as.numeric(colnames(TestStartIn)[length(colnames(TestStartIn)[])]);   #last testing period
  nTest=testEnd-testStart+1;                                                    #number of testing periods -accounting for skipped periods
  nTest0<-dim(TestStartIn)[2]                                                   #Original number of testing periods - not accounting for skipped periods 
  maxAgeMnth=as.numeric(strsplit(rownames(TestStartIn)[nAgeGrp],"-")) ;
  maxAgeYr=maxAgeMnth/12;
  chtStart=testStart-(maxAgeYr)                                                  #date of first birth Cohort
  chtEnd=as.numeric(colnames(TestStartIn)[length(colnames(TestStartIn)[])]); 
  nCht=chtEnd-chtStart+1                                                         #number of birth cohorts in model
                  
 ### remove columns with no test data.
  if (colSkip>0) {
        TestStartIn <- TestStartIn[,(1+colSkip):nTest0]
        CSbyAge <-CSbyAge[,(1+colSkip):nTest0]
        ESbyAge <-ESbyAge[,(1+colSkip):nTest0]
        FSbyAge <-FSbyAge[,(1+colSkip):nTest0]
        HSbyAge <-HSbyAge[,(1+colSkip):nTest0]
        
        CSbyAgePos <-CSbyAgePos[,(1+colSkip):nTest0]
        ESbyAgePos <-ESbyAgePos[,(1+colSkip):nTest0]
        FSbyAgePos <-FSbyAgePos[,(1+colSkip):nTest0]
        HSbyAgePos <-HSbyAgePos[,(1+colSkip):nTest0]
        
        ESNumDeadbyAge<-ESNumDeadbyAge[,(1+colSkip):nTest0]
        FSNumDeadbyAge<-FSNumDeadbyAge[,(1+colSkip):nTest0]
        HSNumDeadbyAge<-HSNumDeadbyAge[,(1+colSkip):nTest0]
  }
################################### PARAMETERS #################################

rAtc= dim(AgeToCht)[1]; cAtc=dim(AgeToCht)[2]
### If AgeToCht matrix too small then add in extra rows: ASSUMES SAME SPLIT AS OTHER ADULT CATTLE
if (rAtc > nAgeGrp) {AgeToCht<-AgeToCht[1:nAgeGrp,]}
if (rAtc < nAgeGrp) {
   AgeToCht0<-AgeToCht
   nAtc=nAgeGrp-rAtc
   AgeToCht<-matrix(0,rAtc+nAtc,cAtc+nAtc)
   AgeToCht[1:rAtc,1:cAtc]=AgeToCht0
   for (i in 1:nAtc) { AgeToCht[(i+1):(i+rAtc),cAtc+i]=AgeToCht0[,cAtc]}
   }
   
  RiskRate<-matrix(0,nCht,nTest)        #set up matrix for Risk for CS & FS
  ageQ=matrix(1:(testEnd-chtStart+1))
  minAgeTest=0                          # animals below this age will not be tested and thus should not have a risk of positive detection - set to 0 as using all data for baseline model
         
  RiskRateHS<-matrix(0,nCht,nTest)      #set up matrix for risk for HS & FS
  ageQ_HS=matrix(1:(testEnd-chtStart+1))
  ##EU25 proportion test positive
  pCSTestPosEU25<-  read.csv(paste(wd,"/","EU25pCSTestPos.csv",sep=""),row.names=1) [(28-nCht+1):28,(1+colSkip):nTest0]
  pHSTestPosEU25<-  read.csv(paste(wd,"/","EU25pHSTestPos.csv",sep=""),row.names=1) [(28-nCht+1):28,(1+colSkip):nTest0]
  ## PCS and PHS alternative estimates for Emergence by Age using EU25 Test positive data
  ## 3 estimates: col 1= average of all years ; col2=average of first 3 years ; col3 = average of last 4 years
  ## where there are no test positives values for a given age group, the value from the average of all years or a close age group is used
  pCSEstEU25<-  read.csv(paste(wd,"/","EU25pCSEsts.csv",sep=""),row.names=1)
  pHSEstEU25<-read.csv(paste(wd,"/","EU25pHSEsts.csv",sep=""),row.names=1)
################################################################################
#### Parameters for GetMLEs  
################################################################################
  DoLowerCI=0;          #set=1 if you want to estimate lower confidence intervals in the GetMLE's function
  Specificity = 1.0;    #specificity of test
  MaxSensitivity = 0.99;# sensitivity of test for CS & FS
  
################################################################################ 
######### Design Prevalence - Parameters
################################################################################
  dP=1/UserInput[1,14]   # define what the required design prevalence is: i.e. the prevalence at you want to be able to detect
  tau=UserInput[1,15]    # Power - fixed value
  upDpLim=1e-2;        #upper limit for solver routine when finding Dp in the fixed testing scenario.
  dpSolve0=1e-5        #inital value for solver routine when finding Dp in the fixed testing scenario.

################################################################################
#### EMERGENCE  - Parameters
################################################################################
 fP=100                      # max number of years to forward predict - extra num years=extra num new birth cohorts
 incRate=UserInput[1,28]     # rate at which trend increases per year - same rate each year - if change rate need to be a bit more clever
 if ( UserInput[1,29]==TRUE && is.numeric(UserInput[1,30])) incRate<-UserInput[1,30]    # use user defined increasing trend - only if proper value used
 ifelse(UserInput[1,52]==TRUE,trendMeth<-1,trendMeth<-2)  # define which method to use for emergence: 1 =95% CI method, 2=number of cases method
#Get rows corresponding to ages of animals you want to test for detection of emergence
emerAgeStart=gsub(">","",UserInput[1,54])    #start Age to test emergence
ifelse(UserInput[1,55]=="NO LIMIT",emerAgeEnd<-ageHead[length(ageHead)], emerAgeEnd<-gsub("<","",UserInput[1,55]) )      #end age to test emergence  (if no limit set equal to max row in ageHead)
 AgeDetEmer2<-pmatch(emerAgeStart,ageHead):pmatch(emerAgeEnd,ageHead)    #get corresponding rows
 
 trendInc=matrix(0,1,fP+1)   # set up matrix : baseline
 trendIncTrue=matrix(0,1,fP+1)   # set up matrix  -infected number missed counter : baseline

 trendIncS=matrix(0,1,fP+1)   # set up matrix : scenario
 trendIncTrueS=matrix(0,1,fP+1)   # set up matrix  -infected number missed counter: scenario
################################################################################
#### Output folders
################################################################################ 
  
today<-Sys.Date()
today<-format(today,format="%d.%m.%y")
mainOutDir<-file.path(paste("C:\\C-TSEMM\\Output Data Files\\",today,"(tau=",tau,")",sep=""))
outSuff="_v1"
  ifelse (file.exists(mainOutDir),setwd(mainOutDir), {
          dir.create(file.path(mainOutDir))
          setwd(mainOutDir) }  )      #create output file if doesn't exist
  subDir<-paste(MSName,sep="")
  iDir=1
  while (file.exists(subDir)){ 
      iDir=iDir+1;
      subDir <-  paste(MSName,"_v",iDir,sep="")
       outSuff <-  paste("_v",iDir,sep="")    #suffix to add to end of output files
      }     
  dir.create(file.path(mainOutDir,subDir))
  outDir<-file.path(mainOutDir,subDir)
  dir.create( file.path(mainOutDir,subDir,"Validation"))
  outDirPrev<-file.path(mainOutDir,subDir,"Validation")
  dir.create( file.path(mainOutDir,subDir,"Emer2013"))
  dir.create( file.path(mainOutDir,subDir,"Emer2013/base"))
  dir.create( file.path(mainOutDir,subDir,"Emer2013/scen"))
 setwd(wd)

 ##############################################################################
 ################################################################################################
