################################################################################
#### EFSA BSE Surveillance model: BSE_Functions.r  ############################
#### File containing Model Functions   #####################################
################################################################################

#####Estimate probability of onset of clinical symptoms at age a
AgeOnsetDist<-function(animalAge,Poly,knot,Onset){
       ageVal<-dlnorm(animalAge,Onset[1],Onset[2])
       ifelse( animalAge>maxAge,  PolyTmp1<-0, PolyTmp1<-max(0,ageVal))     #if Animal>maxAge assume AgeOnset=0.
}
#-------------------------------------------------------------------------------
####AgeOnsetDist including test sensitivity - i.e. prob. of detection if tested at age 'offest' given clinical onset at age t
AOTLC<-function(t,offset,Poly,knot,Onset) {
  TestSensitivity<-exp(testA+testB*(t-offset))/(1.0+exp(testA+testB*(t-offset)))
  ageVal<-dlnorm(t,Onset[1],Onset[2])
  ifelse( t>maxAge,  AgeOnsetD<-0, AgeOnsetD<-max(0,ageVal))
  ifelse(t<0,0,AgeOnsetD* TestSensitivity )
}
################################################################################ 
#### convert Age data to Cohort Data
################################################################################
AgeConv<-function(byAge)  {
  byCht=matrix(0,nCht,nTest)
 for (j in 1:nTest)    {
    for (i in 1:nCht)  {   
        byCht[i,j]=  sum(AgeToCht[,ageMatNA[i,j]+1]*(byAge[,j]))
            }  # end i
            }  # end j
byCht[is.na(byCht)]=0 #set NA values =0 - is NA if testing period earlier than year of birth
byCht   }
################################################################################ 
#### Emergence method 2: additional project update
################################################################################     
chtConv1203<-function(byCht)  {                       #defining the cohort to age matrices just for this section
      byAge=matrix(0,nAgeGrp2,1)                      #set up age matrix
      n1=size(byCht)[2]                               #define start number                   
      byCht2=byCht[n1:max(1,(n1-nAgeGrp2))]           #reverse order of byCht but only go back maximum of nAgeGrp2 years
      for (i in 1:nAgeGrp2)  {         
        byAge[i,]=  sum(byCht2*ChtToAge3[i,1:length(byCht2)])     #multiply sum of byCht by chtToAge matrix to convert cohort to age.
            }  # end i
       byAge  }  # end function

 AgeConv1203<-function(byAge)  {                      #defining the age to cohort matrices just for this section
  byCht=matrix(0,nCht,nTest)
 for (j in 1:nTest)    {
    for (i in 1:nCht)  {   
        byCht[i,j]=  sum(AgeToCht2[1:nrow(byAge),ageMatNA[i,j]+1]*(byAge[,j]))
            }  # end i
            }  # end j
  byCht[is.na(byCht)]=0                               #set NA values =0 - is NA if testing period earlier than year of birth
  byCht   }


################################################################################
##Function to plot colSums of model output and upper, lower CIS's 
################################################################################
plotCI<-function(avg,lci,uci,obs,figName)       {
  #pdf(file=paste(wdOut,sep="","/",figName))#, height=3.5, width=5)       #if you want output as a pdf file
     png(file=file.path(outDirPrev,figName),width=580,height=580  )       #if you want output as a png file.
  colour<-c("magenta","seagreen","royalblue","tomato")
  g_range <- range(avg,lci,uci)
  plot(obs,ylim=g_range,ylab="Number Positive",xlab="Testing Period",axes=F)
  axis(1, lab=F)
  text(1:length(colnames(TestStartIn)), par("usr")[3]-(par("usr")[4]-par("usr")[3])/100, srt=45, adj=1,
          labels=colnames(TestStartIn),
          xpd=T, cex=0.8)
  axis(2, las=1, cex.axis=0.8)
  box()

  lines(obs,col=colour[1],lwd=2)   
  lines(avg,col=colour[2],lwd=2)
  lines(lci,col=colour[3],lwd=2)
  lines(uci,col=colour[4],lwd=2)

  legend(7, g_range[2], c("Observed","baseline","lower CI","upper CI"), cex=1, col=colour,lwd=2, pch=21, lty=1);
  dev.off()
}
################################################################################
##Function to plot baseline vs Scenario test positives
################################################################################
plotNumMissed<-function(baseline,scenario,figName)       {
  tiff(file=file.path(outDirPrev,figName),width=480,height=480  )
  colour<-c("seagreen","royalblue")
  g_range <- range(baseline,scenario)
  plot(baseline,ylim=g_range,ylab="Number Positive",xlab="Testing Period",axes=F)
  axis(1, lab=F)
  text(1:length(colnames(TestStartIn)), par("usr")[3] - 2, srt=45, adj=1,
          labels=colnames(TestStartIn),
          xpd=T, cex=0.8)
  axis(2, las=1, cex.axis=0.8)
  box()

  lines(baseline,col=colour[1],lwd=2)
  lines(scenario,col=colour[2],lwd=2)

  legend(7, g_range[2], c("baseline","scenario"), cex=1, col=colour,lwd=2, pch=21, lty=1);
  dev.off()
}
################################################################################
## ammend byAge inputs for Scenario - no testing on animals below ageCap and pctDecrease in animals tested above ageCap
################################################################################
ScenarioIn<-function(stream,ageCap,ageCapEnd,prTest) {
  streamDim<-dim(stream)
  if(is.na(ageCapEnd)) ageCapEnd=streamDim[1]
  streamS<-stream*0
  streamS[(ageCap+1):ageCapEnd,]=stream[(ageCap+1):ageCapEnd,]*(prTest)
  streamS 
  }
############################################################################
##### log likelihood Function for Get MLEs
################################################################
SumPrevLogK<-function( paras,aa,PC,PD,PS,Wbl )  {
    ifelse(Wbl==1, { A=paras[1]; B=paras[2];C=paras[3]      #Estimated Weibull distribution given current parameter values
      Prev=C* ((1:nCht) / A)^(B-1)* exp(-((1:nCht) / A)^B); 
      K=  exp(paras[4])/(1+exp(paras[4])) },                #Estimated value of differential slaughter, K, given current parameter vlaue; paras[3]     
  {    
      Prev=paras[1]*exp(paras[2]*((1:nCht)-1));             #Estimated Exponential distribution given current parameter values; paras[1] & paras[2]
      K=  exp(paras[3])/(1+exp(paras[3]))                   #Estimated value of differential slaughter, K, given current parameter vlaue; paras[3]     
  })
      Prev2=t(repmat(Prev,nTest,1))                         #turn Prev into a nCht*nTest matrix - so trend value for each [cht,testPeriod] pair
      pclin=(1-K)*Prev2*PC*MaxSensitivity                   #Estimate prob of test positive for CS and FS
          pclin[pclin<=0]=10^-200                              ### p_CLINICAL <=0 possible due to rounding errors
           pclin[pclin>=0.99999]=0.99999                       ### p_CLINICAL >=1 possible due to rounding errors
      psurv = K * Prev2 * ( 1 - PS) * PD;                   #Estimate prob of test positive for HS and ES
          psurv[psurv<=0]=10^-200                              ### p_CLINICAL <=0 possible due to rounding errors
           psurv[psurv>=0.99999]=0.99999                       ### p_CLINICAL >=1 possible due to rounding errors 
      LcMat= -( (aa$CSFSNumTest -aa$CSFSTestPos)*log(1-pclin) + aa$CSFSTestPos * log(pclin) )           #Log likelihood for CS and FS
      LsMat= -( (aa$HSESNumTest - aa$HSESTestPos       )* log(1-psurv)+ aa$HSESTestPos * log (psurv) )  #Log likelihood for HS and ES
      sum((LcMat+LsMat),na.rm=TRUE)                          #Sum of Log Liklihood - need to minimise this.
}  

###############################################################################
#################### FUNCTION TEST SENSITIVITY  ##################################
TestSensitivity<-function(t,offset) {
  exp(testA+testB*(t-offset))/(1.0+exp(testA+testB*(t-offset)))
}

