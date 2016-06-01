################################################################################
####### EFSA BSE Surveillance model: GetMLEs Function             #############
################################################################################
GetMLEs<-function( paras,aa,bbU,bbL,Wbl) {                             
  if(estTrend==1) {      #use MLE method to to find best values for A,B&K, if estTrend=1, else use pre-defined values of paras (still need to get Logliklihood )           
             prevsOptim<-nlminb(paras,SumPrevLogK,gradient=NULL,hessian=NULL,aa,PC,PD,PS,Wbl)          
  ifelse(useOnsetCI==TRUE,{                    
            #### CI: Using Age of Onset CI values
             prevsOptimUci<-optim(paras,SumPrevLogK,gr=NULL,aa,PCU,PDU,PSU,Wbl)
            prevsOptimLci<-optim(paras,SumPrevLogK,gr=NULL,aa,PCL,PDL,PSL,Wbl) 
            },{
            #### CI: Poisson about input values                                               
            prevsOptimUci<-optim(paras,SumPrevLogK,gr=NULL,bbU,PC,PD,PS,Wbl)   #### Save output
            prevsOptimLci<-optim(paras,SumPrevLogK,gr=NULL,bbL,PC,PD,PS,Wbl)
            })
  parasB<-prevsOptim$par ; fvalsB<-prevsOptim$objective 
  parasU<-prevsOptimUci$par ; fvalsUci<-prevsOptimUci$value 
  parasL<-prevsOptimLci$par ; fvalsLci<-prevsOptimLci$value  
  }
  if(estTrend==2) {     #fixed input values
            parasIn<- as.matrix( read.csv(paste(wd,"/",EUNum,"paras.csv",sep=""),row.names=1)  )
            parasB<-as.numeric(parasIn[1,]) ; fvalsB<-NA 
            parasU<-as.numeric(parasIn[2,])  ; fvalsUci<-NA 
            parasL<-as.numeric(parasIn[3,])  ; fvalsLci<-NA          
                  }
#----------------- OUTPUT ------------------------------------------------------
list(paras=parasB,parasL=parasL,parasU=parasU,fvalsB=fvalsB,fvalsUci=fvalsUci,fvalsLci=fvalsLci)
}  #end GetMLEs
