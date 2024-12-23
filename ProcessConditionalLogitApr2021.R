# 
# This file contains R functions to process all conditions within a given dataset
# it takes in 1) a dataframe of the form:
#     associations1_10<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double())
#     2) a dataset to process for example: medications_file_ageproper.conditionsSens1_10
#     3) Which column to begin processing e.g. (36)


# updated 2024 to adjust for BMI/SES and Medications in reviewer requested sensitivity analyses


ProcessConditionsMedications<-function(DF, DS, startcolumn, errorDS,warningDS,verbose=F){
  
  for(i in seq(startcolumn, length(DS), 1))
  {
    if(verbose){print(paste("i = ", i))}
    CS<-unlist(DS$CASE_STATUS)
    if(verbose){print(table(CS))}
    condition<-unlist(DS[,i])
    if(verbose){print(table(condition))}
    SETID<-unlist(DS$SETID)
    indexage<-unlist(DS$indexage)
    crosstabs<-table(CS,condition)
    if(verbose){print(crosstabs)}
    
    if(dim(crosstabs)[2]>1){
      exposedcase<-crosstabs[2,2]
      exposedcontrol<-crosstabs[1,2]
      if(exposedcase > 0 & exposedcontrol > 0){
        tryCatch({
          clr<-clogit(CS~condition+indexage+strata(SETID),coxph.control(iter.max = 5000))
          if(verbose){print(paste("conditional logisitic regression  converged", i, sep=" "))}
          condition<-unlist(colnames(DS))[i]
          if(verbose){print("extracting OR")}
          OR<-summary(clr)$coef[1,2]
          ORlow<-exp(summary(clr)$coef[1,1]-1.96*summary(clr)$coef[1,3])
          ORhigh<-exp(summary(clr)$coef[1,1]+1.96*summary(clr)$coef[1,3])
          pvalue<-summary(clr)$coef[1,5]
          
          if(verbose){print(paste(condition,format(round(OR,3),nsmall=3),format(round(ORlow,3),nsmall=3),format(round(ORhigh,3),nsmall=3),format(round(pvalue,7),nsmall=7),sep = " "))}
          DF<-rbind(DF,cbind(condition,exposedcase,exposedcontrol,OR=format(round(OR,3),nsmall=3),ORlow=format(round(ORlow,3),nsmall=3),ORhigh=format(round(ORhigh,3),nsmall=3),pval=format(round(pvalue,7),nsmall=7)))
     
          
          
             },
        error=function(e){
          condition<-unlist(colnames(DS))[i]
          errorstring<-cbind(condition,message(e))
          print(paste("The following error message appeared in clogit for condition: ",condition,message(e),sep= " "))
          message(e)
          errorsDS<-rbind(errorstring,condition)
        },
        warning=function(e){
          condition<-unlist(colnames(DS))[i]
          print(paste("The following warning message appeared in clogit for condition: ",condition,message(e),sep= " "))
          message(e)
          condition<-unlist(colnames(DS))[i]
          OR<-NA
          ORlow<-NA
          ORhigh<-NA
          pvalue<-NA
          if(verbose){print(paste(condition,exposedcase,exposedcontrol,OR=format(round(OR,3),nsmall=3),ORlow=format(round(ORlow,3),nsmall=3),ORhigh=format(round(ORhigh,3),nsmall=3),pval=format(round(pvalue,7),nsmall=7),sep = " "))}
          DF<-rbind(DF,cbind(condition,exposedcase,exposedcontrol,OR=format(round(OR,3),nsmall=3),ORlow=format(round(ORlow,3),nsmall=3),ORhigh=format(round(ORhigh,3),nsmall=3),pval=format(round(pvalue,7),nsmall=7)))
        })
      }else{
        print(paste0(colnames(DS)[i],": The number of observed conditions in cases = 0; Not estimatable"))
        condition<-unlist(colnames(DS))[i]
        OR<-NA
        ORlow<-NA
        ORhigh<-NA
        pvalue<-NA
        if(verbose){print(paste(condition,exposedcase,exposedcontrol,format(round(OR,3),nsmall=3),format(round(ORlow,3),nsmall=3),format(round(ORhigh,3),nsmall=3),format(round(pvalue,7),nsmall=7),sep = " "))}
        DF<-rbind(DF,cbind(condition,exposedcase,exposedcontrol,OR=format(round(OR,3),nsmall=3),ORlow=format(round(ORlow,3),nsmall=3),ORhigh=format(round(ORhigh,3),nsmall=3),pval=format(round(pvalue,7),nsmall=7)))
      }
    }else{
      print(paste0(colnames(DS)[i],": The number of observed conditions in cases and controls = 0; Not estimatable"))
      exposedcase<-0
      exposedcontrol<-0
      condition<-unlist(colnames(DS))[i]
      OR<-NA
      ORlow<-NA
      ORhigh<-NA
      pvalue<-NA
      if(verbose){print(paste(condition,exposedcase,exposedcontrol,format(round(OR,3),nsmall=3),format(round(ORlow,3),nsmall=3),format(round(ORhigh,3),nsmall=3),format(round(pvalue,7),nsmall=7),sep = " "))}
      DF<-rbind(DF,cbind(condition,exposedcase,exposedcontrol,OR=format(round(OR,3),nsmall=3),ORlow=format(round(ORlow,3),nsmall=3),ORhigh=format(round(ORhigh,3),nsmall=3),pval=format(round(pvalue,7),nsmall=7)))
    }
  }
  return(DF)
}

# All following functions use similar structure as function above, with different adjustments in the model

# sensitivity analyses for reviewers
# function that performs conditional logistic regression additionally mutually adjusted for both BMI (factor) and 
#    SES (deciles of multiple deprivation index; continuous linear term 1-10)
# Function takes in an empty data frame (DF), dataset (DS), start column containing condition (startcolumn)
# returns a full data frame (DF) with condition, unadjusted associations and BMI+SES adjusted associations
ProcessConditions_BMISESAdjusted<-function(DF, DS, startcolumn, errorDS, warningsDS,verbose=F){
  for(i in seq(startcolumn, length(DS), 1))
  {
    if(verbose){print(paste("i = ", i))}
    CS<-unlist(DS$CASE_STATUS)
    condition<-unlist(DS[,i])
    SETID<-unlist(DS$SETID)
    indexage<-unlist(DS$indexage)
    IMD<-unlist(DS$IMD)
    BMI<-unlist(DS$BMI_M12_ALL.f)
    crosstabs<-table(CS,condition)
    if(verbose){print(crosstabs)}
    exposedcase<-0
    exposedcontrol<-0
    OR<-0.0
    ORlow<-0.0
    ORhigh<-0.0
    p<-0.0
    
    OR1<-0.0
    OR1low<-0.0
    OR1high<-0.0
    p1<-0.0    
    if(dim(crosstabs)[2]>1){
      exposedcase<-crosstabs[2,2]
      exposedcontrol<-crosstabs[1,2]
      if(exposedcase > 0 & exposedcontrol > 0){
        tryCatch({

          clr<-clogit(CS~condition+indexage+strata(SETID),coxph.control(iter.max = 5000))

          if(verbose){print(paste("age adjusted conditional logisitic regression  converged", i, sep=" "))}
          if(verbose){print("extracting OR")}
          OR<-summary(clr)$coef[1,2]
          if(verbose){print(summary(clr))}
          ORlow<-exp(summary(clr)$coef[1,1]-1.96*summary(clr)$coef[1,3])
          ORhigh<-exp(summary(clr)$coef[1,1]+1.96*summary(clr)$coef[1,3])
          p<-summary(clr)$coef[1,5]
            
          clr1<-clogit(CS~condition+indexage+IMD+BMI+strata(SETID),coxph.control(iter.max = 5000))
          if(verbose){print(paste("age, BMI, SES conditional logisitic regression  converged", i, sep=" "))}

          if(verbose){print("extracting adjusted OR")}
          OR1<-summary(clr1)$coef[1,2]
          OR1low<-exp(summary(clr1)$coef[1,1]-1.96*summary(clr1)$coef[1,3])
          OR1high<-exp(summary(clr1)$coef[1,1]+1.96*summary(clr1)$coef[1,3])
          p1<-summary(clr1)$coef[1,5]
          if(verbose){print(summary(clr1))}
          
          condition<-unlist(colnames(DS))[i]        
          if(verbose){print(paste(condition,OR=format(round(OR,3),nsmall=3),ORlow=format(round(ORlow,3),nsmall=3),
                      ORhigh=format(round(ORhigh,3),nsmall=3),p=format(round(p,10),nsmall=10),
                      OR1=format(round(OR1,3),nsmall=3),OR1low=format(round(OR1low,3),nsmall=3),
                      O1high=format(round(OR1high,3),nsmall=3),p1=format(round(p1,10),nsmall=10),
                      sep = " "))}
          DF<-rbind(DF,cbind(condition,exposedcase,exposedcontrol,
                             OR=format(round(OR,3),nsmall=3),
                             ORlow=format(round(ORlow,3),nsmall=3),
                             ORhigh=format(round(ORhigh,3),nsmall=3),
                             p=format(round(p,10),nsmall=10),
                             OR1=format(round(OR1,3),nsmall=3),
                             OR1low=format(round(OR1low,3),nsmall=3),
                             OR1high=format(round(OR1high,3),nsmall=3),
                             p1=format(round(p1,10),nsmall=10)))
        },
        error=function(e){
          condition<-unlist(colnames(DS))[i]
          errorstring<-cbind(condition,message(e))
          print(paste("The following error message appeared in clogit for condition: ",condition,message(e),sep= " "))
          message(e)
          errorsDS<-rbind(errorstring,condition)
        },
        warning=function(e){
          condition<-unlist(colnames(DS))[i]
          print(paste("The following warning message appeared in clogit for condition: ",condition,message(e),sep= " "))
          message(e)
          condition<-unlist(colnames(DS))[i]
          if(verbose){print(paste(condition,exposedcase,exposedcontrol,
                      OR=format(round(OR,3),nsmall=3),
                      ORlow=format(round(ORlow,3),nsmall=3),
                      ORhigh=format(round(ORhigh,3),nsmall=3),
                      p=format(round(p,10),nsmall=10),
                      OR1=format(round(OR1,3),nsmall=3),
                      OR1low=format(round(OR1low,3),nsmall=3),
                      OR1high=format(round(OR1high,3),nsmall=3),
                      p1=format(round(p1,10),nsmall=10),sep = " "))}      
          
          DF<-rbind(DF,cbind(condition,exposedcase,exposedcontrol,
                             OR=format(round(OR,3),nsmall=3),
                             ORlow=format(round(ORlow,3),nsmall=3),
                             ORhigh=format(round(ORhigh,3),nsmall=3),
                             p=format(round(p,10),nsmall=10),
                             OR1=format(round(OR1,3),nsmall=3),
                             OR1low=format(round(OR1low,3),nsmall=3),
                             OR1high=format(round(OR1high,3),nsmall=3),
                             p1=format(round(p1,10),nsmall=10)))
        })
      }else{
        print(paste0(colnames(DS)[i],": The number of observed conditions in cases = 0; Not estimatable"))
        condition<-unlist(colnames(DS))[i]
        if(verbose){print(paste(condition,exposedcase,exposedcontrol,
                    OR=format(round(OR,3),nsmall=3),
                    ORlow=format(round(ORlow,3),nsmall=3),
                    ORhigh=format(round(ORhigh,3),nsmall=3),
                    p=format(round(p,10),nsmall=10),
                    OR1=format(round(OR1,3),nsmall=3),
                    OR1low=format(round(OR1low,3),nsmall=3),
                    OR1high=format(round(OR1high,3),nsmall=3),
                    p1=format(round(p1,10),nsmall=10),sep = " "))}
        
        DF<-rbind(DF,cbind(condition,exposedcase,exposedcontrol,
                 OR=format(round(OR,3),nsmall=3),
                 ORlow=format(round(ORlow,3),nsmall=3),
                 ORhigh=format(round(ORhigh,3),nsmall=3),
                 p=format(round(p,10),nsmall=10),
                 OR1=format(round(OR1,3),nsmall=3),
                 OR1low=format(round(OR1low,3),nsmall=3),
                 OR1high=format(round(OR1high,3),nsmall=3),
                 p1=format(round(p1,10),nsmall=10)))
        }
    }else{
      condition<-unlist(colnames(DS))[i]
      print(paste0(colnames(DS)[i],": The number of observed conditions in cases and controls = 0; Not estimatable"))
      if(verbose){print(paste(condition,exposedcase,exposedcontrol,
                  OR=format(round(OR,3),nsmall=3),
                  ORlow=format(round(ORlow,3),nsmall=3),
                  ORhigh=format(round(ORhigh,3),nsmall=3),
                  p=format(round(p,10),nsmall=10),
                  OR1=format(round(OR1,3),nsmall=3),
                  OR1low=format(round(OR1low,3),nsmall=3),
                  OR1high=format(round(OR1high,3),nsmall=3),
                  p1=format(round(p1,10),nsmall=10),sep = " "))}
      
      DF<-rbind(DF,cbind(condition,exposedcase,exposedcontrol,
                         OR=format(round(OR,3),nsmall=3),
                         ORlow=format(round(ORlow,3),nsmall=3),
                         ORhigh=format(round(ORhigh,3),nsmall=3),
                         p=format(round(p,10),nsmall=10),
                         OR1=format(round(OR1,3),nsmall=3),
                         OR1low=format(round(OR1low,3),nsmall=3),
                         OR1high=format(round(OR1high,3),nsmall=3),
                         p1=format(round(p1,10),nsmall=10)))   }
  }
  return(DF)
}

###### May 2024 for revision.
# Function to adjust for medication use that might be associated with LCINS (mediate the condition-LCINS association)
## examples, oral corticosteroids, immunosuppressing medications, PPIs, H2Rec
EstimateAdjAssociation<-function(DF,Case,indexage,SETID,Condition,Medication,Labels,verbose=F){
  if(verbose){print(dim(DF))}
  table(DF[[Condition]])
  if(verbose){print({Labels}[{{Condition}}])}
  ConditionReadable<-{Labels}[{{Condition}}]
  if(verbose){print({Labels}[{{Medication}}])}
  MedicationReadable<-{Labels}[{{Medication}}]
  clr_age<-clogit(DF[[Case]]~DF[[Condition]]+DF[[indexage]]+strata(DF[[SETID]]),
                  coxph.control(iter.max = 5000))
  
  OR<-summary(clr_age)$coef[1,2]
  ORlow<-exp(summary(clr_age)$coef[1,1]-1.96*summary(clr_age)$coef[1,3])
  ORhigh<-exp(summary(clr_age)$coef[1,1]+1.96*summary(clr_age)$coef[1,3])
  pvalue<-summary(clr_age)$coef[1,5]
  
  clr_med<-clogit(DF[[Case]]~DF[[Condition]]+DF[[Medication]]+DF[[indexage]]+strata(DF[[SETID]]),
                  coxph.control(iter.max = 5000))
  
  OR1<-summary(clr_med)$coef[1,2]
  OR1low<-exp(summary(clr_med)$coef[1,1]-1.96*summary(clr_med)$coef[1,3])
  OR1high<-exp(summary(clr_med)$coef[1,1]+1.96*summary(clr_med)$coef[1,3])
  pval1<-summary(clr_med)$coef[1,5]
  
  OR2<-summary(clr_med)$coef[2,2]
  OR2low<-exp(summary(clr_med)$coef[2,1]-1.96*summary(clr_med)$coef[2,3])
  OR2high<-exp(summary(clr_med)$coef[2,1]+1.96*summary(clr_med)$coef[2,3])
  pval2<-summary(clr_med)$coef[2,5]
  if(verbose){print(paste({{Condition}},{{Medication}},ConditionReadable, MedicationReadable, format(round(OR,3),nsmall=3),format(round(ORlow,3),nsmall=3),format(round(ORhigh,3),nsmall=3),
              format(round(pvalue,7),nsmall=7),
              paste(format(round(OR,2),nsmall=2)," (",format(round(ORlow,2),nsmall=2),",",format(round(ORhigh,2),nsmall=2),")",sep=""),
              format(round(OR1,3),nsmall=3),format(round(OR1low,3),nsmall=3),format(round(OR1high,3),nsmall=3),
              format(round(pval1,7),nsmall=10),
              paste(format(round(OR1,2),nsmall=2)," (",format(round(OR1low,2),nsmall=2),",",format(round(OR1high,2),nsmall=2),")",sep=""),
              format(round(OR2,3),nsmall=3),format(round(OR2low,3),nsmall=3),format(round(OR2high,3),nsmall=3),
              paste(format(round(OR2,2),nsmall=2)," (",format(round(OR2low,2),nsmall=2),",",format(round(OR2high,2),nsmall=2),")",sep=""),
              format(round(pval2,7),nsmall=7),sep = " "))}
  
  #return associations of unadjusted, medication-adjusted
  #add in formatted OR, 95%CI
  associations<-data.frame(cbind(Condition,Medication,ConditionReadable, MedicationReadable,
                                 format(round(OR,3),nsmall=3),format(round(ORlow,3),nsmall=3),format(round(ORhigh,3),nsmall=3),
                           paste(format(round(OR,2),nsmall=2)," (",format(round(ORlow,2),nsmall=2),",",format(round(ORhigh,2),nsmall=2),")",sep=""),
                                 format(round(pvalue,7),nsmall=10),
                                 format(round(OR1,3),nsmall=3),
                                 format(round(OR1low,3),nsmall=3),
                                 format(round(OR1high,3),nsmall=3),
                          paste(format(round(OR1,2),nsmall=2)," (",format(round(OR1low,2),nsmall=2),",",format(round(OR1high,2),nsmall=2),")",sep=""),
                                 format(round(pval1,10),nsmall=10),
                                format(round(OR2,3),nsmall=3),
                                 format(round(OR2low,3),nsmall=3),
                                 format(round(OR2high,3),nsmall=3),
                          paste(format(round(OR2,2),nsmall=2)," (",format(round(OR2low,2),nsmall=2),",",format(round(OR2high,2),nsmall=2),")",sep=""),
                                 format(round(pval2,10),nsmall=10)))
  return (associations)
}
