# Main Aurum (validation data) analyses. 
# run conditional logistic regression analyses and print out to text file

#set your working directory
WD<-"./"
setwd(WD)

#read in data file
load( file = paste0(WD,"AurumDeDupped.May2024.RData"))

###################################################  ###################################################
#perform conditional logistic regression analyses with non-duplicated population
#(first 1-10 years; then 10-32 years)

NonDupAurumFullMatches$CASE_STATUS<-NonDupAurumFullMatches$case_status
lt10_idx<-which(grepl( "_lt10" , colnames( NonDupAurumFullMatches ) ) )

#rearrange the columns to put created demographic variables (lines 298:321 are created demographic variables)
NonDupAurumFullMatchesSens1_10<-NonDupAurumFullMatches[,c(1:25,298:321,lt10_idx)]

#now remove date columns
d_idx<-which(grepl("_date",colnames( NonDupAurumFullMatchesSens1_10)))
NonDupAurumFullMatchesSens1_10<-NonDupAurumFullMatchesSens1_10[,-d_idx]

#remove some medication columns and other columns that are not needed for these analyses
# remove columns 70-75, 89-116
NonDupAurumFullMatchesSens1_10<-NonDupAurumFullMatchesSens1_10[,-c(70:75,89:116)]

####overall sensitive conditions (1+ diagnosis code)
#### contains conditions and medications
errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",110)
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double())
DS<-NonDupAurumFullMatchesSens1_10
DF<-ProcessConditionsMedications(DF,DS,49,errorsDS,warningDS)

#this file is the basis of table 4. 
write.table(DF,file=paste0(WD,"NonDupAurumFullMatchesSens1_10.Apr182024.txt"),row.names = F,sep="\t")

##### #####################now stratify by sex
###############  sensitive conditions for females (GENDER==2))
DS<-NonDupAurumFullMatchesSens1_10
errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",13)
female_idx<-which(DS$GENDER==2)

#extract out specific columns to validate. The dataset includes medications and "specific" condition definitions
validating_idx<-which(colnames(DS)%in%c("tb_lt10","flu_lt10","uresp_lt10","gi_gastrexact_lt10","gi_gerd_lt10",
                                        "ai_dmt1_lt10","ai_lupus_lt10","ii_ai_psor_lt10","ai_oth2_lt10","ai_rose_exact_lt10",
                                        "copd_lt10","anemiaclose_lt10","gi_ulcerall_lt10"))
ConditionLabel=rbind("TB","Influenza","Upper Respiratory Infections","Gastritis/NIGC","GERD",
                     "DMT1","Lupus","Psoriasis", "Other AI NOS", "Rosascea","COPD","Anemia","GI Ulcer")
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double())
DSFemale<-DS[female_idx,c(1:48,validating_idx)]
DF<-ProcessConditionsMedications(DF,DSFemale,49,errorsDS,warningDS)
DF<-cbind(ConditionLabel,DF)
write.table(DF,file=paste0(WD,"NonDupAurumFullMatchesSens1_10.Female.Oct2024.txt"),row.names = F,sep="\t")

############### sensitive conditions for males (GENDER==1))
DS<-NonDupAurumFullMatchesSens1_10
errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",13)
male_idx<-which(DS$GENDER==1)

#extract out specific columns to validate. The dataset includes medications and "specific" condition definitions
validating_idx<-which(colnames(DS)%in%c("tb_lt10","flu_lt10","uresp_lt10","gi_gastrexact_lt10","gi_gerd_lt10",
                                        "ai_dmt1_lt10","ai_lupus_lt10","ii_ai_psor_lt10","ai_oth2_lt10","ai_rose_exact_lt10",
                                        "copd_lt10","anemiaclose_lt10","gi_ulcerall_lt10"))
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double())
DSMale<-DS[male_idx,c(1:48,validating_idx)]
DF<-ProcessConditionsMedications(DF,DSMale,49,errorsDS,warningDS)
DF<-cbind(ConditionLabel,DF)

write.table(DF,file=paste0(WD,"NonDupAurumFullMatchesSens1_10.Male.Oct2024.txt"),row.names = F,sep="\t")

####################################################################
# now do 10-32 years
###########################################################################################################

# now do conditions 10-32 years (sensitive definition)
# first identify people with at least 10 years of registration
NonDupAurumFullMatches.longregistration<-
  NonDupAurumFullMatches[which(NonDupAurumFullMatches$difference_year>10),]

#keep all gt10 definitions
gt10_idx<-which(grepl( "_gt10" , colnames( NonDupAurumFullMatches.longregistration ) ) )

#rearrange generated demographic variables ( (lines 298:321 are created demographic variables))
NonDupAurumFullMatches.conditionsSens10_32<-NonDupAurumFullMatches.longregistration[,c(1:25,298:321,gt10_idx)]

#remove "date"  fields
d_idx<-which(grepl( "_date" , colnames( NonDupAurumFullMatches.conditionsSens10_32 ) ) )
NonDupAurumFullMatches.conditionsSens10_32<-NonDupAurumFullMatches.conditionsSens10_32[,-d_idx]

#remove some medication columns and other columns that are not needed for these analyses
# remove columns 70-75, 89-116
NonDupAurumFullMatches.conditionsSens10_32<-NonDupAurumFullMatches.conditionsSens10_32[,-c(70:75,89:116)]

# now run the analysis
errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",110)
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double())
DS<-NonDupAurumFullMatches.conditionsSens10_32
DF<-ProcessConditionsMedications(DF,DS,49,errorsDS,warningDS)
write.table(DF,file=paste0(WD,"NonDupAurumFullMatchesSens10_32.Apr182024.txt"),row.names = F,sep="\t")

