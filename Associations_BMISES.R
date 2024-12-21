####  Associations_BMISES.R
#     SES/BMI adjusted associations. Perform in Aurum and GOLD. 
library(sqldf)
library(sas7bdat) 
library(haven)
library(reshape2)
library(survival)
library(epitools)

# SET WORKING DIRECTORY
WD<-"./"
setwd(WD)
source(paste0(WD,"NSLC/ProcessConditionalLogitApr2021.R"))

#read in R file 
load( file = paste0(WD,"AurumDeDupped.May2024.RData"))

#Estimate unadjusted and jointly adjusted associations for both BMI and SES.

##### ########### ########### ########### ########### ########### ########### ########### ######
##### ######create your own working directory to linkage files sent by CPRD
aurumlinkedWD<-"NSLC/LinkedData/Results/Aurum_linked/Final/"
Aurum_IMD<-read.delim(file=paste0(aurumlinkedWD,"practice_imd_18_231R.txt"))

#link IMD file with analytic data file. 
NonDupAurumFullMatchesTmp<-NonDupAurumFullMatches[,c(1:319)]
linkedSES_query<-"SELECT distinct b.IMD, b.IMD_pracid,a.* FROM NonDupAurumFullMatchesTmp as a
                          INNER JOIN (SELECT pracid as IMD_pracid, e2015_imd_10 as IMD from Aurum_IMD) as b
                ON a.pracid = b.IMD_pracid"
linkedAurum_SES<-sqldf(linkedSES_query,method = "name__class")

#### ###################################examine associations adjusted for SES/BMI
#### SES is in deciles and is coded as continuous variable
### BMI is a categorical variable coded as (Underweight, normal, overweight, obese (included morbidly obese), missing)
### BMI is classified as missing if it was measured > 15 years before selection

#conditional logistic regression analyses with de-duplicated Aurum population
#(first 1-10 years; then 10-32 years)
linkedAurum_SES$CASE_STATUS<-linkedAurum_SES$case_status

#demographic variables, extra demographic/case variables, validation variables 1-10 years prior
linkedAurum_SES.conditionsSens1_10.small<-linkedAurum_SES[,c(1:26,301:322,28:107)]

lt10_idx<-which(grepl( "_lt10" , colnames( linkedAurum_SES.conditionsSens1_10.small ) ) )

#now remove date columns
d_idx<-which(grepl("_date",colnames( linkedAurum_SES.conditionsSens1_10.small)))
linkedAurum_SES.conditionsSens1_10.small<-linkedAurum_SES.conditionsSens1_10.small[,-d_idx]

# run on practices that were linked to SES - adjust for SES (in deciles) and BMI (categories including missing category)
# print out associations for all conditions to validate including conditions with possible multiple etiologies (e.g., anemia, iron-related anemia )
errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",50)
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double(),OR1=double(),low1=double(),high1=double(),pvalue1=double())
DS<-linkedAurum_SES.conditionsSens1_10.small
DF<-ProcessConditions_BMISESAdjusted(DF,DS,49,errorsDS, warningDS)
write.table(DF,file=paste0(WD,"NonDupAurumFullMatchesSESBMIAdjSens1_10.Jul3.txt"),row.names = F,sep="\t")

############################################################################################################
###### now do adjusted associations for GOLD - we lose a lot of individuals because many were not linkable.
#     This was expected
############################################################################################################
load( file = paste0(WD, "FullDSAgeProper.Feb2022.RData"))
#load GOLD SES data, and link to GOLD SES
#18_231R_linkage_eligibility_gold.txt file that contains all eligibility.

#identify which GOLD participants could be linked to SES information
GOLDLinkWD<-"NSLC/LinkedData/Results/GOLD_linked/Final/"
GOLD_All_Linked<-read.delim(file=paste0(GOLDLinkWD,'18_231R_linkage_eligibility_gold.txt'))

#practice_imd_18_231R.txt contains index of multiple deprivation for practices
GOLD_IMD<-read.delim(file=paste0(GOLDLinkWD,"practice_imd_18_231R.txt"))

#link IMD file with analytic data file
linkedSES_query<-"SELECT distinct b.IMD, b.IMD_pracid,a.* FROM medications_file_ageproper as a
                          INNER JOIN (SELECT pracid as IMD_pracid, e2015_imd_10 as IMD from GOLD_IMD) as b
                ON a.pracid = b.IMD_pracid"
linkedGOLD_SES<-sqldf(linkedSES_query,method = "name__class")

# Save smaller data set with SES for another project (Network analysis)
save(linkedGOLD_SES, file = paste0(WD,"GOLD_SESLinked.Dec2024.RData"))
write.table(linkedGOLD_SES, file=paste0(WD,"GOLD_SESLinked.Dec2024.txt"),col.names=TRUE,row.names = F, sep="\t")


#process sensitive 1-10 conditions
#extract of the significant columns (associations)
lt10_idx<-which(grepl( "_lt10" , colnames( linkedGOLD_SES ) ) )
length(lt10_idx)
#518

Year_Idx<-which(grepl( "REGISTRATIONYEAR", colnames( linkedGOLD_SES )) | 
                  grepl( "INDEXYEAR" , colnames( linkedGOLD_SES )))
BMI_Idx<-which(grepl( "BMI_M12_ALL.f", colnames( linkedGOLD_SES )))
linkedGOLD_SES.conditionsSens1_10<-linkedGOLD_SES[,c(1:34,Year_Idx,BMI_Idx,lt10_idx)]

###########    now remove date columns from analytical file
d_idx<-which(grepl("_d",colnames( linkedGOLD_SES.conditionsSens1_10)))
linkedGOLD_SES.conditionsSens1_10<-linkedGOLD_SES.conditionsSens1_10[,-d_idx]

# remove alternate definitions
alt_lt10_idx<-which(grepl( "_alt_lt10" , colnames( linkedGOLD_SES.conditionsSens1_10 ) ) )
linkedGOLD_SES.conditionsSens1_10<-linkedGOLD_SES.conditionsSens1_10[,-alt_lt10_idx]

#just run on significant associations in gold
sigcols_gold_idx<-
  which(colnames(linkedGOLD_SES.conditionsSens1_10)%in%c("InfectInflam_lt10",
                                                         "InfectInflam_TB_lt10",
                                                         "InfectInflam_Flu_lt10",
                                                         "InfectInflam_Resp_lt10",
                                                         "InfectInflam_UResp_lt10",
                                                         "InfectInflam_AI_lt10",
                                                         "InfectInflam_AI_Psor_lt10",
                                                         "InfectInflam_AI_other_lt10",
                                                         "InfectInflam_AI_DMT1_lt10",
                                                         "InfectInflam_COPD_lt10",
                                                         "Anemia_lt10",
                                                         "InfectInflam_AI_Lupus_lt10",
                                                         "InfectInflam_AI_oth2_lt10",
                                                         "InfectInflam_AI_Rose_lt10",
                                                         "InfectInflam_GI_Ulc_lt10",
                                                         "InfectInflam_GI_Gastr_lt10",
                                                         "InfectInflam_GI_GERD_lt10",
                                                         "InfectInflam_GI_NIIGC_lt10",
                                                         "InfectInflam_GI_GastrNIIGC_lt10"
  ))

linkedGOLD_SES.conditionsSens1_10.small<-linkedGOLD_SES.conditionsSens1_10[,c(1:38,sigcols_gold_idx)]

### rerun analyses on subset that is linkable
# run on practices that were linked to SES - adjust for SES and BMI
errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",50)
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double(),OR1=double(),low1=double(),high1=double(),pvalue1=double())
DS<-linkedGOLD_SES.conditionsSens1_10.small

#conditions of interest start on column 44
DF<-ProcessConditions_BMISESAdjusted(DF,DS,44,errorsDS, warningDS)

#this is a supplemental table
write.table(DF,file=paste0(WD,"GOLDLinkedSESBMIAdjSens1_10.Jul2024.txt"),row.names = F,sep="\t")


