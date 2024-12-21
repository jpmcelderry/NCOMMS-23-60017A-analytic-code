# this file runs all conditional logistic associations
#                 it uses ProcessConditionsMedications() function that is found in
#                 ProcessConditionalLogitApr2021.R

# outputs several text files that form the basis of Table 2 and 
# supplementary tables of all conditional logistic regression results.

library(sqldf)
library(sas7bdat) 
library(haven)
library(reshape2)
library(survival)
library(epitools)

#set working directory
WD<-"./"
setwd(WD)
source(paste0(WD,"NSLC/ProcessConditionalLogitApr2021.R",sep=""))

load( file = paste0(WD,"FullDSAgeProper.Feb2022.RData"))

#write table 1
# main table 1 (GOLD)
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=3), c("",
                                                           "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%2.1f %%)", FREQ, PCT))))
}
gold_table1<-table1(~ AGE_CAT_STRING6 + Sex.f + CRD_CAT_STRING5 +INDEXYEAR_CAT_STRING5 + difference_year +
                      BMI_M12_ALL.f +  ALCOHOL_STATUS_MR_1YRS.f + REGION.f| Case_Status.f, 
                    topclass = "Rtable1-zebra",
                    render.categorical=my.render.cat,render.continuous=c(.="Mean (SD)", .="Median [Min,Q1,Q3 Max]"),
                    data=medications_file_ageproper)

###### ########################################################################
# write out table to something easily imported into excel
library(rvest)
tmp1<-read_html(gold_table1)
gold_table1.table<-html_table(html_nodes(tmp1, "table")[[1]])

#write table 1
write.table(gold_table1.table,file=paste0(WD,"GoldTable1_May2024.txt"),sep="\t")


##################################################################################
##################### ##################### ##################### ##################### 
# process sensitive 1-10 conditions. This is the primary analysis of associations between
# conditions diagnosed 1-10 years before selection date (index date) and LCINS

# identify columns corresponding to the definitions of associations we want to estimate
lt10_idx<-which(grepl( "_lt10" , colnames( medications_file_ageproper ) ) )

#rearrange some columns, move created variables next to other demographic variables
Year_Idx<-which(grepl( "REGISTRATIONYEAR", colnames( medications_file_ageproper )) | 
                  grepl( "INDEXYEAR" , colnames( medications_file_ageproper )))
medications_file_ageproper.conditionsSens1_10<-medications_file_ageproper[,c(1:32,Year_Idx,lt10_idx)]

# remove alternate definitions (2+ DX in time frame)
alt_lt10_idx<-which(grepl( "_alt_lt10" , colnames( medications_file_ageproper.conditionsSens1_10 ) ) )
medications_file_ageproper.conditionsSens1_10<-medications_file_ageproper.conditionsSens1_10[,-alt_lt10_idx]

#now remove date columns because they are not needed for the analyses
d_idx<-which(grepl("_d",colnames( medications_file_ageproper.conditionsSens1_10)))
medications_file_ageproper.conditionsSens1_10<-medications_file_ageproper.conditionsSens1_10[,-d_idx]

####overall sensitive condition associations. Primary. Subset is table 2 ("main categories in hierarchy")
errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",125)
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double())
DS<-medications_file_ageproper.conditionsSens1_10

#analysis starts at column 40 and continues until the end
DF<-ProcessConditionsMedications(DF,DS,40,errorsDS,warningDS)

#create tab delimited file with ALL associations across full hierarchy.
# these raw results were used to generate much of table 2 and supplemental table2
write.table(DF,file=paste0(WD,"Conditions1-10Sens.Feb62022.txt"),row.names = F,sep="\t")

############### associations in females (GENDER==2)).

errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",125)
female_idx<-which(medications_file_ageproper.conditionsSens1_10$GENDER==2)
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double())
DS<-medications_file_ageproper.conditionsSens1_10[female_idx,]

DF<-ProcessConditionsMedications(DF,DS,40,errorsDS,warningDS)

write.table(DF,file=paste0(WD,"ConditionsFemale1-10Sens.Feb62022.txt"),row.names = F,sep="\t")

############### sensitive conditions for males (GENDER==1))

errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",125)
male_idx<-which(medications_file_ageproper.conditionsSens1_10$GENDER==1)
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double())
DS<-medications_file_ageproper.conditionsSens1_10[male_idx,]

DF<-ProcessConditionsMedications(DF,DS,40,errorsDS,warningDS)

write.table(DF,file=paste0(WD,"ConditionsMale1-10Sens.Feb62022.txt"),row.names = F,sep="\t")

###########################################################################################################
###########################################################################################################
###########################################################################################################
# ANALYSES for conditions 10-32 years (sensitive definition)
# MUST first identify people with at least 10 years of registration
medications_file_ageproper.longregistration<-medications_file_ageproper[which(medications_file_ageproper$difference_year>10),]

alt_gt10_idx<-which(grepl( "alt_gt10" , colnames( medications_file_ageproper.longregistration ) ) )
# remove all alternate definitions
medications_file_ageproper.conditionsSens10_32<-medications_file_ageproper.longregistration[,-alt_gt10_idx]

#keep all gt10 definitions
gt10_idx<-which(grepl( "_gt10" , colnames( medications_file_ageproper.conditionsSens10_32 ) ) )
medications_file_ageproper.conditionsSens10_32<-medications_file_ageproper.conditionsSens10_32[,c(1:32,gt10_idx)]

#remove "date"  fields
d_idx<-which(grepl( "_d" , colnames( medications_file_ageproper.conditionsSens10_32 ) ) )
medications_file_ageproper.conditionsSens10_32<-medications_file_ageproper.conditionsSens10_32[,-d_idx]

# remove medications
medications_file_ageproper.conditionsSens10_32<-medications_file_ageproper.conditionsSens10_32[,-c(160:260)]

## ## ## ## now run the analysis. variables of interest start on column 36
errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",125)
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double())
DS<-medications_file_ageproper.conditionsSens10_32

DF<-ProcessConditionsMedications(DF,DS,36,errorsDS,warningDS)

write.table(DF,file=paste0(WD,"Conditions10-32Sens.Feb62022.txt"),row.names = F,sep = "\t")


# now do for females (conditions 10-32 years)
###############  sensitive conditions for females (GENDER==2))

errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",125)
female_idx<-which(medications_file_ageproper.conditionsSens10_32$GENDER==2)
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double())
DS<-medications_file_ageproper.conditionsSens10_32[female_idx,]
DF<-ProcessConditionsMedications(DF,DS,36,errorsDS,warningDS)

write.table(DF,file=paste0(WD, "ConditionsFemale10-32Sens.Feb62022.txt"),row.names = F,sep="\t")

# now do for males(conditions 10-32 years)
###############  sensitive conditions for males (GENDER==1))

errorsDS<-data.frame(error=character())
warningDS<-rep("No Error",125)
male_idx<-which(medications_file_ageproper.conditionsSens10_32$GENDER==1)
DF<-data.frame(Condition=character(),ECase=integer(), EControl=integer(),OR=double(),low=double(),high=double(),pvalue=double())
DS<-medications_file_ageproper.conditionsSens10_32[male_idx,]
dim(DS)
DF<-ProcessConditionsMedications(DF,DS,36,errorsDS,warningDS)

write.table(DF,file=paste0(WD,"ConditionsMale10-32Sens.Feb62022.txt"),row.names = F,sep="\t")


