
# this file creates a text file that will be read into SAS for hierarchy analysis
#             level 3 conditions as level 2 (e.g., autoimmune-lupus)
#             ignore all level 2 conditions which had a level 3 condition
library(sqldf)
library(sas7bdat) 
library(haven)
library(reshape2)
library(survival)

##### 
WD<-"./"
setwd(WD)
source(paste0(WD,"ProcessConditionalLogitApr2021.R"))

# load and read in primay DS if not already there; 
#       this dataset was created in FirstProcessData.R 

load( file = paste0(WD,"FullDSAgeProper.Feb2022.RData")) #what I called it after generating in FirstProcessData.R

#pull out only primary conditions
# MD June 2022 - add BMI to hierarchy file

options(digits=5)
##################################################################################
##################### ##################### ##################### ##################### 
# sensitive conditions (1-10 years); not "alt_lt10" or "_d"
lt10_idx<-which(grepl( "_lt10" , colnames( medications_file_ageproper ) ) )

Year_Idx<-which(grepl( "REGISTRATIONYEAR", colnames( medications_file_ageproper )) | 
                  grepl( "INDEXYEAR" , colnames( medications_file_ageproper )))

#printing out BMI in case adjustment needed in hierarchy analysis. It ended up being done in conditional logistic
# regresion analyses
BMI_Idx<-
  which(grepl( "BMI_M12_ALL.f", colnames( medications_file_ageproper )))

medications_file_ageproper.conditionsSens1_10<-
  medications_file_ageproper[,c(1:34,Year_Idx,BMI_Idx,lt10_idx)]

# find alternate definitions and remove
alt_lt10_idx<-
  which(grepl( "_alt_lt10" , colnames( medications_file_ageproper.conditionsSens1_10 ) ) )
medications_file_ageproper.conditionsSens1_10<-
  medications_file_ageproper.conditionsSens1_10[,-alt_lt10_idx]

###########    now remove date columns from analytical file
d_idx<-which(grepl("_d",colnames( medications_file_ageproper.conditionsSens1_10)))
medications_file_ageproper.conditionsSens1_10<-medications_file_ageproper.conditionsSens1_10[,-d_idx]

####################################################################################################################

########### 
########### CREATE Hierarchty analysis for for Infections and Inflammation Category

#identify all columns associated with the infections and inflammation category
II_idx<-which(grepl( "InfectInflam" , colnames( medications_file_ageproper.conditionsSens1_10 ) ) )

### make sure the following variables are present in the generated file:
#indexage REGION GENDER REGISTRATIONYEAR INDEXYEAR CASE_STATUS SETID 
# colums 1:7,13,17,33:34 are demographic and necessary variables
# column 169 is "InfectInflam_Other_lt10" for conditions that were II but didn't fit into any category
# do not include level 2 categories that have level 3 categories (InfectInflam_GI, InfectInflam_AI). This is a two level anaylysis
medications_file_ageproper.Sens1_10HierarchyII<-medications_file_ageproper.conditionsSens1_10[c(1:7,13,17,33:34,169,II_idx)]

#remove categories that had been merged "InfectInflam_GI_Gastr_lt10", InfectInflam_GI_NIIGC_lt10, InfectInflam_AI_other_lt10
# see processing for when the merging of subcategories occurred. they had been deemed too similar or there had been an 
# error in the read codes that was fixed.
CondRemove_Idx<-which(grepl( "InfectInflam_GI_Gastr_lt10", colnames( medications_file_ageproper.Sens1_10HierarchyII )) | 
                  grepl( "InfectInflam_GI_NIIGC_lt10" , colnames( medications_file_ageproper.Sens1_10HierarchyII))|
                  grepl( "InfectInflam_AI_other_lt10" , colnames( medications_file_ageproper.Sens1_10HierarchyII )))

medications_file_ageproper.Sens1_10HierarchyII<-medications_file_ageproper.Sens1_10HierarchyII[,-CondRemove_Idx]
## write out infections and inflammation hierarchy for analysis
#updated June 2022 in include BMI
write.table(medications_file_ageproper.Sens1_10HierarchyII,file=paste0(WD,"HierarchyInfectInflamDSBMIJune2022.txt"),col.names=TRUE,sep="\t")

#######################################################################################################################
##################### ##################### ##################### ##################### ##################### ##################### 
# MD Jan 2022
# this may not be necessary, but should do regardless for supplemental information
# sensitive conditions (10-32 years); for Infections and Inflammation
II_idx<-which(grepl( "InfectInflam" , colnames( medications_file_ageproper ) ) )

Year_Idx<-which(grepl( "REGISTRATIONYEAR", colnames( medications_file_ageproper )) | 
                  grepl( "INDEXYEAR" , colnames( medications_file_ageproper )) |
                  grepl( "difference_year" , colnames( medications_file_ageproper )))


medications_file_ageproper.conditionsSens10_32<-medications_file_ageproper[,c(1:32,Year_Idx,II_idx)]


# remove lt10
lt10_idx<-which(grepl( "_lt10" , colnames( medications_file_ageproper.conditionsSens10_32 ) ) )
medications_file_ageproper.conditionsSens10_32<-medications_file_ageproper.conditionsSens10_32[,-lt10_idx]

# remove dates _d
d_idx<-which(grepl( "_d" , colnames( medications_file_ageproper.conditionsSens10_32 ) ) )
medications_file_ageproper.conditionsSens10_32<-medications_file_ageproper.conditionsSens10_32[,-d_idx]

# remove alternate (more specific) definition _alt_gt10
alt_gt10_idx<-which(grepl( "_alt_gt10" , colnames( medications_file_ageproper.conditionsSens10_32 ) ) )
medications_file_ageproper.conditionsSens10_32<-medications_file_ageproper.conditionsSens10_32[,-alt_gt10_idx]

CondRemove_Idx<-which(grepl( "InfectInflam_GI_Gastr_gt10", colnames( medications_file_ageproper.conditionsSens10_32 )) | 
                        grepl( "InfectInflam_GI_NIIGC_gt10" , colnames( medications_file_ageproper.conditionsSens10_32))|
                        grepl( "InfectInflam_AI_other_gt10" , colnames( medications_file_ageproper.conditionsSens10_32 )))

medications_file_ageproper.conditionsSens10_32<-medications_file_ageproper.conditionsSens10_32[,-CondRemove_Idx]

# now restrict to persons with >=10 years CPRD enrollment
medications_file_ageproper.conditionsSens10_32<-
  medications_file_ageproper.conditionsSens10_32[which(medications_file_ageproper.conditionsSens10_32$difference_year>10),]

medications_file_ageproper.Sens10_32HierarchyII<-
  medications_file_ageproper.conditionsSens10_32[c(1:7,13,17,32:33,37:90)]

## write out infections and inflammation hierarchy for analysis
write.table(medications_file_ageproper.Sens10_32HierarchyII,file=paste0(WD,"HierarchyInfectInflam10_32DSJan2022.txt"),col.names=TRUE,sep="\t")


