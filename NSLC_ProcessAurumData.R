# this file:
#     1) reads in initial SAS file
#     2) creates some additional variables
#     3) removes individuals who are not 30-89 and anyone with a "smoking date"
#     4) Removes individuals who also appeared in GOLD (CPRD gave us a text file of deduplicated persons)
#.    5) removes cases without controls in de-duplication process
#     6) creates Table 1 and output to text file that is easily imported into excel
#     7) writes out R file for use in primary analyses (AurumPrimaryAnalyses.R); 
#               Sensitivty analyses (Associations_BMISES.R)
# 


library(sqldf)
library(sas7bdat) 
library(haven)
library(reshape2)
library(survival)
library(epitools)
# read in sas file. set your own working directory
WD<-"./"
setwd(WD)
source("ProcessConditionalLogitApr2021.R")

############################ read in file that IMS gave us
medications_file<-read_sas("covariate_file_051022.sas7bdat")

# MOB is missing for nearly all the population. set to Jan 1, YOB for those missing MOB
medications_file$month_d<-ifelse(is.na(medications_file$MOB), 1, medications_file$MOB)
medications_file$DOB <- as.Date(paste(medications_file$YOB, medications_file$month_d, 1,sep="-"), "%Y-%m-%d")

### create index age variable
indexage<-floor(as.numeric(difftime(as.Date(medications_file$INDEXDATE,origin ="1960-01-01"),medications_file$DOB,units='days'))/365.25)
medications_file$indexage<-indexage

#create difference variable(time between diagnosis date and registration start)
medications_file$difference<-as.numeric(difftime(as.Date(medications_file$INDEXDATE,origin="1960-01-01"),
                                                 as.Date(medications_file$REGSTARTDATE,origin="1960-01-01"),
                                                 units='days'))


#create number of years prior to indexdate
medications_file$difference_year<-medications_file$difference/365.25

# exclude individuals <30 and >89
medications_file_ageproper<-medications_file[which(medications_file$indexage>=30 & medications_file$indexage<90),]

options(digits=3)

######## remove anyone who had a smoking diagnosis after selection date
smokers.idx<-which(!is.na(medications_file_ageproper$first_date_of_smoking))
medications_file_ageproper<-medications_file_ageproper[-smokers.idx,]

#create INDEXYEAR variable
medications_file_ageproper$INDEXYEAR <- as.numeric(format(as.Date(medications_file_ageproper$INDEXDATE,origin ="1960-01-01"),'%Y'))

#create REGISTRATIONYEAR variable
medications_file_ageproper$REGISTRATIONYEAR <- as.numeric(format(as.Date(medications_file_ageproper$REGSTARTDATE,origin ="1960-01-01"),'%Y'))

# create categories of registration
medications_file_ageproper$CRD_CAT_STRING<-ifelse(medications_file_ageproper$REGISTRATIONYEAR < 1921, "1908-1920", 
                                                  ifelse(medications_file_ageproper$REGISTRATIONYEAR < 1931, "1921-1930",
                                                         ifelse(medications_file_ageproper$REGISTRATIONYEAR < 1941, "1931-1940",
                                                                ifelse(medications_file_ageproper$REGISTRATIONYEAR < 1951, "1941-1950",
                                                                       ifelse(medications_file_ageproper$REGISTRATIONYEAR < 1961, "1951-1960",
                                                                              ifelse(medications_file_ageproper$REGISTRATIONYEAR < 1971, "1961-1970", 
                                                                                     ifelse(medications_file_ageproper$REGISTRATIONYEAR < 1981, "1971-1980",
                                                                                            ifelse(medications_file_ageproper$REGISTRATIONYEAR < 1991, "1981-1990",
                                                                                                   ifelse(medications_file_ageproper$REGISTRATIONYEAR < 2001, "1991-2000",
                                                                                                          ifelse(medications_file_ageproper$REGISTRATIONYEAR < 2011, "2001-2010","2011-2019"))))))))))

# create categories of diagnosis year
medications_file_ageproper$INDEXYEAR_CAT_STRING<-ifelse(medications_file_ageproper$INDEXYEAR < 1995, "1989-1994", 
                                                        ifelse(medications_file_ageproper$INDEXYEAR < 2001, "1995-2000",
                                                               ifelse(medications_file_ageproper$INDEXYEAR < 2007, "2001-2006",
                                                                      ifelse(medications_file_ageproper$INDEXYEAR < 2013, "2007-2012","2013-2018"))))


####################################################################
# if bmi assessment date > 15 years ago then set to missing. create new variable
# below is the obesity flag we were given. Need to update and make those with measurements > 15 YEARS ago missing

medications_file_ageproper$BMI_DIFFYR<-ifelse(is.na(medications_file_ageproper$OBESITY), NA, 
                                              medications_file_ageproper$OBESITY_TIME_TO_ENDOFFUP/365.25)

medications_file_ageproper$BMI_12M_Rev1<-ifelse(is.na(medications_file_ageproper$BMI_DIFFYR), NA, 
                                                ifelse(medications_file_ageproper$BMI_DIFFYR > 15, NA, 
                                                       medications_file_ageproper$OBESITY))


####################################################################.   #####################################
# if alcohol assessment date > 15 years ago then set to missing. create new variable

medications_file_ageproper$ALCOHOL_DIFFYR<-ifelse(is.na(medications_file_ageproper$ALCOHOL_STATUS_CAT), NA, 
                                                  medications_file_ageproper$ALCOHOL_STATUS_TIME_TO_ENDOFFUP/365.25)


medications_file_ageproper$ALCOHOL_STATUS_1YRS_Rev1<-ifelse(is.na(medications_file_ageproper$ALCOHOL_DIFFYR), NA, 
                                                            ifelse(medications_file_ageproper$ALCOHOL_DIFFYR > 15, NA, 
                                                                   medications_file_ageproper$ALCOHOL_STATUS_CAT))

##############################################################################
# update height variable if > 120 or < 220 then / 100. some were poorly input

medications_file_ageproper$HEIGHT_M_Rev1<-ifelse(is.na(medications_file_ageproper$HEIGHT_IN_M), NA, 
                                                 ifelse( (medications_file_ageproper$HEIGHT_IN_M >120 & medications_file_ageproper$HEIGHT_IN_M<220), medications_file_ageproper$HEIGHT_IN_M/100,
                                                         medications_file_ageproper$HEIGHT_IN_M))

hist(medications_file_ageproper$HEIGHT_M_Rev1)


###################.          create more variables for table 1
#COPD variable is for checking if estimates in controls matched expected % of never smoking COPD cases in population
medications_file_ageproper$copd.f<-ifelse(medications_file_ageproper$copd_lt10==0, "COPD: No","COPD: Yes")
medications_file_ageproper$Case_Status.f<-ifelse(medications_file_ageproper$case_status==0,"Control","Case")
medications_file_ageproper$Sex.f<-ifelse(medications_file_ageproper$GENDER==1,"Male","Female")
medications_file_ageproper$BMI_M12_ALL.f<-ifelse(medications_file_ageproper$BMI_12M_Rev1==1, "Underweight",
                                                 ifelse(medications_file_ageproper$BMI_12M_Rev1==2, "Normal",
                                                        ifelse(medications_file_ageproper$BMI_12M_Rev1==3, "Overweight",
                                                               ifelse(medications_file_ageproper$BMI_12M_Rev1==4, "Obese",
                                                                      ifelse(medications_file_ageproper$BMI_12M_Rev1==5, "Morbidly Obese", "Missing")))))

medications_file_ageproper$ALCOHOL_STATUS_MR_1YRS.f<-ifelse(medications_file_ageproper$ALCOHOL_STATUS_1YRS_Rev1==0,"Never",
                                                            ifelse(medications_file_ageproper$ALCOHOL_STATUS_1YRS_Rev1==1,"Former",
                                                                   ifelse(medications_file_ageproper$ALCOHOL_STATUS_1YRS_Rev1==2,"Current","Missing")))

medications_file_ageproper$REGION.f<-ifelse(medications_file_ageproper$REGION==1, "North East",
                                            ifelse(medications_file_ageproper$REGION==2, "North West",
                                                   ifelse(medications_file_ageproper$REGION==3, "Yorkshire & The Humber",
                                                          ifelse(medications_file_ageproper$REGION==4, "East Midlands",
                                                                 ifelse(medications_file_ageproper$REGION==5, "West Midlands",
                                                                        ifelse(medications_file_ageproper$REGION==6, "East of England",
                                                                               ifelse(medications_file_ageproper$REGION==7, "South West",
                                                                                      ifelse(medications_file_ageproper$REGION==8, "South Central",
                                                                                             ifelse(medications_file_ageproper$REGION==9, "London",
                                                                                                    ifelse(medications_file_ageproper$REGION==10, "South East Coast",
                                                                                                           ifelse(medications_file_ageproper$REGION==11, "Northern Ireland",
                                                                                                                  ifelse(medications_file_ageproper$REGION==12, "Scotland",
                                                                                                                         ifelse(medications_file_ageproper$REGION==13, "Wales","Unknown")))))))))))))


########################################################################################################################################################
# #######
# create dataset of individuals who are not duplicated in GOLD (28416)
# and also have any linkage
# See below for specifics
#

# from log file that CPRD gave me
#Total number of patients	28416	Total number of patients after de-duplication
#Total number of unique patients	28416	
#Number of acceptable patients	28416	
#Number of patients eligible for linkage	28416	
#Numbers eligible for linkage to HES	28416	
#load Aurum linkable and dedupped data
#18_231R_linkage_eligibility_gold.txt file that contains all eligibility.

Aurum_All_Linked_File<-'LinkedData/Results/Aurum_linked/Final/18_231R_linkage_eligibility_aurum.txt'
Aurum_All_Linked<-read.delim(file=Aurum_All_Linked_File)


AurumLinked_IMD<-'LinkedData/Results/Aurum_linked/Final/practice_imd_18_231R.txt'
Aurum_IMD<-read.delim(file=AurumLinked_IMD)

########### ###########link only keep non-duplicates in analytic data file #############################
linked_query<-"SELECT distinct a.* FROM medications_file_ageproper as a where a.patid in (Select distinct patid from Aurum_All_Linked)";
NonDuplicatedAurum<-sqldf(linked_query,method = "name__class")

##################################################################################
# FINAL population is below! 
#         Need to remove controls who are missing a matched case!!!!!!
#         first identify cases that don't have controls (SETID)

case_idx<-which(NonDuplicatedAurum$case_status==1)

#cases below
case_nondup<-NonDuplicatedAurum[case_idx,]
# controls 
ctrl_nondup<-NonDuplicatedAurum[-case_idx,]

#only keep individuals with a match!!!!!!
NonDupAurumFullMatches_idx<-which(NonDuplicatedAurum$SETID%in%case_nondup$SETID)

NonDupAurumFullMatches<-NonDuplicatedAurum[NonDupAurumFullMatches_idx,]

NonDupAurumFullMatches$AGE_CAT_STRING6<- ifelse(NonDupAurumFullMatches$indexage < 40, "30-39",
                                                ifelse(NonDupAurumFullMatches$indexage < 50, "40-49",
                                                       ifelse(NonDupAurumFullMatches$indexage < 60, "50-59",
                                                              ifelse(NonDupAurumFullMatches$indexage < 70, "60-69",
                                                                     ifelse(NonDupAurumFullMatches$indexage < 80, "70-79", "80-89")))))

NonDupAurumFullMatches$CRD_CAT_STRING5<-ifelse(NonDupAurumFullMatches$REGISTRATIONYEAR < 1981, "1980 and earlier", 
                                               ifelse(NonDupAurumFullMatches$REGISTRATIONYEAR < 1991, "1981-1990",
                                                      ifelse(NonDupAurumFullMatches$REGISTRATIONYEAR < 2001, "1991-2000",
                                                             ifelse(NonDupAurumFullMatches$REGISTRATIONYEAR < 2011, "2001-2010","2011-2018"))))

# create categories of diagnosis year
NonDupAurumFullMatches$INDEXYEAR_CAT_STRING5<-ifelse(NonDupAurumFullMatches$INDEXYEAR < 1995, "1989-1994", 
                                                     ifelse(NonDupAurumFullMatches$INDEXYEAR < 2001, "1995-2000",
                                                            ifelse(NonDupAurumFullMatches$INDEXYEAR < 2007, "2001-2006",
                                                                   ifelse(NonDupAurumFullMatches$INDEXYEAR < 2013, "2007-2012","2013-2019"))))

#Hack to keep things in correct order in table 1
#group obese and morbidly obese together in aurum
NonDupAurumFullMatches$Sex.f<-ifelse(NonDupAurumFullMatches$GENDER==1,"1.Male","0.Female")
NonDupAurumFullMatches$BMI_M12_ALL.f<-ifelse(is.na(NonDupAurumFullMatches$BMI_12M_Rev1),"5.Missing",
                                             ifelse(NonDupAurumFullMatches$BMI_12M_Rev1==1, "0.Underweight",
                                                    ifelse(NonDupAurumFullMatches$BMI_12M_Rev1==2, "1.Normal",
                                                           ifelse(NonDupAurumFullMatches$BMI_12M_Rev1==3, "2.Overweight","3.Obese")))) #combine obese and morbidly obese

NonDupAurumFullMatches$ALCOHOL_STATUS_MR_1YRS.f<-ifelse(NonDupAurumFullMatches$ALCOHOL_STATUS_1YRS_Rev1==0,"2.Never",
                                                        ifelse(NonDupAurumFullMatches$ALCOHOL_STATUS_1YRS_Rev1==1,"1.Former",
                                                               ifelse(NonDupAurumFullMatches$ALCOHOL_STATUS_1YRS_Rev1==2,"0.Current","3.Missing")))


label(NonDupAurumFullMatches$REGION.f)<-"Region, n (%)"
label(NonDupAurumFullMatches$Case_Status.f) <- "Case status"
label(NonDupAurumFullMatches$INDEXYEAR_CAT_STRING5)<-"Year of diagnosis/selection, n (%)"
label(NonDupAurumFullMatches$Sex.f)<-"Sex, n(%)"
label(NonDupAurumFullMatches$AGE_CAT_STRING6)<-"Age (years) at selection/diagnosis, n (%)"
label(NonDupAurumFullMatches$CRD_CAT_STRING5)<-"Year of practice registration, n (%)"

label(NonDupAurumFullMatches$BMI_M12_ALL.f)<-"Body Mass Index (kg/m2), n (%)"  
label(NonDupAurumFullMatches$ALCOHOL_STATUS_MR_1YRS.f)<-"Alcohol Consumption, n (%)"  
label(NonDupAurumFullMatches$difference_year)<-"Duration of registration in practice, year (IQR)"
#label(NonDuplicatedAurum$IMD)<-"Practice level index of deprivation"

labels <- list(
  variables=list(sex=render.varlabel(NonDupAurumFullMatches$Sex.f),
                 crd=render.varlabel(NonDupAurumFullMatches$CRD_CAT_STRING5),
                 indexyear=render.varlabel(NonDupAurumFullMatches$INDEXYEAR_CAT_STRING5),
                 age=render.varlabel(NonDupAurumFullMatches$AGE_CAT_STRING6),
                 bmi=render.varlabel(NonDupAurumFullMatches$BMI_M12_ALL.f),
                 bmidiff=render.varlabel(NonDupAurumFullMatches$BMI_DIFFYR_Rev1),
                 alcohol=render.varlabel(NonDupAurumFullMatches$ALCOHOL_STATUS_MR_1YRS.f),
                 enrollmenttime=render.varlabel(NonDupAurumFullMatches$difference_year)))

# main table 1
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=3), c("",
                                                           "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%2.1f %%)", FREQ, PCT))))
}
aurum_table1<-table1(~ AGE_CAT_STRING6 + Sex.f + CRD_CAT_STRING5 +INDEXYEAR_CAT_STRING5 + difference_year +
                       BMI_M12_ALL.f +  ALCOHOL_STATUS_MR_1YRS.f + REGION.f| Case_Status.f, 
                     topclass = "Rtable1-zebra",
                     render.categorical=my.render.cat,render.continuous=c(.="Mean (SD)", .="Median [Min,Q1,Q3 Max]"),
                     data=NonDupAurumFullMatches)

###### write out table to something easily imported into excel
library(rvest)
tmp1<-read_html(aurum_table1)
aurum_table1.table<-html_table(html_nodes(tmp1, "table")[[1]])

#write table 1
write.table(aurum_table1.table,file=paste0(WD,"AurumDeDuppedTable1_May2024.txt"),sep="\t")
#save file for analyses

save(NonDupAurumFullMatches, file = paste0(WD,"AurumDeDupped.May2024.RData"))


