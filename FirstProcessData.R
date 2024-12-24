# #######################. FirstProcessData.R # #######################
# 1) create/refine important variables for Table 1 
# 2) create new dataset that only includes individuals who are 39-89 and do not have evidence of smoking after selection.
    # A few people were identified with a smoking date after selection. It is highly unlikely that someone would begin
    # smoking >= 30 years of age
# 3) update variables like BMI/alcohol that had measurements > 15 years before selection/diagnosis

# input SAS dataset
# output R dataset used for later analyses
# UPDATE FEB 2022 - add in non-categorized subcategories for II, II-AI , II-GI 
      # These variables are necessary for Hierarchical analyses (CleanHierarchy.R) and 

# Primary dataset saved on final line

library(sqldf)
library(sas7bdat) 
library(haven)
library(reshape2)
library(epitools)

WD<-"./"
setwd(WD)

# read in original SAS data file.
# new R file is called medications_file as it also contains medication use information.
# final R dataset written at end of this file will be used in analyses
medications_file<-read_sas("covariate_file_020922.sas7bdat") 
nsmall=2
# create age grouping
medications_file$AGE_CAT_STRING<-
  ifelse(medications_file$indexage < 30, "< 30", 
         ifelse(medications_file$indexage < 40, "30-39",
                ifelse(medications_file$indexage < 50, "40-49",
                       ifelse(medications_file$indexage < 60, "50-59",
                              ifelse(medications_file$indexage < 70, "60-69",
                                     ifelse(medications_file$indexage < 80, "70-79", 
                                            ifelse(medications_file$indexage < 86, "80-85",
                                                   ifelse(medications_file$indexage < 90, "86-89", "90+"))))))))

# then create a UTSEntry variable as in above
medications_file$UTSEntry<-ifelse(is.na(medications_file$UTS),0,
                                  ifelse(medications_file$UTS>medications_file$CRD,0,1))

medications_file$UTSIndex<-ifelse(is.na(medications_file$UTS),0,
                                  ifelse(medications_file$UTS>medications_file$INDEXDATE,0,1))

#how long the person was registered in the primary care clinic before the diagnosis or selection (indexdate)
#CRD is the registration date
difference <- as.numeric(difftime(medications_file$INDEXDATE,medications_file$CRD,units='days'))
medications_file$difference<-as.numeric(difftime(medications_file$INDEXDATE,medications_file$CRD,units='days'))


## restrict to persons >=30, < 90. Rename R dataset to reflect new age range
medications_file_ageproper<-medications_file[which(medications_file$indexage>=30 & medications_file$indexage<90),]
options(digits=3)

#remove the small percentage of people who had a smoking date identified later than selection or diagnosis
smokers.idx<-which(!is.na(medications_file_ageproper$first_date_of_smoking))
medications_file_ageproper<-medications_file_ageproper[-smokers.idx,]


#format index date and practice registration years
 medications_file_ageproper$INDEXYEAR <- as.numeric(format(medications_file_ageproper$INDEXDATE,'%Y'))
 medications_file_ageproper$REGISTRATIONYEAR <- as.numeric(format(medications_file_ageproper$CRD,'%Y'))
 
 # create categories of registration
 medications_file_ageproper$CRD_CAT_STRING<-
   ifelse(medications_file_ageproper$REGISTRATIONYEAR < 1921, "1908-1920", 
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
 medications_file_ageproper$INDEXYEAR_CAT_STRING<-
   ifelse(medications_file_ageproper$INDEXYEAR < 1995, "1989-1994", 
          ifelse(medications_file_ageproper$INDEXYEAR < 2001, "1995-2000",
                 ifelse(medications_file_ageproper$INDEXYEAR < 2007, "2001-2006",
                        ifelse(medications_file_ageproper$INDEXYEAR < 2013, "2007-2012","2013-2018"))))
 
#create number of years registered in clinic prior to indexdate
medications_file_ageproper$difference_year<-medications_file_ageproper$difference/365.25


#create variable in years to indicate when BMI and alcohol consumption measurements were taken
medications_file_ageproper$BMI_DIFFYR<-ifelse(is.na(medications_file_ageproper$BMI_M12_ALL_DIFF), NA, 
                                              medications_file_ageproper$BMI_M12_ALL_DIFF/365.25)

medications_file_ageproper$ALCOHOL_DIFFYR<-ifelse(is.na(medications_file_ageproper$ALCOHOL_DIFF_1YRS), NA, 
                                              medications_file_ageproper$ALCOHOL_DIFF_1YRS/365.25)

#if BMI or alcohol use was assessed > 15 years before selection/diagnosis, consider it missing
medications_file_ageproper$BMI_DIFFYR_Rev1<-ifelse(is.na(medications_file_ageproper$BMI_DIFFYR) | medications_file_ageproper$BMI_DIFFYR > 15, NA, 
                                              medications_file_ageproper$BMI_DIFFYR)

medications_file_ageproper$ALCOHOL_DIFFYR_Rev1<-ifelse(is.na(medications_file_ageproper$ALCOHOL_DIFFYR) | medications_file_ageproper$ALCOHOL_DIFFYR > 15, NA, 
                                                   medications_file_ageproper$ALCOHOL_DIFFYR)
########################################################################################################################################
# if bmi assessment date > 15 years ago then set to missing. create new variable
medications_file_ageproper$BMI_12M_Rev1<-ifelse(is.na(medications_file_ageproper$BMI_DIFFYR), NA, 
                                                ifelse(medications_file_ageproper$BMI_DIFFYR > 15, NA, 
                                                       medications_file_ageproper$BMI_M12_ALL))


########################################################################################################################################
# if alcohol assessment date > 15 years ago then set to missing. create new variable
medications_file_ageproper$ALCOHOL_STATUS_1YRS_Rev1<-ifelse(is.na(medications_file_ageproper$ALCOHOL_DIFFYR), NA, 
                                                ifelse(medications_file_ageproper$ALCOHOL_DIFFYR > 15, NA, 
                                                       medications_file_ageproper$ALCOHOL_STATUS_MR_1YRS))

####################################################################
# update height variable if > 120 or < 220 then / 100

medications_file_ageproper$HEIGHT_M_Rev1<-ifelse(is.na(medications_file_ageproper$HEIGHT_M), NA, 
                                                 ifelse( (medications_file_ageproper$HEIGHT_M >120 & medications_file_ageproper$HEIGHT_M <220), medications_file_ageproper$HEIGHT_M/100,
                                                         medications_file_ageproper$HEIGHT_M))

# create nice table 1 using table1 library
# One level of stratification
library(table1)
medications_file_ageproper$REGION.f<-
  ifelse(medications_file_ageproper$REGION==1, "North East",
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


medications_file_ageproper$AGE_CAT_STRING6<-
  ifelse(medications_file_ageproper$indexage < 40, "30-39",
         ifelse(medications_file_ageproper$indexage < 50, "40-49",
                ifelse(medications_file_ageproper$indexage < 60, "50-59",
                       ifelse(medications_file_ageproper$indexage < 70, "60-69",
                              ifelse(medications_file_ageproper$indexage < 80, "70-79", "80-89")))))




medications_file_ageproper$CRD_CAT_STRING5<-
  ifelse(medications_file_ageproper$REGISTRATIONYEAR < 1981, "1980 and earlier", 
         ifelse(medications_file_ageproper$REGISTRATIONYEAR < 1991, "1981-1990",
                ifelse(medications_file_ageproper$REGISTRATIONYEAR < 2001, "1991-2000",
                       ifelse(medications_file_ageproper$REGISTRATIONYEAR < 2011, "2001-2010","2011-2018"))))

# create categories of diagnosis year
medications_file_ageproper$INDEXYEAR_CAT_STRING5<-
  ifelse(medications_file_ageproper$INDEXYEAR < 1995, "1989-1994", 
         ifelse(medications_file_ageproper$INDEXYEAR < 2001, "1995-2000",
                ifelse(medications_file_ageproper$INDEXYEAR < 2007, "2001-2006",
                       ifelse(medications_file_ageproper$INDEXYEAR < 2013, "2007-2012","2013-2019"))))


medications_file_ageproper$UTSIndex.f<-ifelse(medications_file_ageproper$UTSIndex==0,"No","Yes")
medications_file_ageproper$UTSEntry.f<-ifelse(medications_file_ageproper$UTSEntry==0,"No","Yes")
medications_file_ageproper$copd.f<-ifelse(medications_file_ageproper$InfectInflam_COPD_lt10==0, "COPD: No","COPD: Yes")
medications_file_ageproper$Case_Status.f<-ifelse(medications_file_ageproper$CASE_STATUS==0,"Control","Case")
medications_file_ageproper$Sex.f<-ifelse(medications_file_ageproper$GENDER==1,"1.Male","0.Female")
medications_file_ageproper$BMI_M12_ALL.f<-
  ifelse(is.na(medications_file_ageproper$BMI_12M_Rev1),"Missing",
         ifelse(medications_file_ageproper$BMI_12M_Rev1==1, "0.Underweight",
                ifelse(medications_file_ageproper$BMI_12M_Rev1==2, "1.Normal",
                       ifelse(medications_file_ageproper$BMI_12M_Rev1==3, "2.Overweight",
                              ifelse(medications_file_ageproper$BMI_12M_Rev1==4, "3.Obese","Missing")))))

medications_file_ageproper$ALCOHOL_STATUS_MR_1YRS.f<-
  ifelse(medications_file_ageproper$ALCOHOL_STATUS_1YRS_Rev1==0,"2.Never",
         ifelse(medications_file_ageproper$ALCOHOL_STATUS_1YRS_Rev1==1,"1.Former",
                ifelse(medications_file_ageproper$ALCOHOL_STATUS_1YRS_Rev1==2,"0.Current","3.Missing")))

################## create labels for Table 1
label(medications_file_ageproper$REGION.f)<-"Region, n (%)"
label(medications_file_ageproper$Case_Status.f) <- "Case status"
label(medications_file_ageproper$copd.f)<-"COPD (1-10 yrs prior to selection)"
label(medications_file_ageproper$CRD_CAT_STRING5)<-"Year of practice registration, n (%)"

label(medications_file_ageproper$Sex.f)<-"Sex, n(%)"
label(medications_file_ageproper$indexage)<-"Age at Diagnosis/Selection"
label(medications_file_ageproper$INDEXYEAR_CAT_STRING5)<-"Year of diagnosis/selection, n (%)"
label(medications_file_ageproper$AGE_CAT_STRING6)<-"Age (years) at selection/diagnosis, n (%)"
label(medications_file_ageproper$BMI_M12_ALL.f)<-"Body Mass Index (kg/m2), n (%)"  
label(medications_file_ageproper$ALCOHOL_STATUS_MR_1YRS.f)<-"Alcohol Consumption, n (%)"  
label(medications_file_ageproper$difference_year)<-"Duration of registration in practice, year (IQR)"

labels <- list(
  variables=list(sex=render.varlabel(medications_file_ageproper$Sex.f),
                 crd=render.varlabel(medications_file_ageproper$CRD_CAT_STRING5),
                 indexyear=render.varlabel(medications_file_ageproper$INDEXYEAR_CAT_STRING5),
                 age=render.varlabel(medications_file_ageproper$AGE_CAT_STRING6),
                 bmi=render.varlabel(medications_file_ageproper$BMI_M12_ALL.f),
                 bmidiff=render.varlabel(medications_file_ageproper$BMI_DIFFYR_Rev1),
                 alcohol=render.varlabel(medications_file_ageproper$ALCOHOL_STATUS_MR_1YRS.f),
                 enrollmenttime=render.varlabel(medications_file_ageproper$difference_year)))

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

#######################################################################################
##################
#
#
# UPDATE: FEB 2022 - Related to hierarchy subcategory variables that did not have a 
# specific category. Some hierarchy categories weren't explicitly generated because they were non-specific
# Combine some specific disease groupings that were decided to be too similar to have distinct category
# Do for 1-10 and 10-32 sensitive defs

#there were read codes for autoimmune conditions (within infections and inflammation category) that were not specific autoimmune disease entities
medications_file_ageproper$InfectInflam_AI_NoCAT_lt10<-
  ifelse(medications_file_ageproper$InfectInflam_AI_lt10==1 &
           (medications_file_ageproper$InfectInflam_AI_IBD_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_RA_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_Psor_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_LThy_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_DMT1_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_Lupus_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_IDP_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_oth2_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_Alo_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_PolyR_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_Rose_lt10==0), 1, 0)

#there were read codes for inflammatory GI conditions (within infections and inflammation category) that were not specific autoimmune disease entities

medications_file_ageproper$InfectInflam_GI_NoCAT_lt10<-
  ifelse(medications_file_ageproper$InfectInflam_GI_lt10==1 &
           (medications_file_ageproper$InfectInflam_GI_Diver_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Haem_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Dys_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Ulc_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Other_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Gastr_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Ose_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Inf_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Chol_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_GERD_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_HP_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_NIIGC_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Aden_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Calc_lt10==0), 1, 0)

#create combined Gastritis and NIIGC categories. We decided that they weren't conceptually different enough
medications_file_ageproper$InfectInflam_GI_GastrNIIGC_lt10<-ifelse(medications_file_ageproper$InfectInflam_GI_Gastr_lt10==1 | medications_file_ageproper$InfectInflam_GI_NIIGC_lt10==1, 1, 0)

medications_file_ageproper$InfectInflam_Other_lt10<-
  ifelse(medications_file_ageproper$InfectInflam_lt10==1 & 
           (medications_file_ageproper$InfectInflam_Pneum_lt10==0 &
              medications_file_ageproper$InfectInflam_Heart_lt10==0 &
              medications_file_ageproper$InfectInflam_HN_lt10==0 &
              medications_file_ageproper$InfectInflam_Menin_lt10==0 &
              medications_file_ageproper$InfectInflam_Hep_lt10==0 &
              medications_file_ageproper$InfectInflam_TB_lt10==0 &
              medications_file_ageproper$InfectInflam_Herpes_lt10==0 &
              medications_file_ageproper$InfectInflam_UT_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_NoCAT_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Diver_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Haem_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Dys_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Ulc_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Other_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Gastr_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Ose_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Inf_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Chol_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_GERD_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_HP_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_NIIGC_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Aden_lt10==0 &
              medications_file_ageproper$InfectInflam_GI_Calc_lt10==0 &
              medications_file_ageproper$InfectInflam_OBGYN_lt10==0 &
              medications_file_ageproper$InfectInflam_Limb_lt10==0 &
              medications_file_ageproper$InfectInflam_Skin_lt10==0 &
              medications_file_ageproper$InfectInflam_Bone_lt10==0 &
              medications_file_ageproper$InfectInflam_Septi_lt10==0 &
              medications_file_ageproper$InfectInflam_Enceph_lt10==0 &
              medications_file_ageproper$InfectInflam_Flu_lt10==0 &
              medications_file_ageproper$InfectInflam_Malaria_lt10==0 &
              medications_file_ageproper$InfectInflam_Shingles_lt10==0 &
              medications_file_ageproper$InfectInflam_Dental_lt10==0 &
              medications_file_ageproper$InfectInflam_Resp_lt10==0 &
              medications_file_ageproper$InfectInflam_UResp_lt10==0 &
              medications_file_ageproper$InfectInflam_Arthrit_lt10==0 &
              medications_file_ageproper$InfectInflam_Gout_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_NoCAT_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_IBD_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_RA_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_Psor_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_LThy_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_DMT1_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_Lupus_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_IDP_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_oth2_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_Alo_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_PolyR_lt10==0 &
              medications_file_ageproper$InfectInflam_AI_Rose_lt10==0 &
              medications_file_ageproper$InfectInflam_COPD_lt10==0), 1, 0)

###### same as above but create variables for exposure in 10-32 years before selection
medications_file_ageproper$InfectInflam_AI_NoCAT_gt10<-
  ifelse(medications_file_ageproper$InfectInflam_AI_gt10==1 &
           (medications_file_ageproper$InfectInflam_AI_IBD_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_RA_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_Psor_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_LThy_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_DMT1_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_Lupus_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_IDP_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_oth2_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_Alo_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_PolyR_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_Rose_gt10==0), 1, 0)

medications_file_ageproper$InfectInflam_GI_NoCAT_gt10<-
  ifelse(medications_file_ageproper$InfectInflam_GI_gt10==1 &
           (medications_file_ageproper$InfectInflam_GI_Diver_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Haem_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Dys_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Ulc_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Other_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Gastr_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Ose_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Inf_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Chol_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_GERD_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_HP_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_NIIGC_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Aden_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Calc_gt10==0), 1, 0)

#create combined Gastritis and NIIGC categories.  They aren't conceptually different enough
medications_file_ageproper$InfectInflam_GI_GastrNIIGC_gt10<-
  ifelse(medications_file_ageproper$InfectInflam_GI_Gastr_gt10==1 | medications_file_ageproper$InfectInflam_GI_NIIGC_gt10==1, 1, 0)


medications_file_ageproper$InfectInflam_Other_gt10<-
  ifelse(medications_file_ageproper$InfectInflam_gt10==1 & 
           (medications_file_ageproper$InfectInflam_Pneum_gt10==0 &
              medications_file_ageproper$InfectInflam_Heart_gt10==0 &
              medications_file_ageproper$InfectInflam_HN_gt10==0 &
              medications_file_ageproper$InfectInflam_Menin_gt10==0 &
              medications_file_ageproper$InfectInflam_Hep_gt10==0 &
              medications_file_ageproper$InfectInflam_TB_gt10==0 &
              medications_file_ageproper$InfectInflam_Herpes_gt10==0 &
              medications_file_ageproper$InfectInflam_UT_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_NoCAT_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Diver_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Haem_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Dys_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Ulc_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Other_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Gastr_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Ose_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Inf_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Chol_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_GERD_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_HP_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_NIIGC_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Aden_gt10==0 &
              medications_file_ageproper$InfectInflam_GI_Calc_gt10==0 &
              medications_file_ageproper$InfectInflam_OBGYN_gt10==0 &
              medications_file_ageproper$InfectInflam_Limb_gt10==0 &
              medications_file_ageproper$InfectInflam_Skin_gt10==0 &
              medications_file_ageproper$InfectInflam_Bone_gt10==0 &
              medications_file_ageproper$InfectInflam_Septi_gt10==0 &
              medications_file_ageproper$InfectInflam_Enceph_gt10==0 &
              medications_file_ageproper$InfectInflam_Flu_gt10==0 &
              medications_file_ageproper$InfectInflam_Malaria_gt10==0 &
              medications_file_ageproper$InfectInflam_Shingles_gt10==0 &
              medications_file_ageproper$InfectInflam_Dental_gt10==0 &
              medications_file_ageproper$InfectInflam_Resp_gt10==0 &
              medications_file_ageproper$InfectInflam_UResp_gt10==0 &
              medications_file_ageproper$InfectInflam_Arthrit_gt10==0 &
              medications_file_ageproper$InfectInflam_Gout_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_NoCAT_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_IBD_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_RA_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_Psor_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_LThy_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_DMT1_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_Lupus_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_IDP_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_oth2_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_Alo_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_PolyR_gt10==0 &
              medications_file_ageproper$InfectInflam_AI_Rose_gt10==0 &
              medications_file_ageproper$InfectInflam_COPD_gt10==0), 1, 0)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# CVD 1-10 years
medications_file_ageproper$CVD_CHD_Other_lt10<-ifelse(medications_file_ageproper$CVD_CHD_lt10==1 & (
  medications_file_ageproper$CVD_CHD_Angina_lt10==0 &
    medications_file_ageproper$CVD_CHD_Infarction_lt10==0), 1, 0)

medications_file_ageproper$CVD_Other_lt10<-
  ifelse(medications_file_ageproper$CVD_lt10==1 & 
           ( medications_file_ageproper$CVD_ValveProblem_lt10==0 &
               medications_file_ageproper$CVD_CHD_Angina_lt10==0 &
               medications_file_ageproper$CVD_CHD_Infarction_lt10==0 &
               medications_file_ageproper$CVD_CHD_Other_lt10==0 &
               medications_file_ageproper$CVD_HF_Precursor_lt10==0 &
               medications_file_ageproper$CVD_HF_Failure_lt10==0 &
               medications_file_ageproper$CVD_Hypertension_lt10==0 &
               medications_file_ageproper$CVD_Arrhythmia_lt10==0 &
               medications_file_ageproper$CVD_Dyslipidemia_lt10==0 &
               medications_file_ageproper$CVD_PVD_lt10==0), 1, 0)

# create 10-32 other categories
medications_file_ageproper$CVD_CHD_Other_gt10<-
  ifelse(medications_file_ageproper$CVD_CHD_gt10==1 & (
    medications_file_ageproper$CVD_CHD_Angina_gt10==0 &
      medications_file_ageproper$CVD_CHD_Infarction_gt10==0), 1, 0)

medications_file_ageproper$CVD_Other_gt10<-
  ifelse(medications_file_ageproper$CVD_gt10==1 & 
           ( medications_file_ageproper$CVD_ValveProblem_gt10==0 &
               medications_file_ageproper$CVD_CHD_Angina_gt10==0 &
               medications_file_ageproper$CVD_CHD_Infarction_gt10==0 &
               medications_file_ageproper$CVD_CHD_Other_gt10==0 &
               medications_file_ageproper$CVD_HF_Precursor_gt10==0 &
               medications_file_ageproper$CVD_HF_Failure_gt10==0 &
               medications_file_ageproper$CVD_Hypertension_gt10==0 &
               medications_file_ageproper$CVD_Arrhythmia_gt10==0 &
               medications_file_ageproper$CVD_Dyslipidemia_gt10==0 &
               medications_file_ageproper$CVD_PVD_gt10==0), 1, 0)
#save analytic dataset for later use
#save file for analyses
save(medications_file_ageproper, file = paste0(WD,"FullDSAgeProper.Feb2022.RData"))