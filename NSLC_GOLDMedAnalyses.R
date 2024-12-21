# GOLD medication-adjusted analyses
# output a text file with medication-adjusted analyses

WD<-"./"
setwd(WD)
source(paste0(WD,"NSLC/ProcessConditionalLogitApr2021.R"))

# load dataset containing table 'medications_file_ageproper'
load( file = paste0(WD,"FullDSAgeProper.Feb2022.RData"))

library(haven)
library(reshape2)
library(survival)
library(epitools)

#create hash map of variable names and human readable labels
library(r2r)
HumanLabels<-hashmap()

#includes non-topical medications that are used to treat conditions.
HumanLabels[
  c("InfectInflam_COPD_lt10","steroidsoral","macrolides","InfectInflam_GI_GERD_lt10","proton","InfectInflam_GI_Gastr_lt10",
              "h2rec","InfectInflam_AI_Psor_lt10","methotrexate","antacid","allimmune_lt10",
              "InfectInflam_AI_Lupus_lt10","laba","saba","nsaids","corticoster","antimusc","InfectInflam_AI_Rose_lt10","tetracyclines")
  ] <- 
  c("COPD", "Oral Corticosteroids","Macrolides","GERD","PPIs","Gastritis/NIGGC","H2 Receptor Agonists","Psoriasis",
    "Methotrexate","Antacid","Immunosupressing","Lupus","LABA","SABA","NSAIDs","Inhaled Corticosteroids","Antimuscarinic Bronchodilators", "Rosacea","Tetracyclines")


# perform all medication-adjusted analyses
GOLDMedicationAdjusted<-data.frame(Condition=character(),Medication=character(),ConditionReadable=character(), MedicationReadble=character(),
                                   OR=double(),low=double(),high=double(),ORfmt=character(),pvalue=double(),
                                   ORadj=double(),lowadj=double(),highadj=double(),OR1fmt=character(),padj=double(),
                                   ORmed=double(),lowmed=double(),highmed=double(),OR2fmt=character(),pmed=double())

a1<-EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_COPD_lt10","steroidsoral",HumanLabels)

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,a1)

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_COPD_lt10","laba",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_COPD_lt10","saba",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_COPD_lt10","corticoster",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_COPD_lt10","antimusc",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_COPD_lt10","macrolides",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_GI_GERD_lt10","proton",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_GI_GERD_lt10","h2rec",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_GI_GERD_lt10","antacid",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_GI_Gastr_lt10","proton",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_GI_Gastr_lt10","h2rec",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_GI_Gastr_lt10","antacid",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_AI_Psor_lt10","methotrexate",HumanLabels))

medications_file_ageproper$allimmune_lt10<-
  medications_file_ageproper$antiproimmuno|medications_file_ageproper$immunes|medications_file_ageproper$immuneresponse


GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_AI_Lupus_lt10","allimmune_lt10",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_AI_Lupus_lt10","steroidsoral",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_AI_Lupus_lt10","nsaids",HumanLabels))


GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_AI_Rose_lt10","macrolides",HumanLabels))

GOLDMedicationAdjusted<-rbind(GOLDMedicationAdjusted,
                              EstimateAdjAssociation(medications_file_ageproper,"CASE_STATUS","indexage", "SETID","InfectInflam_AI_Rose_lt10","tetracyclines",HumanLabels))


colnames(GOLDMedicationAdjusted)<-c("Condition", "Medication","Condition Label","Medication Lable","ORunadj", "ORlow","ORhigh", "ORFmt","P",
                                  "ORadj", "ORadjLow", "ORadjHigh", "ORCondFmt","Padj",
                                 "ORmed", "ORmedlow","ORmedhigh","ORMedFmt","Pmed")

GOLDMedicationAdjusted <- apply(GOLDMedicationAdjusted,2,as.character)

write.table(as.data.frame(GOLDMedicationAdjusted),paste0(WD,"GoldMedAdjustedRawAug152024.txt"),row.names = F,sep="\t")

