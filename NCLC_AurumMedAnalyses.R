# Aurum medication-adjusted analyses
# output a text file with medication-adjusted analyses

WD<-"./"
setwd(WD)
source(paste0(WD,"ProcessConditionalLogitApr2021.R"))

# load and read in primay DS
load( file = paste0(WD,"AurumDeDupped.May2024.RData"))

###################################################################################
#create hash map of variable names and human readable labels
library(r2r)
HumanLabels<-hashmap()

# no need to look at lupus-medication adjusted association since it was insignificant or macrolides, inhaled corticosteroids
HumanLabels[c("copd_lt10", "steroidsoral_lt10", "antimuscbronc_lt10","SABA_lt10","LABA_lt10",
              "gi_gerd_lt10", "gi_gastrexact_lt10","antacid_lt10","PPI_lt10", "H2recanta_lt10",
              "methotrexate_lt10", "ii_ai_psor_lt10")]<-
  c("COPD", "Oral Corticosteroids","Antimuscarinic Bronchodilators","SABA","LABA",
    "GERD","Gastritis/NIGGC","Antacid","PPIs","H2-receptor Antagonists",
    "Methotrexate","Psoriasis")

AurumMedicationAdjusted<-data.frame(Condition=character(),Medication=character(),ConditionReadable=character(), MedicationReadble=character(),
                                               OR=double(),low=double(),high=double(),ORfmt=character(),pvalue=double(),
                                               ORadj=double(),lowadj=double(),highadj=double(),OR1fmt=character(),padj=double(),
                                               ORmed=double(),lowmed=double(),highmed=double(),OR2fmt=character(),pmed=double())
                                    
a1<-EstimateAdjAssociation(NonDupAurumFullMatches,"case_status","indexage", "SETID","copd_lt10","steroidsoral_lt10",HumanLabels)

AurumMedicationAdjusted<-rbind(AurumMedicationAdjusted,a1)


AurumMedicationAdjusted<-rbind(AurumMedicationAdjusted,
                               EstimateAdjAssociation(NonDupAurumFullMatches,"case_status","indexage", "SETID","copd_lt10","antimuscbronc_lt10",HumanLabels))

AurumMedicationAdjusted<-rbind(AurumMedicationAdjusted,
                               EstimateAdjAssociation(NonDupAurumFullMatches,"case_status","indexage", "SETID","copd_lt10","SABA_lt10",HumanLabels))

AurumMedicationAdjusted<-rbind(AurumMedicationAdjusted,
                               EstimateAdjAssociation(NonDupAurumFullMatches,"case_status","indexage", "SETID","copd_lt10","LABA_lt10",HumanLabels))

AurumMedicationAdjusted<-rbind(AurumMedicationAdjusted,
                               EstimateAdjAssociation(NonDupAurumFullMatches,"case_status","indexage", "SETID","gi_gerd_lt10","antacid_lt10",HumanLabels))

AurumMedicationAdjusted<-rbind(AurumMedicationAdjusted,
                              EstimateAdjAssociation(NonDupAurumFullMatches,"case_status","indexage", "SETID","gi_gerd_lt10","PPI_lt10",HumanLabels))


AurumMedicationAdjusted<-rbind(AurumMedicationAdjusted,
                              EstimateAdjAssociation(NonDupAurumFullMatches,"case_status","indexage", "SETID","gi_gerd_lt10","H2recanta_lt10",HumanLabels))


AurumMedicationAdjusted<-rbind(AurumMedicationAdjusted,
                               EstimateAdjAssociation(NonDupAurumFullMatches,"case_status","indexage", "SETID","gi_gastrexact_lt10","antacid_lt10",HumanLabels))

AurumMedicationAdjusted<-rbind(AurumMedicationAdjusted,
                              EstimateAdjAssociation(NonDupAurumFullMatches,"case_status","indexage", "SETID","gi_gastrexact_lt10","PPI_lt10",HumanLabels))

AurumMedicationAdjusted<-rbind(AurumMedicationAdjusted,
                               EstimateAdjAssociation(NonDupAurumFullMatches,"case_status","indexage", "SETID","gi_gastrexact_lt10","H2recanta_lt10",HumanLabels))

AurumMedicationAdjusted<-rbind(AurumMedicationAdjusted,
                              EstimateAdjAssociation(NonDupAurumFullMatches,"case_status","indexage", "SETID","ii_ai_psor_lt10","methotrexate_lt10",HumanLabels))

# no need to adjust for Lupus or rosacea because not significant.
colnames(AurumMedicationAdjusted)<-c("Condition", "Medication","Condition Label","Medication Lable","ORunadj", "ORlow","ORhigh", "ORFmt","P",
                                                                       "ORadj", "ORadjLow", "ORadjHigh", "ORCondFmt","Padj",
                                                                       "ORmed", "ORmedlow","ORmedhigh","ORMedFmt","Pmed")

AurumMedicationAdjusted <- apply(AurumMedicationAdjusted,2,as.character)

write.table(as.data.frame(AurumMedicationAdjusted),paste0(WD,"AurumMedAdjustedRawAug62024.txt"),row.names = F,sep="\t")

