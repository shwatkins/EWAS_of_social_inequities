############################################################################################################################################################
#
# 
# Use data from EWAS catalog to regroup phenotypes into a smaller number of categories
# 
# copied from Hannah Elliot's github:
# https://github.com/hannah-e/collapse_EWAS_catalog_phenotypes/blob/9b65be66399d0c1d2fd71c2003dbf58e4e5b62ff/functional_analysis_regroup_EWAS_catalogue_phenotypes.R
#
###########################################################################################################################################################

#load data and libraries
library(readr)

# read in EWAS catalog
# two options: 
# option 1: can read off the ewas cat website (caution: sometimes this doesn't work)
#EWAS_catalog <- as.data.frame(read_delim("http://ewascatalog.org/static//docs/ewascatalog-results.txt.gz", delim="\t"))
# option 2: download the ewas catalog and use the file as input
EWAS_catalog <- as.data.frame(read_delim("ewascatalog-results.txt.gz", delim="\t"))

dim(EWAS_catalog)

#### reclassify columns in the object "EWAS_catalog" #######################################################################################################
# the following code categorises the findings in the ewas catalog
# these are not fixed so they can be changed

EWAS_catalog$phenotype<-EWAS_catalog$StudyID

temp<-grepl("age|aging|fetal_vs_adult_liver", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"age"

temp<-grepl("tissue", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"tissue"

temp<-grepl("smok|cotinine", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"smoking"

temp<-grepl("alcohol", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"alcohol"

temp<-grepl("sex", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"sex"

temp<-grepl("ancestry|ethnicity", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"ancestry"

temp<-grepl("cancer|carcinoma|adenoma|melanoma", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"cancer"

temp<-grepl("rheumatoid|ulcerative_colitis|lupus|sjogrens|crohn|inflammatory_bowel_disease|atopy|graves|psoriasis|multiple_sclerosis|cow_milk_allergy|26292806|ige", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"autoimmune"
#26292806 study captures atopy

temp<-grepl("hiv|human_immunodeficiency_virus", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"infection"

temp<-grepl("blood_pressure|chronic_kidney_disease|atrial_fibrillation|ischaemic_stroke|hypertension|myocardial_infarction|coronary_heart_disease|obesity|statin_use|Battram-T_arterial_distensibility|Battram-T_common_carotid_intima-media_thickness|Battram-T_pulse_rate", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"cardiovascular"

temp<-grepl("diabetes|chronic_kidney_disease|obesity|hepatic_fat|statin_use|insulin|glucose|homair|resistin|hba1c|adiponectin|leptin|liver_fat|proinsulin", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"metabolic"

temp<-grepl("creactive|c-reactive_protein|crp", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"inflammation"

temp<-grepl("triglycerides|hdl|highdensity_lipoprotein|lipemia|cholesterol|lipoprotein|ldl|vldl| idl|phospholipids|lp_a|Battram-T_concentration_of_idl_particles|Battram-T_total_lipids_in_idl|Battram-T_total_phosphoglycerides", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"lipid_lipoprotein"

temp<-grepl("perinatal|birth_weight|birthweight|maternal_underweight|plasma_folate|prenatal|pregnancy|preterm_birth|season_of_birth|breastfeeding|31230546|33396735|utero|fetal_intolerance_of_labor|Starling-PS_maternal_serum|gestational_weight_gain|parity|fetal_brain_development", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"perinatal"
#31230546 study captures hypertension in pregnancy, preeclampsia
#33396735 study captures bottle, breast and mixed feeding behaviours

temp<-grepl("birth_weight|birthweight|preterm_birth|31230546|fetal_intolerance_of_labor|gestational_weight_gain|fetal_brain_development", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"prenatal_complications"
#31230546 study captures hypertension in pregnancy, preeclampsia
#33396735 study captures bottle, breast and mixed feeding behaviours
temp<-grepl("perinatal|maternal_underweight|plasma_folate|prenatal|pregnancy|preterm_birth|season_of_birth|breastfeeding|33396735|utero|Starling-PS_maternal_serum|gestational_weight_gain|parity", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"prenatal_exposures"
#31230546 study captures hypertension in pregnancy, preeclampsia
#33396735 study captures bottle, breast and mixed feeding behaviours

temp<-grepl("copd|fev1|fvc|chronic_obstructive_pulmonary_disease|asthma|lung_function|28621160|26681806", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"lung"
#28621160 study captures fev1 and fvc
#26681806 study captures cardiovascular biomarker GDF-15

temp<-grepl("bmi|body_mass_index|waist_circumference|arm_circumference|head_circumference|hip_circumference|fat|mass|weight|height|skinfold|waist|bone_mineral_density|Battram-T_head|Battram-T_hip|Battram-T_leg|Battram-T_pelvis|Battram-T_ribs|Battram-T_spine|Battram-T_total_body|Battram-T_trunk|Battram-T_alpha_neck_angle__hip_measurement|Battram-T_theta_neck_angle__hip_measurement|Battram-T_arm_area|Battram-T_arm_bone_mineral_content", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"anthropometric"

temp<-grepl("dementia|palsy|alzheimer|amyloid_plaques|cognitive|infant_attention|cortical|neurobehavioural_scale|seizures|parkinson|social_communication_deficits|hippocampus_volume|thalamus_volume|apolipoprotein|apoe", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"neurological"

temp<-grepl("ptsd|schizophrenia|depressive_disorder|depressive_symptoms|attention_deficit_hyperactivity_disorder|wellbeing|depression|personality_disorder|tic_disorders|aggressive_behaviour|stress|anxiety|conduct_problems|antidepressant_use|response_to_antidepressants|psychosis|psychotic_experiences", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"mental__health"


temp<-grepl("socioeconomic_position|maternal_education|educational_attainment|31062658", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"SEP_and_education"
#31062658 captures various markers of maternal SEP

temp<-grepl("child_abuse|exposure_to_community_and_domestic_violence|exposure_to_intimate_partner_violence|bullying|perceived_discrimination|victimization|30905381|abuse", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"violence_adversity"

temp<-grepl("healthy_eating|mediterranean_diet|diet_quality|pufa_intake|tea_consumption|fruit_consumption|juice_consumption|30101351|folate_intake|serum_copper", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"diet_"
#30101351 study captures one-carbon metabolism nutritents

temp<-grepl("selenium|cadmium|arsenic|ppdde|substance_use|illicit_drug_use|esterogen_exposure|altitude|noise_pollution|noise_polution|sunlight_duration|Battram-T_18_2_linoleic_acid|docosahexaenoic_acid|lead", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"environment"

temp<-grepl("atmospheric_iron|atmospheric_nickel|atmospheric_vanadium|nitrogen_dioxide_exposure|proximity_to_major_roadways|particulate_matter|air_pollution|nitrogen", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"air_pollutants"

temp<-grepl("24014485|33413638|trimethylaminenoxide", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"metabolites"
#24014485 study captures study of blood serum metabolites
#33413638 NMR of predominantly lipids

temp<-grepl("31282290", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==T]<-"miRNA"
#31282290 study captures EWAS of miRNA expression

#classify everything else as "other":
temp<-grepl("age|tissue|smoking|alcohol|sex|ancestry|cancer|autoimmune|infection|cardiovascular|metabolic|inflammation|prenatal_complications|prenatal_exposures|lung|lipid_lipoprotein|anthropometric|neurological|mental__health|SEP_and_education|violence_adversity|diet_|environment|air_pollutants|metabolites|miRNA", EWAS_catalog$phenotype)
EWAS_catalog$phenotype[temp==F]<-"other"

#### check classifications ################################################################################################################################
dim(EWAS_catalog)
table(EWAS_catalog$phenotype)

#### Create table of unique CpGs per category to use in enrichment analysis ###############################################################################

#cut to required p-value threshold
EWAS_catalog<-EWAS_catalog[EWAS_catalog$P<2.4e-07,] #EDIT TO REDUCE TO REQUIRED P-THRESHOLD
dim(EWAS_catalog)
#create output table
EWAS_catalog_collapsed_freq_table<-table(EWAS_catalog$phenotype)
dim(EWAS_catalog_collapsed_freq_table)

#get number of unique CpGs in each category
for (i in names(EWAS_catalog_collapsed_freq_table)){
  
  temp<-EWAS_catalog[EWAS_catalog$phenotype==i,]
  temp<-length(temp$CpG[duplicated(temp$CpG)==F])
  EWAS_catalog_collapsed_freq_table[names(EWAS_catalog_collapsed_freq_table)==i]<-temp
}

rm(temp)

dim(EWAS_catalog_collapsed_freq_table)
EWAS_catalog_collapsed_freq_table

save(EWAS_catalog,EWAS_catalog_collapsed_freq_table, file="ewas_catalog_for_enrichment.Rdata")

#### end ################################################################################################################################


