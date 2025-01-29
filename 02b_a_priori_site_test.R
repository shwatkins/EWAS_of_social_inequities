## a-priori sites
## These are sites that have previously been associated with our exposures
## Not sure if we will 1) 

library(readr)
library(meffil)
library(data.table)

EWAS_catalog <- as.data.frame(read_delim("http://ewascatalog.org/static//docs/ewascatalog-results.txt.gz", delim="\t"))
dim(EWAS_catalog)

EWAS_catalog$phenotype<-EWAS_catalog$StudyID

#######
print("starting SES")
## Nearest equivalent to household poverty:income ratio ##----

# https://pubmed.ncbi.nlm.nih.gov/26295359/ # Needham 2015
## AVP, BDNF, FKBP5, and OXTR, CCL1, CD1D, and F8
needham <- c("AVP", "BDNF", "FKBP5", "OXTR", "CCL1", "CD1D", "F8")

# https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-021-01189-0 
# Cerruti 2021 - review

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6462839/ # Bush 2018 
## all results are FDR adjusted with no p-reporting so hard to threshold
## additionally there are some very large Ns of sites identified;
## no adjustment for cell counts or sex so perhaps that's why.

# https://www.annualreviews.org/doi/pdf/10.1146/annurev-publhealth-052020-105613 # Martin 2022 - review

# https://academic.oup.com/hmg/advance-article/doi/10.1093/hmg/ddac171/6650946 # Liu 2022 - job loss, income reduction, low family income, financial hardship, major financial problemscg10498926
liu <- unique(c("cg10639395","cg14613617","cg18913843","cg26112574","cg00043076","cg03764134","cg08641963",
                "cg21112452","cg06957310","cg04695077","cg09179691","cg05767720","cg12715235","cg19621008",
                "cg26114043","cg11308211","cg21609106","cg07031083","cg23685969","cg12608703","cg09087363",
                "cg11256802","cg03305840","cg14023073","cg02345321","cg00863271","cg27513979","cg19543471",
                "cg11883141","cg19257274","cg27006129","cg17100528","cg00611614","cg01022012","cg26360968",
                "cg24121967","cg05060427","cg19498110","cg24938210","cg19260606","cg13671374","cg10808198",
                "cg22275306","cg04913057"))

# https://pubmed.ncbi.nlm.nih.gov/30771258/ # McDade 2022
# results presented by gene not cpg, and high numbers 
# of associations with SES (about 2.5k)
# they don't provide a full list and just include the top 20 genes 
# in the main paper:
# CD44,ZNF827,MAD1L1,UBE4A,NLRC5,SPARC,SFRS8,EZH2,MEFV,CHST15,PAQR5,OTUD6B,NDRG1,WDFY1,EBF4,SLC7A7,DTX3L,PARP9,SLC24A4,VTI1A
mcdade <- c("CD44","ZNF827","MAD1L1","UBE4A","NLRC5","SPARC","SFRS8",
            "EZH2","MEFV","CHST15","PAQR5","OTUD6B","NDRG1","WDFY1",
            "EBF4","SLC7A7","DTX3L","PARP9","SLC24A4","VTI1A")

# https://academic.oup.com/ije/article/44/4/1320/672172 # Stringhini 2015
# candidate gene analysis using genes selected from 3 papers "exploring SES and gene regulation patterns 
# (in humans and primates), on the basis of their involvement in SES-related inflammation".
# Associations survived FDR threshold:
# household's occupational position EWAS:
# NFATC1, CXCL2, PTGS2, MAP2K5, MAP3K6, IL1A, GPR132, TNFRSF11A, ADM, and OLR1
# SES trajectories EWAS:
# NFATC1, MAP3K6, IL1A, GPR132, CXCL2, MAP2K5
stringhini <- c("NFATC1", "CXCL2", "PTGS2", "MAP2K5", "MAP3K6",
                "IL1A", "GPR132", "TNFRSF11A", "ADM", "OLR1", 
                "NFATC1", "MAP3K6", "IL1A", "GPR132", "CXCL2", "MAP2K5")

# https://srcd.onlinelibrary.wiley.com/doi/full/10.1111/cdev.12486 # Beach 2016
# # EWAS results not reported in the paper

temp<-grepl("socioeconomic|socio-economic", EWAS_catalog$phenotype) 
# Laubach 2019, Loucks 2016, Santos Jr 2019
socioeconomic <- EWAS_catalog[temp==T,]
dim(socioeconomic)
length(unique(socioeconomic$phenotype))
table(socioeconomic$phenotype)
# remove reuben study because this looks at area level
socioeconomic <- socioeconomic[!grepl("32478847", socioeconomic$phenotype),]
# reduce to selected p-value threshold
socioeconomic <- socioeconomic[socioeconomic$P<2.4e-7,]

# get genes
genes <- unique(c(needham,mcdade,stringhini))
# only select genes common to the 3 studies
print(length(genes))
## identify all CpG sites annotated to these genes
features <- rbind(meffil.featureset("epic"),
                  meffil.featureset("450k"))
features <- features[match(unique(features$name), features$name),]
features$genes <- strsplit(features$gene.symbol, ";")
feature.idx <- sapply(genes, function(gene) grep(gene, features$gene.symbol), simplify=F)
feature.idx <- sapply(genes, function(gene) {
  idx <- feature.idx[[gene]]
  idx[which(sapply(features$genes[idx], function(genes) gene %in% genes))]
}, simplify=F)

## identify all CpG sites within 100kb of these genes
## http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz
##  downloaded Oct 9, 2020
ucsc <- read.table("/ncbiRefSeq.txt.gz", stringsAsFactors=F,comment.char="",header=F,sep="\t")
colnames(ucsc) <- c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
ucsc.idx <- sapply(genes, function(gene) which(ucsc$name2 == gene), simplify=F)
ucsc <- ucsc[unlist(ucsc.idx), c("chrom","txStart","txEnd","name2")]
colnames(ucsc) <- c("chr","start","end","gene")

ucsc$feature.idx <- sapply(1:nrow(ucsc), function(i) {
  which(features$chr == ucsc$chr[i]
        & features$pos >= ucsc$start[i]-1000
        & features$pos <= ucsc$end[i]+1000)
})
## number of CpG sites within 100kb but not annotated to the genes
length(setdiff(unlist(ucsc$feature.idx), unlist(feature.idx)))
## 1388
## the within 100kb list includes all sites annotated to the genes
length(setdiff(unlist(feature.idx), unlist(ucsc$feature.idx)))
## 0
## the number of sites within 100kb of a gene
length(unique(unlist(ucsc$feature.idx)))
## 1594
gene.sites <- features$name[unique(c(unlist(feature.idx), unlist(ucsc$feature.idx)))]

print("socioeconomic ewas cat n_occur:")
n_occur <- data.frame(table(socioeconomic$CpG))
n_occur <- n_occur[n_occur$Freq > 1,]
test <- c(liu,socioeconomic$CpG,gene.sites)
n_occur <- data.frame(table(test))
n_occur <- n_occur[n_occur$Freq > 1,]
print("socioeconomic combined n_occur:")
ses_freq <- as.character(n_occur$test)

ses_cpgs <- Reduce(intersect, list(gene.sites,liu,socioeconomic$CpG))
print(length(ses_cpgs))
print("length SES cpgs:")
length(ses_cpgs)

########
print("starting neighbourhood")
## Neighborhood disadvantage (closest equivalent to ICE) ##----

# https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2766579 # Reuben 2020
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5687339/ # Smith 2017
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6430282/ # Giurgescu 2019 - review
# https://academic.oup.com/hmg/advance-article/doi/10.1093/hmg/ddac171/6650946 # Liu 2022 - neighborhood disadvantage
liu <- unique(c("cg13075263","cg14062745","cg03912614","cg05465572","cg13294509","cg04007726",
                "cg20102336","cg08638097","cg12651540","cg27025925","cg00450784","cg11632238",
                "cg24526499","cg20559403","cg23405172","cg14212190","cg03192060"))

temp<-grepl("disadvantage|residence|Neighborhood|neighborhood", EWAS_catalog$phenotype) # Dunn 2019
area_socio_measure <- EWAS_catalog[temp==T,]
# reduce to selected p-value threshold
area_socio_measure <- area_socio_measure[area_socio_measure$P<2.4e-7,]

print("area_socio_measure ewas cat n_occur:")
n_occur <- data.frame(table(area_socio_measure$CpG))
n_occur <- n_occur[n_occur$Freq > 1,]
test <- c(liu,area_socio_measure$CpG)
n_occur <- data.frame(table(test))
n_occur <- n_occur[n_occur$Freq > 1,]
neighborhood_cpgs <- n_occur$test
print("length neighborhood cpgs:")
print(length(neighborhood_cpgs))

############
print("starting air pollution")
## air pollution ##----
# https://www.nature.com/articles/s41598-021-83577-3 # Prunicki 2021
# candidate gene analysis but candidates are based on old literature
# use mass cytometry by time of flight (CyTOF) to more sensitively measure 30-40 cell markers simultaneously,

# https://www.tandfonline.com/doi/full/10.1080/15592294.2021.1900028 # Chi 2021
# DMRs using bumphunter
# also some EWAS assocs
chi_2021_pm2.5 <- c("cg05926640","cg04310517")
chi_2021_nox <- c("cg11756214")


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5131503/ # Chi 2016
# candidate site analysis using CpGs previously associated with gene expression
# in MESA. 2,713 expression-associated methylation sites (eMS) 
# "adjusted for age, sex, race/ethnicity (black, Hispanic, white), household 
# income, education, neighborhood socioeconomic status, smoking (smoking status 
# and pack-years), secondhand smoke, body mass index, recent infection, methyl 
# nutrient intake (continuous folate, vitamin B12, vitamin B6, methionine, zinc),
# physical activity, study site, microarray chip, and chip position (DNA 
# methylation analysis only)
chi_pm2.5 <- c("cg20455854","cg07855639","cg07598385","cg17360854","cg23599683")
# but no associations with NOx

# https://pubmed.ncbi.nlm.nih.gov/29597345/ # Plusquin 2018
# PM10 assoc longitudinally with 1 CpG in children
plusquin <- c("cg21785536")

# https://www.sciencedirect.com/science/article/pii/S0160412018308882?via%3Dihub # Li 2018
# randomised trial of air purifier to reduce PM2.5
# have taken associations from table 4 that pass GWS (9e-08)
li_2018 <- read.csv("/user/work/sw14321/from_bc4/NIH_dnam_inequities/Li_2018_airpollution_cpgs.csv",stringsAsFactors = F)
li_2018 <- li_2018[li_2018$pval<9e-8,]
li_2018 <- li_2018$Probe.ID

temp<-grepl("air_pollution|nitrogen_dioxide_exposure|proximity_to_major_roadways|particulate_matter|air_pollution|nitrogen|pm2.5", EWAS_catalog$phenotype)
air_pollution <- EWAS_catalog[temp==T,]
dim(air_pollution)
length(unique(air_pollution$phenotype))
table(air_pollution$phenotype)
# 26731791 - Panni 2016 - particles
# 27058926 - Kingsley 2016 - maternal proximity to roadways
# 29410382 - Lichtenfels 2018 particulate; NO2
# 31148503 - Gruzieva 2019 particulate_matter pm2.5;pm10
# 31208937 - Gondalia 2019 particulate_matter pm2.5;pm10
# 32484729 - Eze 2020 NO2;particulate
# 33818294 - Chi 2022 (MESA) pm2.5
# reduce to selected p-value threshold
air_pollution <- air_pollution[air_pollution$P<2.4e-7,]
table(air_pollution$phenotype)
# remove Lichtenfels as there's something wrong with the p values
air_pollution <- air_pollution[!air_pollution$phenotype=="29410382_Lichtenfels-AJ_nitrogen_dioxide_exposure",]
Lichtenfels_nox <- c("cg04908668","cg14938677","cg00344801","cg18379295","cg25769469","cg02234653","cg08500171")

# Jan 2025 resubmission addition: two papers

koenigsberg_nox <- c("cg24537688", "cg24269657","cg07635198","cg08097847",
                     "cg01593570","cg13538431","cg23717809","cg06738602",
                     "cg00212245","cg01906102","cg04902542","cg05411829",
                     "cg08052751","cg11514293","cg10011083","cg16722016",
                     "cg18580650","cg11130461","cg24544356","cg17608585")

Holliday_nox <- c("cg01885635")

# taken from table S2
Wang_2022_pm2.5 <- read.csv("/user/work/sw14321/from_bc4/NIH_dnam_inequities/data/Wang_2022_pm25.csv",header = F)
Wang_2022_pm2.5 <- as.character(Wang_2022_pm2.5[,1])
# taken from table S3
Wang_2022_blackcarbon <- read.csv("/user/work/sw14321/from_bc4/NIH_dnam_inequities/data/Wang_2022_black_carbon.csv",header = F)
Wang_2022_blackcarbon <- as.character(Wang_2022_blackcarbon[,1])
# from 2023 paper: NO3 and organic carbon are the most simmilar. Taken from table s2.
Wang_2023 <- c("cg08045932","cg27637895","cg25551691","cg23010344","cg07376297",
               "cg24697733","cg03295933","cg24229819","cg21271452","cg03070236",
               "cg24548108","cg06754557","cg01962727","cg06462347","cg00207921",
               "cg01995099","cg07691004")

print("air_pollution ewas cat n_occur:")
n_occur <- data.frame(table(air_pollution$CpG))
n_occur <- n_occur[n_occur$Freq > 1,]
test <- c(chi_pm2.5,plusquin,li_2018,air_pollution$CpG,Lichtenfels_nox,chi_2021_pm2.5,chi_2021_nox,
          koenigsberg_nox,Holliday_nox,Wang_2022_pm2.5,Wang_2022_blackcarbon,Wang_2023)
nox <- c(chi_2021_nox,Lichtenfels_nox,koenigsberg_nox,Holliday_nox)
n_occur <- data.frame(table(test))
n_occur <- n_occur[n_occur$Freq > 1,]
airpollution_cpgs <- as.character(n_occur$test)

################
print("starting eduction")
## Education ##----

temp<-grepl("maternal_education|educational_attainment|31062658|education", EWAS_catalog$phenotype)
education <- EWAS_catalog[temp==T,]
education <- education[education$P<2.4e-7,]

print("education ewas cat n_occur:")
n_occur <- data.frame(table(education$CpG))
n_occur <- n_occur[n_occur$Freq > 1,]
print(n_occur)
education_cpgs <- as.character(n_occur$Var1)
print("length education cpgs:")
print(length(education_cpgs))

#################

## discrimination/racial discrimination ##----

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5741522/ # Barcelona de Mendoza 2018
# 9 sites pass FDR correction for MLD scale only:
mld <- c("cg09437479","cg16143945","cg03317714","cg19958721","cg05711042","cg02578462","cg07726178",
         "cg15209133","cg04928761")
# alternatively only two pass the geneome wide threshold:
mld_gws <- c("cg09437479","cg16143945")

# https://pubmed.ncbi.nlm.nih.gov/30144780/ # Santos 2018
# candidate gene analysis of EDS. Genes pyrosequenced.
# 3 candidate genes associated with stree were selected
santos_2018 <- c("NR3C1","BDNF")

temp<-grepl("discrim", EWAS_catalog$phenotype) # van der Laan 2020
discrimination <- EWAS_catalog[temp==T,]
discrimination <- discrimination[discrimination$P<2.4e-7,]

# get genes
genes <- santos_2018
## identify all CpG sites annotated to these genes
features <- rbind(meffil.featureset("epic"),
                  meffil.featureset("450k"))
features <- features[match(unique(features$name), features$name),]
features$genes <- strsplit(features$gene.symbol, ";")
feature.idx <- sapply(genes, function(gene) grep(gene, features$gene.symbol), simplify=F)
feature.idx <- sapply(genes, function(gene) {
  idx <- feature.idx[[gene]]
  idx[which(sapply(features$genes[idx], function(genes) gene %in% genes))]
}, simplify=F)

## identify all CpG sites within 100kb of these genes
## http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz
##  downloaded Oct 9, 2020
ucsc <- read.table("/ncbiRefSeq.txt.gz", stringsAsFactors=F,comment.char="",header=F,sep="\t")
colnames(ucsc) <- c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
ucsc.idx <- sapply(genes, function(gene) which(ucsc$name2 == gene), simplify=F)
ucsc <- ucsc[unlist(ucsc.idx), c("chrom","txStart","txEnd","name2")]
colnames(ucsc) <- c("chr","start","end","gene")

ucsc$feature.idx <- sapply(1:nrow(ucsc), function(i) {
  which(features$chr == ucsc$chr[i]
        & features$pos >= ucsc$start[i]-1000
        & features$pos <= ucsc$end[i]+1000)
})
## number of CpG sites within 100kb but not annotated to the genes
length(setdiff(unlist(ucsc$feature.idx), unlist(feature.idx)))
## 1388
## the within 100kb list includes all sites annotated to the genes
length(setdiff(unlist(feature.idx), unlist(ucsc$feature.idx)))
## 0
## the number of sites within 100kb of a gene
length(unique(unlist(ucsc$feature.idx)))
## 1594
gene.sites <- features$name[unique(c(unlist(feature.idx), unlist(ucsc$feature.idx)))]
#writeLines(gene.sites, con="gene-sites.txt")

print("discrimination ewas cat n_occur:")
n_occur <- data.frame(table(discrimination$CpG))
n_occur <- n_occur[n_occur$Freq > 1,]
test <- c(mld,discrimination$CpG,gene.sites)
n_occur <- data.frame(table(test))
n_occur <- n_occur[n_occur$Freq > 1,]
discrimination_cpgs <- unique(c(mld,discrimination$CpG,gene.sites))

print("length discrimination cpgs:")
print(length(discrimination_cpgs))

#########
print("starting stress")
## stress ##----

temp<-grepl("stress", EWAS_catalog$phenotype)
stress <- EWAS_catalog[temp==T,]
stress <- stress[stress$P<2.4e-7,]

stress_cpgs <- unique(stress$CpG)

print("stress ewas cat n_occur:")
n_occur <- data.frame(table(stress$CpG))
n_occur <- n_occur[n_occur$Freq > 1,]
stress_cpgs <- as.character(n_occur$Var1)

print("length stress cpgs:")
print(length(stress_cpgs))

#############

## CRP ##----

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9065016/ # Weilscher 2022
tt.crp=read.csv("/Weilscher_crp_cpgs.csv",header = T)
weilscher_cpgs <- as.character(tt.crp$ID)

temp<-grepl("c-reactive|CRP|crp", EWAS_catalog$phenotype)
crp <- EWAS_catalog[temp==T,]
crp <- crp[crp$P<2.4e-7,]

test <- c(weilscher_cpgs,crp$CpG)
n_occur <- data.frame(table(test))
n_occur <- n_occur[n_occur$Freq > 1,]
crp_cpgs <- as.character(n_occur$test)
print("length crp cpgs:")
print(length(crp_cpgs))

#################

## mental health ##----

temp<-grepl("ptsd|schizophrenia|depressive_disorder|depressive_symptoms|attention_deficit_hyperactivity_disorder|wellbeing|depression|personality_disorder|tic_disorders|aggressive_behaviour|stress|anxiety|conduct_problems|antidepressant_use|response_to_antidepressants|psychosis|psychotic_experiences", EWAS_catalog$phenotype)
mental_health <- EWAS_catalog[temp==T,]
mental_health <- mental_health[mental_health$P<2.4e-7,]
mentalhealth_cpgs <- unique(mental_health$CpG)

print("mental_health ewas cat n_occur:")
n_occur <- data.frame(table(mental_health$CpG))
n_occur <- n_occur[n_occur$Freq > 1,]
n_occur <- as.data.frame(n_occur)
n_occur <- n_occur[order(n_occur$Freq),]
print(n_occur[1:10,])
mentalhealth_cpgs <- as.character(n_occur$Var1)

print("length mental health cpgs:")
print(length(mentalhealth_cpgs))

##########################

###########################


# now check for these sites in each EWAS

check_a_priori_ewas <- function(x){
  print("start of function, data input:")
  temp.ewas <- x
  temp.ewas <- temp.ewas$analyses$all$table
  temp.ewas$name <- as.character(rownames(temp.ewas))
  temp.ewas$ewas_set <- "NA"
  # SES
  ewas.ses <- temp.ewas[temp.ewas$name %in% ses_cpgs,]
  if(nrow(ewas.ses>0)){
    ewas.ses$ewas_set <- "SES"
  }
  ewas.ses <- ewas.ses[ewas.ses$p.value<0.05/length(ses_cpgs),]
  # neighborhood
  ewas.neighborhood <- temp.ewas[temp.ewas$name %in% neighborhood_cpgs,]
  if(nrow(ewas.neighborhood>0)){
    ewas.neighborhood$ewas_set <- "Neighborhood"
  }
  ewas.neighborhood <- ewas.neighborhood[ewas.neighborhood$p.value<0.05/length(neighborhood_cpgs),]
  # air pollution
  ewas.pollution <- temp.ewas[temp.ewas$name %in% airpollution_cpgs,]
  if(nrow(ewas.pollution>0)){
    ewas.pollution$ewas_set <- "Pollution"
  }
  ewas.pollution <- ewas.pollution[ewas.pollution$p.value<0.05/length(airpollution_cpgs),]
  # NOx
  ewas.nox <- temp.ewas[temp.ewas$name %in% nox,]
  if(nrow(ewas.nox>0)){
    ewas.nox$ewas_set <- "Pollution"
  }
  ewas.nox <- ewas.nox[ewas.nox$p.value<0.05/length(ewas.nox),]
  
  # education
  ewas.education <- temp.ewas[temp.ewas$name %in% education_cpgs,]
  if(nrow(ewas.education>0)){
    ewas.education$ewas_set <- "Education"
  }
  ewas.education <- ewas.education[ewas.education$p.value<0.05/length(education_cpgs),]
  # discrimination
  ewas.discrimination <- temp.ewas[temp.ewas$name %in% discrimination_cpgs,]
  if(nrow(ewas.discrimination>0)){
    ewas.discrimination$ewas_set <- "Discrimination"
  }
  ewas.discrimination <- ewas.discrimination[ewas.discrimination$p.value<0.05/length(discrimination_cpgs),]
  # stress
  ewas.stress <- temp.ewas[temp.ewas$name %in% stress_cpgs,]
  if(nrow(ewas.stress>0)){
    ewas.stress$ewas_set <- "Stress"
  }
  ewas.stress <- ewas.stress[ewas.stress$p.value<0.05/length(stress_cpgs),]
  # crp
  ewas.crp <- temp.ewas[temp.ewas$name %in% crp_cpgs,]
  if(nrow(ewas.crp>0)){
    ewas.crp$ewas_set <- "CRP"
  }
  ewas.crp <- ewas.crp[order(ewas.crp$p.value),]
  ewas.crp <- ewas.crp[ewas.crp$p.value<0.05/length(crp_cpgs),]
  # mental health
  ewas.mentalhealth <- temp.ewas[temp.ewas$name %in% mentalhealth_cpgs,]
  if(nrow(ewas.mentalhealth>0)){
    ewas.mentalhealth$ewas_set <- "Mentalhealth"
  }
  ewas.mentalhealth <- ewas.mentalhealth[order(ewas.mentalhealth$p.value),]
  ewas.mentalhealth <- ewas.mentalhealth[ewas.mentalhealth$p.value<0.05/length(mentalhealth_cpgs),]
  results_list <- list(ewas.ses,ewas.neighborhood,ewas.pollution,ewas.nox,ewas.education,ewas.discrimination,ewas.stress,ewas.crp,ewas.mentalhealth)
  names(results_list)<-c("ewas.ses","ewas.neighborhood","ewas.pollution","ewas.nox","ewas.education","ewas.discrimination","ewas.stress","ewas.crp","ewas.mentalhealth")
  return(results_list)
}


####

## MESA 

raceethnicity <- c("whiteNH","blackNH","hispanic")
exposure <- c("jim_crow_birth_state",
              "parent_educ_low","parent_educ_mid",
              "educ_low","educ_mid",
              "hh_income_pov_ratio",
              "ICEraceinc",
              "LAC_lik_0_1_yr_exam5_wsp",
              "NOx_lik_0_1_yr_exam5_wght",
              "mds_0","mds_3plus")

list1 <- list(exposure)
mesa_results_list <- list(list1,list1,list1)
names(mesa_results_list) <- raceethnicity

for(i in raceethnicity){
  print("####################")
  print(i)
  for(j in exposure){
    print("####################")
    print(i)
    print(j)
    if(i=="whiteNH"&j=="mds_3plus"){
      print("skip mds3plus for white NH")
    } else {
      load(file=paste0(output.dir,"/",j,"_",i,"_mesa_ewas_svs.Robj"))
      a_priori_sites <- check_a_priori_ewas(ewas.mesa)
      mesa_results_list[[i]][[j]] <- a_priori_sites
      for(k in names(a_priori_sites)){
        print(k)
        print(dim(a_priori_sites[[k]]))
      }
    }
  }
}

save(mesa_results_list,file=paste0(output.dir,"/mesa_a_priori_sites.Rdata"))

##########

## MBMS

print("#######   MBMS  ##########")

raceethnicity <- c("whiteNH","blackNH")
exposure <- c("jim_crow_birth_state",
              "parent_educ_low","parent_educ_mid",
              "educ_low","educ_mid",
              "hh_income_pov_ratio",
              "ICEraceinc",
              "black_carbon_yearly_avg",
              "nitrous_oxides_pollution_proximity_index",
              "eod_0","eod_3plus")

list1 <- list(exposure)
mbms_results_list <- list(list1,list1)
names(mbms_results_list) <- raceethnicity

for(i in raceethnicity){
  print("####################")
  print(i)
  for(j in exposure){
    print("####################")
    print(j)
    load(file=paste0(output.dir,"/",j,"_",i,"_mbms_ewas.Robj"))
    a_priori_sites <- check_a_priori_ewas(ewas.mbms)
    mbms_results_list[[i]][[j]] <- a_priori_sites
    for(k in names(a_priori_sites)){
      print(k)
      print(dim(a_priori_sites[[k]]))
    }
  }
}
save(mbms_results_list,file=paste0(output.dir,"mbms_a_priori_sites.Rdata"))

mesa_whitenh <- mesa_results_list$whiteNH[2:9]
mesa_whitenh_apriori_results <- data.frame(matrix(nrow = 0,ncol=14))
for(i in names(mesa_whitenh)){
  temp <- mesa_whitenh[[i]]
  temp_df <- rbindlist(temp)
  temp_df$exposure <- i
  colnames(mesa_whitenh_apriori_results) <- colnames(temp_df)
  mesa_whitenh_apriori_results <- rbind(mesa_whitenh_apriori_results,temp_df)
}
write.csv(mesa_whitenh_apriori_results,file="mesa_whitenh_apriori_output.csv")

mesa_blacknh <- mesa_results_list$blackNH[2:10]
mesa_blacknh_apriori_results <- data.frame(matrix(nrow = 0,ncol=14))
for(i in names(mesa_blacknh)){
  temp <- mesa_blacknh[[i]]
  temp_df <- rbindlist(temp)
  temp_df$exposure <- i
  colnames(mesa_blacknh_apriori_results) <- colnames(temp_df)
  mesa_blacknh_apriori_results <- rbind(mesa_blacknh_apriori_results,temp_df)
}
write.csv(mesa_blacknh_apriori_results,file="mesa_blacknh_apriori_output.csv")


mbms_whitenh <- mbms_results_list$whiteNH[2:10]
mbms_whitenh_apriori_results <- data.frame(matrix(nrow = 0,ncol=14))
for(i in names(mbms_whitenh)){
  temp <- mbms_whitenh[[i]]
  temp_df <- rbindlist(temp)
  temp_df$exposure <- i
  colnames(mbms_whitenh_apriori_results) <- colnames(temp_df)
  mbms_whitenh_apriori_results <- rbind(mbms_whitenh_apriori_results,temp_df)
}
write.csv(mbms_whitenh_apriori_results,file="mbms_whitenh_apriori_output.csv")

mbms_blacknh <- mbms_results_list$blackNH[2:10]
mbms_blacknh_apriori_results <- data.frame(matrix(nrow = 0,ncol=14))
for(i in names(mbms_blacknh)){
  temp <- mbms_blacknh[[i]]
  temp_df <- rbindlist(temp)
  temp_df$exposure <- i
  colnames(mbms_blacknh_apriori_results) <- colnames(temp_df)
  mbms_blacknh_apriori_results <- rbind(mbms_blacknh_apriori_results,temp_df)
}
write.csv(mbms_blacknh_apriori_results,file="mbms_blacknh_apriori_output.csv")
