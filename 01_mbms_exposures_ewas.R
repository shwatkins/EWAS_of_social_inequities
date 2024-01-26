# EWAS of exposures
# MBMS
# Stratified by race/ethnicity
# started 27/01/22

# Set directories:

project.dir <- ""

data.dir <- file.path(project.dir, "")
output.dir <- file.path(project.dir, "")

output.dir

## http://github.com/perishky/meffil
library(meffil)
options(mc.cores=22)

library(readxl)

# this file contains info on sentrix ID plus study ID, recruitment site, amount of DNA, and a few other variables
pheno <- read.csv(file=paste0(data.dir,"/mbms_for_bristol_2021_08_20.csv"),header = T)
pheno <- pheno[!is.na(pheno$sentrix_id),]
rownames(pheno) <- pheno$sentrix_id

# load in meth so we can subset to participants who have dnam data
load(file=paste0(data.dir,"/norm.beta_mbms_slide_fixed.rda"))
meth <- ret
meth <- as.data.frame(meth)
meth.m <- as.matrix(meth)
meth_samples <- as.character(colnames(meth))
pheno <- pheno[meth_samples,]
pheno$ParticipantID <- as.character(pheno$ParticipantID)

# sort out variables
vars <- c("smoke_3cat","race","gender","eod","eds","hh_income_pov_ratio",
          "ICEraceinc","black_carbon_yearly_avg",
          "nitrous_oxides_pollution_proximity_index","educ","jim_crow_birth_state",
          "birth_state_policy_liberalism", "diastolic_bp",
          "systolic_bp","framingham_cvd_score","diabetic","metabolic_syndrome_score")

# this file contains the cleaned EWAS exposures and outcomes (not sentrix ID)
ewas_vars <- read.csv(file=paste0(data.dir,"/mbms_measures_for_ewas_jan22.csv"),header = T)
ewas_vars$id <- as.character(ewas_vars$id)

# merge those two pheno files together by MBMS study ID
pheno <- merge(pheno,ewas_vars,by.x="ParticipantID",by.y="id")
rownames(pheno) <- pheno$sentrix_id

# sort out the few variables that were duplicated between those two pheno files
pheno$race <- pheno$race.x
pheno$gender <- pheno$gender.x
drops <- c("race.x","gender.x","race.y","gender.y","educ_r")
pheno <- pheno[,!colnames(pheno) %in% drops]

# sort out smoking var to current/former/never
pheno$smoke_3cat <- pheno$Smoke_now_8hrs
class(pheno$smoke_3cat)
pheno$smoke_3cat <- as.character(pheno$smoke_3cat)
pheno$smoke_3cat[pheno$smoke_3cat == "Current, in last 8hrs"] <- "Current smoker"
pheno$smoke_3cat[pheno$smoke_3cat == "Current, not in 8hrs"] <- "Current smoker"
table(pheno$smoke_3cat)
pheno$smoke_3cat <- as.factor(pheno$smoke_3cat)
table(pheno$smoke_3cat)
pheno$smoke_3cat <- ordered(pheno$smoke_3cat,levels=c("Never smoker","Ex-smoker","Current smoker"))
table(pheno$smoke_3cat)

pheno$jim_crow_birth_state <- as.factor(as.character(pheno$jim_crow_birth_state))

# "4+ yrs college" is the reference level for educ
# creat two binary educ columns, each with a non-reference level set to NA
pheno$educ <- as.character(pheno$educ)
pheno$educ_low <- NA
pheno$educ_low[pheno$educ == "4+ yrs college"] <- 0
pheno$educ_low[pheno$educ == "< HS"] <- 1
print(table(pheno$educ_low))
pheno$educ_mid <- NA
pheno$educ_mid[pheno$educ == "4+ yrs college"] <- 0
pheno$educ_mid[pheno$educ == ">=HS and <4 yrs college"] <- 1
print(table(pheno$educ_mid))

# reduce phenotype df to participants with DNAm data
pheno <- pheno[meth_samples,]

# check participants are in the same order in oheno and meth dfs
identical(rownames(pheno),colnames(meth))

# add in cateorical EOD and parental education (a later dataset version)
eod_categorical <- read.csv(file=paste0(data.dir,"/mbms_measures_for_ewas_cat_EOD.csv"),header = T)
# reduce to ID, EOD and parental educ:
eod_categorical <- eod_categorical[,c(1,5,22)]
class(eod_categorical$eod)
colnames(eod_categorical) <- c("id","eod_cat","parent_educ")

pheno <- merge(pheno,eod_categorical,by.x="ParticipantID",by.y="id")
print(head(pheno))
rownames(pheno) <- pheno$sentrix_id

# "1-2" is the reference level for eod
# creat two binary eod columns, each with a non-reference level set to NA
pheno$eod_0 <- NA
pheno$eod_0[pheno$eod_cat == "1-2"] <- 0
pheno$eod_0[pheno$eod_cat == "0"] <- 1
pheno$eod_3plus <- NA
pheno$eod_3plus[pheno$eod_cat == "1-2"] <- 0
pheno$eod_3plus[pheno$eod_cat == "3+"] <- 1

# "4+ yrs college" is the reference level for parent educ
# creat two binary educ columns, each with a non-reference level set to NA
pheno$parent_educ <- as.character(pheno$parent_educ)
pheno$parent_educ_low <- NA
pheno$parent_educ_low[pheno$parent_educ == "4+ yrs college"] <- 0
pheno$parent_educ_low[pheno$parent_educ == "< HS"] <- 1
pheno$parent_educ_mid <- NA
pheno$parent_educ_mid[pheno$parent_educ == "4+ yrs college"] <- 0
pheno$parent_educ_mid[pheno$parent_educ == ">= HS and <4 yrs college"] <- 1

pheno <- pheno[meth_samples,]

# add in cell counts
cellcounts<-meffil.estimate.cell.counts.from.betas(meth.m,cell.type.reference="blood gse35069 complete",verbose=T)
cellcounts<-data.frame(IID=row.names(cellcounts),cellcounts)
cellcounts <- cellcounts[,-1]
# make sure ID are row names!
identical(rownames(pheno), rownames(cellcounts))
cellcount.names <- as.character(colnames(cellcounts))
pheno <- merge(pheno,cellcounts,by="row.names")
rownames(pheno) <- pheno$Row.names
pheno <- pheno[,-1]

# load in batch vars
load(file=paste0(data.dir,"/samplesheet_3.Rdata"))
samplesheet_3 <- samplesheet_3[,c("participant", "Slide", "plate","dna","good")]
samplesheet_3$Slide <- as.factor(as.character(samplesheet_3$Slide))
sum(is.na(samplesheet_3))
identical(rownames(samplesheet_3), colnames(meth)) 
identical(rownames(samplesheet_3), rownames(pheno))
pheno <- merge(pheno,samplesheet_3,by="row.names")
rownames(pheno) <- pheno$Row.names
pheno <- pheno[,-1]

ewas_covars <- c("AgeYRS","gender","Slide")
ewas_covars <- c(ewas_covars,cellcount.names)
ewas_covars_basic <- c("AgeYRS","Slide")
ewas_covars_cellcounts <- c(ewas_covars_basic,cellcount.names)

identical(rownames(pheno),colnames(meth))

# remove X and Y chromosomes
probeDetails <- meffil.featureset("epic")
rownames(probeDetails) <- probeDetails$name
probeDetails <- na.omit(probeDetails)
rownames(probeDetails) <- probeDetails$name
xchr <- probeDetails[probeDetails$chromosome=="chrX",]
ychr <- probeDetails[probeDetails$chromosome=="chrY",]
xchr <- as.character(xchr$name)
ychr <- as.character(ychr$name)
XY <- c(xchr, ychr)
dat_cpgs <- as.character(rownames(meth))
XY <- XY[XY %in% dat_cpgs]
XY <- na.omit(XY)
meth.all <- as.data.frame(meth)[! rownames(meth) %in% XY,]
meth.all <- as.matrix(meth.all)

##### prepare ewas #####

exposures <- c("jim_crow_birth_state",
               "parent_educ_low","parent_educ_mid",
               "educ_low","educ_mid",
               "hh_income_pov_ratio",
               "ICEraceinc",
               "black_carbon_yearly_avg",
               "nitrous_oxides_pollution_proximity_index",
               "eod_0","eod_3plus")

# set meffil EWAS parameters
ewas.parameters <- meffil.ewas.parameters(sig.threshold=9e-8,  ## EWAS p-value threshold
                                          max.plots=10, ## plot at most 100 CpG sites
                                          qq.inflation.method="median",  ## measure inflation using median
                                          model="all") ## select default EWAS model; 

#########

# stratify participants by racialized group
# this is because most of our exposures don't mean the same thing for the two social groups

pheno.b <- pheno[pheno$race == "Black N.H.",]
blackNH_participants <- as.character(rownames(pheno.b))
meth.b <- meth.all[,blackNH_participants]
identical(rownames(pheno.b),colnames(meth.b))
meth.b <- as.matrix(meth.b)

pheno.w <- pheno[pheno$race == "White N.H.",]
whiteNH_participants <- as.character(rownames(pheno.w))
meth.w <- meth.all[,whiteNH_participants]
identical(rownames(pheno.w),colnames(meth.w))
meth.w <- as.matrix(meth.w)

#######

# Run EWAS #

# adjusted for age, smoking, and cell counts; as well as SVs (for batch)

set.seed(23)

# black n.h. participants

for(i in exposures){
  ewas_covars_basic <- c("AgeYRS","smoke_3cat","gender")
  ewas_covars_cellcounts <- c(ewas_covars_basic,cellcount.names)
  print(i)
  if(is.numeric(pheno.b[,i])){
    print(summary(pheno.b[,i]))
  } else {
    print(table(pheno.b[,i]))
  }
  print("n of NAs:")
  print(sum(is.na(pheno.b[,i])))
  # get svs
  print("remove missing cases:")
  pheno.temp <- pheno.b[!is.na(pheno.b[,i]),]
  participants.temp <- as.character(rownames(pheno.temp))
  meth.temp <- meth.b[,participants.temp]
  meth.temp <- meffil.handle.outliers(meth.temp, winsorize.pct=0.05, outlier.iqr.factor=NA)
  meth.temp <- meffil:::impute.matrix(meth.temp, margin=1)    
  var.sites <- meffil.most.variable.cpgs(meth.temp, n=50000, sites=NULL, samples=NULL, winsorize.pct=0.05, outlier.iqr.factor=NA)
  beta <- meth.temp[var.sites,,drop=F]
  cov.frame <- model.frame(~., data.frame(pheno.temp[,ewas_covars_cellcounts], stringsAsFactors=F), na.action=na.pass)
  mod0 <- model.matrix(~., cov.frame)
  mod <- cbind(mod0, pheno.temp[,i])
  colnames(mod)[13] <- "variable"
  sva.blacknh <- sva(beta, mod=mod, mod0=mod0)#, n.sv=NULL)
  dat <- as.data.frame(sva.blacknh$sv)
  names(dat) <- paste0("sv",seq(1:ncol(dat)))
  rownames(dat) <- participants.temp
  save(dat,file=paste0(output.dir,i,"_blackNH_mbms_svs.Rdata"))
  write.csv(dat,file=paste0(output.dir,i,"_blackNH_mbms_svs.csv"),quote = F)
  # now add the first 5 SVs to the EWAS
  dat <- dat[,1:5]
  pheno.temp <- cbind(pheno.temp,dat)
  print(dim(pheno.temp))
  print(head(pheno.temp))
  ewas_covars_basic <- c("AgeYRS","smoke_3cat","gender","sv1","sv2","sv3","sv4","sv5")
  ewas_covars_cellcounts <- c(ewas_covars_basic,cellcount.names)

  ewas.mbms <- meffil.ewas(meth.temp, variable=pheno.temp[,i], covariates=pheno.temp[,ewas_covars_cellcounts], isva=F, random.seed=23) 
  save(ewas.mbms, file=paste0(output.dir,i,"_blackNH_mbms_ewas.Robj"))
  ewas.summary<-meffil.ewas.summary(ewas.mbms,meth.temp,parameters=ewas.parameters)                              
  meffil.ewas.report(ewas.summary, output.file=paste0(output.dir,i,"_blackNH_mbms_summary.html"))
  
}

# white n.h. participants #

for(i in exposures){
  ewas_covars_basic <- c("AgeYRS","smoke_3cat","gender")
  ewas_covars_cellcounts <- c(ewas_covars_basic,cellcount.names)
  print(i)
  if(is.numeric(pheno.w[,i])){
    print(summary(pheno.w[,i]))
  } else {
    print(table(pheno.w[,i]))
  }
  print("n of NAs:")
  print(sum(is.na(pheno.w[,i])))
  # get svs
  print("remove missing cases:")
  pheno.temp <- pheno.w[!is.na(pheno.w[,i]),]
  participants.temp <- as.character(rownames(pheno.temp))
  meth.temp <- meth.w[,participants.temp]
  meth.temp <- meffil.handle.outliers(meth.temp, winsorize.pct=0.05, outlier.iqr.factor=NA)
  meth.temp <- meffil:::impute.matrix(meth.temp, margin=1)    
  var.sites <- meffil.most.variable.cpgs(meth.temp, n=50000, sites=NULL, samples=NULL, winsorize.pct=0.05, outlier.iqr.factor=NA)
  beta <- meth.temp[var.sites,,drop=F]
  cov.frame <- model.frame(~., data.frame(pheno.temp[,ewas_covars_cellcounts], stringsAsFactors=F), na.action=na.pass)
  mod0 <- model.matrix(~., cov.frame)
  mod <- cbind(mod0, pheno.temp[,i])
  sva.whitenh <- sva(beta, mod=mod, mod0=mod0, n.sv=NULL)
  dat <- as.data.frame(sva.whitenh$sv)
  names(dat) <- paste0("sv",seq(1:ncol(dat)))
  rownames(dat) <- participants.temp
  save(dat,file=paste0(output.dir,i,"_whiteNH_mbms_svs.Rdata"))
  write.csv(dat,file=paste0(output.dir,i,"_whiteNH_mbms_svs.csv"),quote = F)
  # now add the first 5 SVs to the EWAS
  # but for some of the white NH EWAS there are less than 5 SVs
  if(ncol(dat)>4){
    dat <- dat[,1:5]
    sv_names <- as.character(colnames(dat))
    ewas_covars_basic <- c("AgeYRS","smoke_3cat","gender")
    ewas_covars_basic <- c(ewas_covars_basic,sv_names)
  } else {
    dat <- dat[,1:ncol(dat)]
    sv_names <- as.character(colnames(dat))
    ewas_covars_basic <- c("AgeYRS","smoke_3cat","gender")
    ewas_covars_basic <- c(ewas_covars_basic,sv_names)
  }
  pheno.temp <- cbind(pheno.temp,dat)
  ewas_covars_cellcounts <- c(ewas_covars_basic,cellcount.names)

  ewas.mbms <- meffil.ewas(meth.temp, variable=pheno.temp[,i], covariates=pheno.temp[,ewas_covars_cellcounts], isva=F, random.seed=23) 
  save(ewas.mbms, file=paste0(output.dir,i,"_whiteNH_mbms_ewas.Robj"))
  ewas.summary<-meffil.ewas.summary(ewas.mbms,meth.temp,parameters=ewas.parameters)                              
  meffil.ewas.report(ewas.summary, output.file=paste0(output.dir,i,"_whiteNH_mbms_summary.html"))
}



