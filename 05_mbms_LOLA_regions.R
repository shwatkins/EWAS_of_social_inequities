# Run the LOLA region set enrichment test for MBMS
# based on a script provided by Josine Min for my PhD
# https://github.com/shwatkins/PhD/blob/master/450k_correlation_analysis/11_F7_LOLA_tfbs.R
# uses region sets created by the LOLA team, available through http://lolaweb.databio.org.

library(data.table)
library(LOLA)
library(simpleCache)
library(reshape2)
library(GenomicAlignments)
library(Rsamtools)
library(biovizBase)
library(meffil)
library(dplyr)
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")
theme_set(theme_bw())

# for plots
library(ggplot2)
library(viridis)
library(scales)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(11,"RdBu"))
colourCount = 125

project.dir <- ""
data.dir <- file.path(project.dir, "","")
output.dir <- file.path(project.dir, "")

##load cell type conversion and colors
cellType_conversions=fread("CellTypes.tsv",drop="collection")
colors=fread("color.tsv")

##load regiondb
regionDB <- loadRegionDB("LOLA/gene_annotation", collection="7regions")
head(regionDB)

# get list of all possible cpgs
load(paste0(output.dir,"hh_income_pov_ratio_blackNH_mbms_ewas.Robj"))
temp <- ewas.mbms$p.value
allcgs <- as.character(rownames(temp))
retaincpg<-allcgs

# prep illumina locations
print("prep illumina locations")
Illumina850 <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations
Illumina850_dt=as.data.table(Illumina850)
Illumina850_dt[,cpgID:=row.names(Illumina850),]
Illumina850_dt <- Illumina850_dt[Illumina850_dt$cpgID%in%retaincpg,]
Illumina850_dt[,cpgstart_pre:=ifelse(strand=="-",pos-100,pos-99),]
Illumina850_dt[,cpgend_pre:=ifelse(strand=="-",pos+100,pos+101),]

#collapse overlaps
print("collapse overlaps")
gr_range = with(Illumina850_dt,GRanges(seqnames=chr,ranges=IRanges(cpgstart_pre,cpgend_pre)))
gr_cpg = with(Illumina850_dt,GRanges(seqnames=chr,ranges=IRanges(pos,pos)))
overlap=as.data.table(findOverlaps(gr_cpg, gr_range))
overlap_red=overlap[,list(subjectHit=min(subjectHits),NsubjectHits=.N),by=queryHits]
Illumina850_dt[,cpgstart:=start(gr_range[overlap_red$subjectHit])]
Illumina850_dt[,cpgend:=end(gr_range[overlap_red$subjectHit])]
Illumina850_dt[,NsubjectHits:=overlap_red$NsubjectHits]
Illumina850_sub=Illumina850_dt[,c("cpgID","cpgstart","cpgend","pos"),with=FALSE]
setnames(Illumina850_sub,c("cpgID","pos"),c("cpg","ill_pos"))

process_LOLA = function (LOLA_res, collections=c("7regions"),cellType_conversions){
  
  LOLA_res=LOLA_res[!is.na(userSet)]
  
  LOLA_res=LOLA_res[collection%in%collections]
  #changed form exp() to 10^ to accomodate change in LOLA
  
  LOLA_res[,p.adjust:=p.adjust(10^(-pValueLog),method="BY"),by=userSet]
  LOLA_res[,mlog10p.adjust:=-log10(p.adjust),]
  
  ##standardize cellTypes
  LOLA_res[cellType==""|is.na(cellType),cellType:="Not defined",]
  LOLA_res=merge(LOLA_res,cellType_conversions,by="cellType",all.x=TRUE)
  ##correct wrong annotation
  LOLA_res[description=="T-cell acute lymphoblastic leukaemia (T-ALL) cell line.",c("Lineage1","Lineage","cellType_corr"):=list(Lineage1="Lymphoid",Lineage="Lymphoid",cellType_corr="T lymphocyte"),]
  LOLA_res[,lineage_count_allstate:=length(unique(filename[!is.na(filename)])),by=c("Lineage","cellState")]
  LOLA_res[,lineage_count_all:=length(unique(filename[!is.na(filename)])),by=c("Lineage")]
  print(head(LOLA_res))

  return(LOLA_res)  
  
}

raceethnicity <- c("whiteNH","blackNH")
exposure <- c("jim_crow_birth_state",
              "parent_educ_low","parent_educ_mid",
              "educ_low","educ_mid",
              "hh_income_pov_ratio",
              "ICEraceinc",
              "black_carbon_yearly_avg",
              "nitrous_oxides_pollution_proximity_index",
              "eod_0","eod_3plus")

#### 

print("loop start")
for(i in raceethnicity){
  print(i)
  regions_list <- list(exposure)
  for(j in exposure){
    print(i)
    print(j)
    # this file is generated in the 'mbms_ewas_catalog_enrichment' script
    load(paste0(output.dir,"/biological_investigation/",j,"_mbms_",i,"_topsites_details.Rdata"))
    top100 <- as.character(results_info$Row.names)
    data <- Illumina850_dt[Illumina850_dt$cpgID%in%top100,]
    ##prepare GoDMC data: merge CpGs and SNPs that are in proximity to eachother to avoid infalting the results, 1kb around cp
    setnames(data,c("cpgID"),c("cpg"))
    data=merge(data,Illumina850_sub,by="cpg",all.x=TRUE)
    #check pos (should be all true)
    print(table(data[,ill_pos==pos,]))
    ####run with external background for CpGs
    hg19_Illumina850_gr=with(Illumina850_dt, GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle(strand),ID=cpgID))
    print("done hg19_Illumina850_gr")
    seq_Illumina850=getSeq(BSgenome.Hsapiens.UCSC.hg19,hg19_Illumina850_gr)
    # add GC, CpG frequency of probes to the data table. Also add wheher or not the
    # probe is in the list of top correlated probes.
    Illumina850_dt[,GC_freq:=letterFrequency(seq_Illumina850, "CG", as.prob=T),]
    Illumina850_dt[,CpG_freq:=dinucleotideFrequency(seq_Illumina850, step=2, as.prob=T)[,"CG"],]
    Illumina850_dt[,isTop100:=ifelse(cpgID%in%rownames(top100),TRUE,FALSE),]
    Illumina850_dt<-Illumina850_dt[Illumina850_dt$cpgID%in%retaincpg,]
    
    F7_cpg_gr=unique(with(Illumina850_dt,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart_pre, end=cpgend_pre),strand=Rle("*"))))
    print("done F7_cpg_gr")
    
    background <- Illumina850_dt
    background <- unique(with(Illumina850_dt,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"))))
    print("done background")
    print(all.equal(data$cpgstart.x,data$cpgstart.y))
    print(all.equal(data$cpgend.x,data$cpgend.y))
    top100hits <- unique(with(data,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart.x, end=cpgend.x),strand=Rle("*"))))
    print("done top100hits")
    top100hits_bed <- as.data.frame(top100hits)
    top100hits_bed <- top100hits_bed[,c(1:3)]

    background_bed <- as.data.frame(background)
    background_bed <- background_bed[,c(1:3)]
    print("done background bed")
    
    dim(Illumina850_dt)
    dim(background)
    lola_res0_matched=runLOLA(top100hits, background, regionDB, cores=1)
    print("done runLOLA")
    lola_res0_matched$logOddsRatio<-log(lola_res0_matched$oddsRatio)
    dat_sort <- lola_res0_matched[order(lola_res0_matched$oddsRatio, decreasing = T),]
    dat_sort$pvalue <- 10^-dat_sort$pValueLog
    dat_sort <- dat_sort[order(dat_sort$pvalue),]
    dat_sort_pval <- dat_sort[dat_sort$pvalue < 0.05,]
    print("done dat_sort")
    
    lola_res0_matched$logOddsRatio
    lola_res0_matched$epigenome = as.character(lapply(strsplit(as.character(lola_res0_matched$description), split=" "),
                                                      tail, n=1))
    
    # change to cell types
    locResults=process_LOLA(LOLA_res=lola_res0_matched,cellType_conversions=cellType_conversions)
    print("done cell type conversion")
    # this is the file we need for further collapsed plots
    locResults$epigenome_brief <- locResults$epigenome
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_cpg_shelves"] <- "CpG island shelves"
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_genes_3UTRs"] <- "3'UTR"
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_cpg_inter"] <- "CpG inter"
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_genes_intergenic"] <- "Intergenic"
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_cpg_shores"] <- "CpG island shores"
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_genes_introns"] <- "Introns"
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_genes_promoters"] <- "Promoters"
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_genes_1to5kb"] <- "1 to 5kb from genes"
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_genes_intronexonboundaries"] <- "Intron-exon boundaries"
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_genes_5UTRs"] <- "5'UTR"
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_cpg_islands"] <- "CpG island"
    locResults$epigenome_brief[locResults$epigenome_brief == "hg19_genes_exons"] <- "Exons"
    
    save(locResults, file = paste0(output.dir,"/",i,"_",j,"_locresults_regions.Rdata"))
    regions_list[[j]] <- locResults
    #}
  }
  save(regions_list,file=paste0(output.dir,"/",i,"_regions_list.Rdata"))
  print(paste(i,"region list saved"))
  regions_list <- regions_list[2:12]
  
  get_blood_ors <- function(x){
    dat <- as.data.frame(x)
    dat_sort <- dat[order(dat$oddsRatio, decreasing = T),]
    dat_sort$pvalue <- 10^-dat_sort$pValueLog
    dat_sort <- dat_sort[order(dat_sort$pvalue),]
    dat_sort_brief <- dat_sort[,c("oddsRatio", "logOddsRatio", "pvalue", "epigenome", "Tissue")]
    return(dat_sort_brief)
  }
  
  blood_region_stuff <- lapply(regions_list, FUN=get_blood_ors)
  # list of length exposures with all region results for each exposure
  names(blood_region_stuff) <- exposure
  print("blood region stuff done")
  
  # get mean ORs
  # y = exposures
  df3_per_exposure <- list(exposure)
  for(y in names(blood_region_stuff)){
    print(y)
    df <- as.data.frame(blood_region_stuff[[y]])
    region_names <- as.character(df$epigenome)
    df$epigenome <- as.character(df$epigenome)
    region_names <- unique(region_names)
    df3 <- list()
    for (z in 1:length(region_names)){
      h <- df[df$epigenome == region_names[z],]
      mean_or <- mean(h$oddsRatio)
      mean_logor <- log(mean_or)
      mean_pval <- mean(h$pvalue)
      epigenome <- region_names[z]
      n_blood_regions <- nrow(h)
      s <- cbind(mean_or, mean_logor, mean_pval, epigenome, n_blood_regions)
      df3[[z]] <- as.data.frame(s)
    }
    df3_per_exposure[[y]] <- df3
    # df3 is a list of the mean OR for each region (so there aren't multiple regions)
  }
  print("df3 done")
  df3_per_exposure <- df3_per_exposure[2:12]
  # mix them together 
  region_df_list <- list()#(exposure)
  for (r in 1:length(df3_per_exposure)){
    region_df <- data.frame(matrix(unlist(df3_per_exposure[[r]]), nrow=length(df3_per_exposure[[r]]), byrow=T), stringsAsFactors = F)
    region_df$group <- names(df3_per_exposure)[[r]]
    colnames(region_df) <- c("oddsRatio", "logOddsRatio", "pvalue", "epigenome", "n_blood_regions", "group")
    region_df_list[[r]] <- region_df
  }
  print("region_df_list done")
  all_groups <- rbindlist(region_df_list)

  all_groups$logOddsRatio <- as.numeric(all_groups$logOddsRatio)
  all_groups$logOddsRatio[all_groups$oddsRatio == 0] <- 0
  all_groups$pvalue <- as.numeric(all_groups$pvalue)
  all_groups$significance <- ""
  all_groups$significance[all_groups$pvalue < 0.05] <- "*"
  all_groups$significance[all_groups$pvalue < 1e-5] <- "**"
  all_groups$significance[all_groups$pvalue < 1e-10] <- "***"
  print("results df tidied")
  ## reduce to only regions with significance
  all_groups_sig <- all_groups[!all_groups$significance == "",]
  class(all_groups_sig$epigenome)
  sig_regions <- all_groups_sig$epigenome
  sig_regions <- unique(sig_regions)
  all_groups_reduced <- all_groups[all_groups$epigenome %in% sig_regions,]
  save(all_groups,all_groups_reduced,file=paste0(output.dir,"/",i,"_regions_ready_to_plot.Rdata"))
  
}

################

# run plots

load(paste0(output.dir,"/whiteNH_regions_ready_to_plot.Rdata"))
all_groups_reduced$group <- factor(all_groups_reduced$group,levels=c("jim_crow_birth_state",
                                                                     "parent_educ_low","parent_educ_mid",
                                                                     "educ_low","educ_mid",
                                                                     "hh_income_pov_ratio",
                                                                     "ICEraceinc",
                                                                     "black_carbon_yearly_avg",
                                                                     "nitrous_oxides_pollution_proximity_index",
                                                                     "eod_0","eod_3plus"),
                                   labels = c("Jim_crow_birth","parent_educ_low","parent_educ_mid","educ_low","educ_mid","hh_inc_pov","ICEraceinc","black_carbon","NOx","eod_0","eod_3plus"))

# use max and min to set plot limits so that all the plots are on the same scale
max(all_groups_reduced$logOddsRatio)
min(all_groups_reduced$logOddsRatio)

jpeg(filename=paste0(output.dir,"/biological_investigation/LOLA/whiteNH_regions_heatmap.jpg"), width=5, height=3, units = "in", res = 600)
ggplot(data = all_groups_reduced, aes(group, epigenome))+
  geom_tile(color = "white", aes(fill=logOddsRatio))+ 
  scale_fill_gradientn(colors = getPalette(colourCount), 
                       values = rescale(c(-2.8,0,2.8)),
                       guide = "colorbar", limits=c(-2.8,2.8),
                       space = "Lab", 
                       name="Log Odds\nRatio") +
  geom_text(aes(label = significance, color = logOddsRatio > 1.5)) +
  labs(title="Genomic region enrichment\nMBMS White NH",
       x ="EWAS analysis", y = "Genomic region")+
  scale_color_manual(guide = FALSE, values = c("black", "white"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))

dev.off()

load(paste0(output.dir,"/blackNH_regions_ready_to_plot.Rdata"))
all_groups_reduced$group <- factor(all_groups_reduced$group,levels=c("jim_crow_birth_state",
                                                                     "parent_educ_low","parent_educ_mid",
                                                                     "educ_low","educ_mid",
                                                                     "hh_income_pov_ratio",
                                                                     "ICEraceinc",
                                                                     "black_carbon_yearly_avg",
                                                                     "nitrous_oxides_pollution_proximity_index",
                                                                     "eod_0","eod_3plus"),
                                   labels = c("Jim_crow_birth","parent_educ_low","parent_educ_mid","educ_low","educ_mid","hh_inc_pov","ICEraceinc","black_carbon","NOx","eod_0","eod_3plus"))

# use max and min to set plot limits so that all the plots are on the same scale
max(all_groups_reduced$logOddsRatio)
min(all_groups_reduced$logOddsRatio)

jpeg(filename=paste0(output.dir,"/blackNH_regions_heatmap.jpg"), width=5, height=3, units = "in", res = 600)
ggplot(data = all_groups_reduced, aes(group, epigenome))+
  geom_tile(color = "white", aes(fill=logOddsRatio))+ 
  scale_fill_gradientn(colors = getPalette(colourCount),  
                       values = rescale(c(-2.8,0,2.8)),
                       guide = "colorbar", limits=c(-2.8,2.8),
                       space = "Lab", 
                       name="Log Odds\nRatio") +
  geom_text(aes(label = significance, color = logOddsRatio > 1.5)) +
  labs(title="Genomic region enrichment\nMBMS Black NH",
       x ="EWAS analysis", y = "Genomic region")+
  scale_color_manual(guide = FALSE, values = c("black", "white"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))

dev.off()


