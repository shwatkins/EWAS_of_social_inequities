# mbms - test biological enrichments using the EWAS catalog

## grab ewas catalog df 
library(data.table)
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(meffil)
library(epitools)
library(viridis)
library(ggpubr)
library(ggplot2)

project.dir <- ""
data.dir <- file.path(project.dir, "","")
output.dir <- file.path(project.dir, "")

probedetails_epic <- meffil.featureset("epic")
rownames(probedetails_epic) <- probedetails_epic$name
probedetails_epic <- na.omit(probedetails_epic)

probedetails_450 <- meffil.featureset("450k")
rownames(probedetails_450) <- probedetails_450$name
probedetails_450 <- na.omit(probedetails_450)

# this file is created by the 'ewas_catalog_category_coding' script
load("/ewas_catalog_for_enrichment.Rdata")
# contains objects: EWAS_catalog,EWAS_catalog_collapsed_freq_table

## now look at top hits stratified by racialized group

## 1. mbms

raceethnicity <- c("whiteNH","blackNH")
exposure <- c("jim_crow_birth_state",
              "parent_educ_low","parent_educ_mid",
              "educ_low","educ_mid",
              "hh_income_pov_ratio",
              "ICEraceinc",
              "black_carbon_yearly_avg",
              "nitrous_oxides_pollution_proximity_index",
              "eod_0","eod_3plus")
list2 <- list(exposure)
plot_list <- list(list2,list2)
names(plot_list) <- raceethnicity
or_output_list <- list(list2,list2)
names(or_output_list) <- raceethnicity

GO_list <- list(list2,list2,list2)
names(GO_list) <- raceethnicity

for(i in raceethnicity){
  print(i)
  for(j in exposure){
    print(j)
    load(file=paste0(output.dir,"/",j,"_",i,"_mbms_ewas.Robj"))
    results <- ewas.mbms$analyses$all$table
    results <- results[order(results$p.value),]
    cpgs <- as.character(rownames(results))
    top_results <- results[1:100,]
    top_results <- top_results[order(top_results$p.value),]
    test.cpgs <- as.character(rownames(top_results))
    save(top_results,file = paste0(output.dir,"/",j,"_mbms_",i,"_topsites.Rdata"))
    ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
    gst <- gometh(sig.cpg = test.cpgs, all.cpg = cpgs, array.type = "EPIC", collection = "GO", plot.bias = F, prior.prob = TRUE)
    print(table(gst$FDR<0.05))
    x <- table(gst$FDR<0.05)
    GO_list[[i]][[j]] <- x
    table_go <- gst[order(gst$FDR),]    
    print(head(table_go))
    # you can also look at KEGG pathway enrichment
    # but currently the link the function uses is broken so it doesn't work
    ####    kegg <- gometh(sig.cpg = test.cpgs, all.cpg = cpgs, collection = "KEGG", plot.bias = F, prior.prob=TRUE)
    ####    # Table of top KEGG results
    ####    print(table(kegg$FDR<0.05))
    ####    table_kegg <- kegg[order(kegg$FDR),]    
    ####    print(head(kegg))
    
    # add in extra info about cpg sites from meffil and save out - this will be used in further scripts
    test.cpgs <- as.character(rownames(top_results))
    results_info <- merge(top_results,probedetails_epic,by.x="row.names",by.y = "row.names")
    results_info <- results_info[order(results_info$p.value),]
    save(results_info,file = paste0(output.dir,"/",j,"_mbms_",i,"_topsites_details.Rdata"))
    write.csv(results_info,file = paste0(output.dir,"/",j,"_mbms_",i,"_topsites_details.csv"),quote = F)
    
    results_info$gene <- sapply(strsplit(results_info$gene.symbol,';'), "[", 1)
    genes <- results_info$gene
    genes.tab <- data.frame(table(genes))
    genes.tab <- genes.tab[order(genes.tab$Freq,decreasing = T),]
    write.csv(genes.tab,file=paste0(output.dir,"/genes_mbms_",i,"_",j,".csv"),quote = F, row.names = F)
    
    # test EWAS catalog enrichment
    # get top 100 cpgs
    top100 <- rownames(results[1:100,])
    top_100_df <- results[1:100,]
    test_ewas_cat <- EWAS_catalog[EWAS_catalog$CpG %in% top100,]
    phenotypes <- unique(test_ewas_cat$phenotype)
    ewas_cat_hits <- list()
    for(k in 1:length(phenotypes)){
      x <- phenotypes[k]
      print(x)
      temp<-test_ewas_cat[test_ewas_cat$phenotype==x,]
      test <- unique(temp$CpG)
      temp <- temp[duplicated(temp$CpG)==F,]
      tab.x <- nrow(temp)
      ewas_cat_hits[[k]] <- tab.x
    }
    print(phenotypes)
    
    names(ewas_cat_hits) <- phenotypes
    
    print(ewas_cat_hits)
    
    hits_for_enrichment <- do.call("rbind", ewas_cat_hits)
    hits_for_enrichment <- as.data.frame(hits_for_enrichment)
    hits_for_enrichment$names <- rownames(hits_for_enrichment)
    names(hits_for_enrichment) <- c("freq","var")
    
    ewas_cat_condensed <- as.data.frame(EWAS_catalog_collapsed_freq_table)
    names(ewas_cat_condensed) <- c("var","freq")
    ewas_cat_condensed <- ewas_cat_condensed[ewas_cat_condensed$var %in% hits_for_enrichment$var,]
    or_list <- list()
    for(s in 1:length(phenotypes)){
      x <- phenotypes[s]
      print(x)
      test_table <- data.frame(matrix(0, nrow = 2, ncol = 2))
      test_table[1,1] <- hits_for_enrichment$freq[hits_for_enrichment$var == x]
      test_table[2,1] <- nrow(top_100_df)-hits_for_enrichment$freq[hits_for_enrichment$var == x]
      test_table[1,2] <- ewas_cat_condensed$freq[ewas_cat_condensed$var == x]
      test_table[2,2] <- 484781-ewas_cat_condensed$freq[ewas_cat_condensed$var == x]
      colnames(test_table) <- c("ewas_hits","full_cat")
      rownames(test_table) <- c("in_cat","not_in_cat")
      test_table <- as.matrix(test_table)
      or_list[[s]] <- oddsratio.fisher(test_table,rev="both")
      print(or_list[[s]])
      
    }
    
    # put all the or outputs in a single df
    or_output <- data.frame(matrix(ncol=6,nrow=0))
    names(or_output) <- c("estimate","lower","upper","midp.exact","fisher.exact","chi.square")
    for(f in 1:length(or_list)){
      or_odds <- as.data.frame(t(or_list[[f]]$measure[2,]))
      or_p <- as.data.frame(t(or_list[[f]]$p.value[2,]))
      or_together <- cbind(or_odds,or_p)
      or_output <- rbind(or_output,or_together)
    }
    rownames(or_output) <- phenotypes
    or_output$phenotype <- as.character(rownames(or_output))
    or_output$phenotype <- as.factor(or_output$phenotype)
    # we can choose remove exposures with silly confidence intervals witht he following line:
    #or_output <- or_output[or_output$upper < 50,]
    or_output_list[[i]][[j]] <- or_output
    print(or_output)
    # make the plot
    adjusted_p <- 0.05/nrow(or_output)
    fp <- ggplot()+
      geom_pointrange(data=or_output, aes(x=estimate, y=phenotype, xmin=lower, xmax=upper)) + 
      geom_pointrange(data=or_output[or_output$fisher.exact<0.05,], aes(x=estimate, y=phenotype, xmin=lower, xmax=upper,colour="#1F968BFF")) + 
      geom_vline(xintercept=1, linetype=2) +  # add a dotted line at x=1 after flip
      labs(title=paste0("mbms, ",i,", ",j), y = "EWAS catalog category", x = "Odds of enrichment (95% CI)")+
      theme_bw()  # use a white background
    
    plot_list[[i]][[j]] <- fp
    
  }
  
}

save(plot_list, file=paste0(output.dir,"/ewascat_enrichment_plot_list.Rdata"))
save(or_output_list, file=paste0(output.dir,"/ewascat_enrichment_or_output_list.Rdata"))
save(GO_list, file=paste0(output.dir,"/mbms_GO_list.Rdata"))


###############


### make composite plot

exposure <- c("jim_crow_birth_state",
              "parent_educ_low","parent_educ_mid",
              "educ_low","educ_mid",
              "hh_income_pov_ratio",
              "ICEraceinc",
              "black_carbon_yearly_avg",
              "nitrous_oxides_pollution_proximity_index",
              "eod_0","eod_3plus")

plot_list <- list()
for(i in exposure){
  print(i)
  or_output <- or_output_list[[1]][[i]]
  print(head(or_output))
  fp <- ggplot(data=or_output, aes(x=estimate, y=phenotype, xmin=lower, xmax=upper)) +
    geom_pointrange()+ 
    geom_pointrange(data=or_output[or_output$fisher.exact<0.05,], aes(x=estimate, y=phenotype, xmin=lower, xmax=upper),colour="#1F968BFF") + 
    geom_vline(xintercept=1, linetype=2) +  # add a dotted line at x=1 after flip
    labs(title=paste0(i), y = "EWAS catalog category", x = "Enrichment odds (95% CI)")+
    theme_bw()  # use a white background
  plot_list[[i]] <- fp
}
plot.out <- ggarrange(plotlist=plot_list,
                      labels = c("A","B","C","D","E","F","G","H","I","J","K"),
                      ncol = 4, nrow = 3)
annotate_figure(plot.out, top = text_grob("MBMS white NH",
                                          color = "black", face = "bold", size = 14))

jpeg(filename = paste0(output.dir,"/ewascat_enrichment_plot_composite_whiteNH.jpg"),width = 13, height = 9, units = "in", res = 600)
print(plot.out)
dev.off()

plot_list <- list()
for(i in exposure){
  print(i)
  or_output <- or_output_list[[2]][[i]]
  print(head(or_output))
  fp <- ggplot(data=or_output, aes(x=estimate, y=phenotype, xmin=lower, xmax=upper)) +
    geom_pointrange()+ 
    geom_pointrange(data=or_output[or_output$fisher.exact<0.05,], aes(x=estimate, y=phenotype, xmin=lower, xmax=upper),colour="#1F968BFF") + 
    geom_vline(xintercept=1, linetype=2) +  # add a dotted line at x=1 after flip
    labs(title=paste0(i), y = "EWAS catalog category", x = "Enrichment odds (95% CI)")+
    #coord_flip() + # flip coordinates (puts labels on y axis)
    theme_bw()  # use a white background
  plot_list[[i]] <- fp
}
plot.out <- ggarrange(plotlist=plot_list,
                      labels = c("A","B","C","D","E","F","G","H","I","J","K"),
                      ncol = 4, nrow = 3)
jpeg(filename = paste0(output.dir,"/ewascat_enrichment_plot_composite_blackNH.jpg"),width = 13, height = 9, units = "in", res = 600)
print(plot.out)
dev.off()
