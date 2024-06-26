
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("R.oo", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggtranscript", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggpubr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)


Assign_TSS_and_linked_peaks = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform out ----
  
  tracking_genes = unlist(strsplit(opt$tracking_genes, split=","))
  
  cat("tracking_genes_\n")
  cat(sprintf(as.character(tracking_genes)))
  cat("\n")
  
  #### Read the K562_Regulatory_Build ----
  
  K562_Regulatory_Build<-readGFF(file=opt$K562_Regulatory_Build)
 

  cat("K562_Regulatory_Build_0\n")
  cat(str(K562_Regulatory_Build))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_Regulatory_Build$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_Regulatory_Build$type))))
  cat("\n")
 
  # gr_K562_Regulatory_Build <- GRanges(
  #   seqnames = as.character(K562_Regulatory_Build$seqid),
  #   name2=as.character(K562_Regulatory_Build$feature_type),
  #   name3=as.character(K562_Regulatory_Build$activity),
  #   ranges=IRanges(
  #     start=as.numeric(K562_Regulatory_Build$start),
  #     end=as.numeric(K562_Regulatory_Build$end),
  #     name=K562_Regulatory_Build$regulatory_feature_stable_id))
  
  #### Read the K562_Promoters_in_TSS ----
  
  K562_Promoters_in_TSS<-readRDS(file=opt$K562_Promoters_in_TSS)
  
  
  cat("K562_Promoters_in_TSS_0\n")
  cat(str(K562_Promoters_in_TSS))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_Promoters_in_TSS$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_Promoters_in_TSS$type))))
  cat("\n")
  
  
  gr_K562_Promoters_in_TSS <- GRanges(
    seqnames = as.character(K562_Promoters_in_TSS$chr),
    name2=as.character(K562_Promoters_in_TSS$Symbol),
    name3=as.character(K562_Promoters_in_TSS$ensembl_gene_id),
    name4=as.character(K562_Promoters_in_TSS$strand),
    name5=as.character(K562_Promoters_in_TSS$feature),
    name6=as.character(K562_Promoters_in_TSS$activity),
    name7=as.character(K562_Promoters_in_TSS$CLASS),
    name8=as.character(K562_Promoters_in_TSS$gene_biotype),
    ranges=IRanges(
      start=as.numeric(K562_Promoters_in_TSS$start),
      end=as.numeric(K562_Promoters_in_TSS$end),
      name=K562_Promoters_in_TSS$value_string))
  
  
  cat("gr_K562_Promoters_in_TSS_0\n")
  cat(str(gr_K562_Promoters_in_TSS))
  cat("\n")
  
  
  #### Read the Peak to genes results file ----
  
  TSS_PoI<-readRDS(file=opt$TSS_PoI)
  

  cat("TSS_PoI_0\n")
  cat(str(TSS_PoI))
  cat("\n")
  cat(sprintf(as.character(names(summary(TSS_PoI$strand)))))
  cat("\n")
  cat(sprintf(as.character(summary(TSS_PoI$strand))))
  cat("\n")
  
  gr_TSS_PoI <- GRanges(
    seqnames = as.character(TSS_PoI$chr),
    name2=as.character(TSS_PoI$TSS_Symbol),
    name3=as.character(TSS_PoI$ensembl_gene_id),
    name4=as.character(TSS_PoI$gene_biotype),
    name5=as.character(TSS_PoI$strand),
    ranges=IRanges(
      start=as.numeric(TSS_PoI$start),
      end=as.numeric(TSS_PoI$end),
      name=TSS_PoI$Peak_ID))
  
  # cat("gr_TSS_PoI_0\n")
  # cat(str(gr_TSS_PoI))
  # cat("\n")
  
  #### Read the Peak to genes results file ----
  
  Linked_PoI<-readRDS(file=opt$Linked_PoI)
  
  
  cat("Linked_PoI_0\n")
  cat(str(Linked_PoI))
  cat("\n")
  cat(sprintf(as.character(names(summary(Linked_PoI$strand)))))
  cat("\n")
  cat(sprintf(as.character(summary(Linked_PoI$strand))))
  cat("\n")
  

  
  gr_Linked_PoI <- GRanges(
    seqnames = as.character(gsub("^chr","",gsub("-.+$","",Linked_PoI$Peak_ID))),
    name2=as.character(Linked_PoI$Symbol),
    ranges=IRanges(
      start=as.integer(gsub("-.+$","",gsub("^[^-]+-","",Linked_PoI$Peak_ID))), 
      end=as.integer(gsub("^[^-]+-[^-]+-","",Linked_PoI$Peak_ID)),
      name=Linked_PoI$Peak_ID))
  
  cat("gr_Linked_PoI$Peak_ID_0\n")
  cat(str(gr_Linked_PoI$Peak_ID))
  cat("\n")
  
 ##### CLASSES ------
  
  
  
  #### Intersect TSS_Peaks to gr_K562_Promoters_in_TSS ----
  
  DEBUG <- 0
  
  m <- findOverlaps(gr_K562_Promoters_in_TSS,gr_TSS_PoI)
  
  if(DEBUG == 1)
  {
    cat("m\n")
    cat(str(m))
    cat("\n")
  }
  
  subjectHits_TSS_PoI<-subjectHits(m)
  
  if(DEBUG == 1)
  {
    cat("subjectHits_TSS_PoI\n")
    cat(str(subjectHits_TSS_PoI))
    cat("\n")
  }
  
  queryHits_Ensembl_promoters_K562_linked_to_genes<-queryHits(m)
  
  if(DEBUG == 1)
  {
    cat("queryHits_Ensembl_promoters_K562_linked_to_genes\n")
    cat(str(queryHits_Ensembl_promoters_K562_linked_to_genes))
    cat("\n")
  }
  
 
  
  Ensembl_promoters_K562_TSS_Peaks <- data.frame(Promoter_chr=as.character(seqnames(gr_K562_Promoters_in_TSS)),
                                                 Promoter_start=as.integer(start(gr_K562_Promoters_in_TSS)),
                                                 Promoter_end=as.integer(end(gr_K562_Promoters_in_TSS)),
                                                  Promoter_Symbol=as.character(gr_K562_Promoters_in_TSS$name2),
                                                          Promoter_ensembl_gene_id=as.character(gr_K562_Promoters_in_TSS$name3),
                                                          Promoter_strand=as.character(gr_K562_Promoters_in_TSS$name4),
                                              feature=as.character(gr_K562_Promoters_in_TSS$name5),
                                              activity=as.character(gr_K562_Promoters_in_TSS$name6),
                                              CLASS=as.character(gr_K562_Promoters_in_TSS$name7),
                                              value_string=names(gr_K562_Promoters_in_TSS), stringsAsFactors = F)
  
  if(DEBUG == 1)
  {
    cat("Ensembl_promoters_K562_TSS_Peaks_0\n")
    cat(str(Ensembl_promoters_K562_TSS_Peaks))
    cat("\n")
  }
  
  Ensembl_promoters_K562_TSS_Peaks_hits<-Ensembl_promoters_K562_TSS_Peaks[queryHits_Ensembl_promoters_K562_linked_to_genes,]
  
  if(DEBUG == 1)
  {
    cat("Ensembl_promoters_K562_TSS_Peaks_hits_0\n")
    cat(str(Ensembl_promoters_K562_TSS_Peaks_hits))
    cat("\n")
  }
  
 
  
  TSS_PoI_df_Ensembl_promoters_K562 <- data.frame(chr=as.character(seqnames(gr_TSS_PoI)),
                                                                        start=as.integer(start(gr_TSS_PoI)),
                                                                        end=as.integer(end(gr_TSS_PoI)),
                                                                        Symbol=as.character(gr_TSS_PoI$name2),
                                                                  ensembl_gene_id=as.character(gr_TSS_PoI$name3),
                                                                  gene_biotype=as.character(gr_TSS_PoI$name4),
                                                                  strand=as.character(gr_TSS_PoI$name5),
                                                            Peak_ID=names(gr_TSS_PoI),
                                                            stringsAsFactors = F)
  
 
  if(DEBUG == 1)
  {
    cat("TSS_PoI_df_Ensembl_promoters_K562_0\n")
    cat(str(TSS_PoI_df_Ensembl_promoters_K562))
    cat("\n")
  }
  
  check<-TSS_PoI_df_Ensembl_promoters_K562[which(TSS_PoI_df_Ensembl_promoters_K562$Symbol%in%tracking_genes),]
  
  
  if(DEBUG == 1)
  {
    cat("check_0\n")
    cat(str(check))
    cat("\n")
    cat(str(unique(check$Symbol)))
    cat("\n")
  }
  
  check_special<-TSS_PoI_df_Ensembl_promoters_K562[which(TSS_PoI_df_Ensembl_promoters_K562$Symbol == 'RUNX1'),]
  
  if(DEBUG == 1)
  {
    cat("check_special_0\n")
    cat(str(check_special))
    cat("\n")
    cat(str(unique(check_special$Symbol)))
    cat("\n")
  }
  
  TSS_PoI_df_Ensembl_promoters_K562_hits<-TSS_PoI_df_Ensembl_promoters_K562[subjectHits_TSS_PoI,]
  
  if(dim(TSS_PoI_df_Ensembl_promoters_K562_hits)[1] >0)
  {
    if(DEBUG == 1)
    {
      cat("TSS_PoI_df_Ensembl_promoters_K562_hits_0\n")
      cat(str(TSS_PoI_df_Ensembl_promoters_K562_hits))
      cat("\n")
    }
    
    TSS_PoI_df_Ensembl_promoters_K562_hits<-cbind(TSS_PoI_df_Ensembl_promoters_K562_hits,Ensembl_promoters_K562_TSS_Peaks_hits)
    
    if(DEBUG == 1)
    {
      cat("TSS_PoI_df_Ensembl_promoters_K562_hits_1\n")
      cat(str(TSS_PoI_df_Ensembl_promoters_K562_hits))
      cat("\n")
    }
    
   

    TSS_PoI_df_Ensembl_promoters_K562_hits_subset<-unique(TSS_PoI_df_Ensembl_promoters_K562_hits[,c(1:8,c(
      which(colnames(TSS_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_ensembl_gene_id'),which(colnames(TSS_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_Symbol'),which(colnames(TSS_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_strand'),
      which(colnames(TSS_PoI_df_Ensembl_promoters_K562_hits) == 'feature'),which(colnames(TSS_PoI_df_Ensembl_promoters_K562_hits) == 'activity'),
      which(colnames(TSS_PoI_df_Ensembl_promoters_K562_hits) == 'CLASS'),which(colnames(TSS_PoI_df_Ensembl_promoters_K562_hits) == 'value_string')))])
    
    if(DEBUG == 1)
    {
      cat("TSS_PoI_df_Ensembl_promoters_K562_hits_subset_0\n")
      cat(str(TSS_PoI_df_Ensembl_promoters_K562_hits_subset))
      cat("\n")
    }
    
    check<-TSS_PoI_df_Ensembl_promoters_K562_hits_subset[which(TSS_PoI_df_Ensembl_promoters_K562_hits_subset$Symbol%in%tracking_genes),]
    
    
    if(DEBUG == 1)
    {
      cat("check_1\n")
      cat(str(check))
      cat("\n")
      cat(str(unique(check$Symbol)))
      cat("\n")
    }
    
    check_special<-TSS_PoI_df_Ensembl_promoters_K562_hits_subset[which(TSS_PoI_df_Ensembl_promoters_K562_hits_subset$Symbol == 'RUNX1'),]
    
    if(DEBUG == 1)
    {
      cat("check_special_1\n")
      cat(str(check_special))
      cat("\n")
      cat(str(unique(check_special$Symbol)))
      cat("\n")
    }
    
    
    #### only keep promoters that match their gene with the gene regulated by the peak ----
    
    indx.match.genes<-which(TSS_PoI_df_Ensembl_promoters_K562_hits$ensembl_gene_id == TSS_PoI_df_Ensembl_promoters_K562_hits$Promoter_ensembl_gene_id &
                              TSS_PoI_df_Ensembl_promoters_K562_hits$chr == TSS_PoI_df_Ensembl_promoters_K562_hits$Promoter_chr &
                              TSS_PoI_df_Ensembl_promoters_K562_hits$strand == TSS_PoI_df_Ensembl_promoters_K562_hits$Promoter_strand)
    
    
    if(DEBUG == 1)
    {
      cat("indx.match.genes_0\n")
      cat(str(indx.match.genes))
      cat("\n")
    }
    
    
    CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits<-TSS_PoI_df_Ensembl_promoters_K562_hits[indx.match.genes,]
    
    
    if(DEBUG == 1)
    {
      cat("CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits_0\n")
      cat(str(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits))
      cat("\n")
    }
    
    indx.dep<-c(which(colnames(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_Symbol'),which(colnames(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_ensembl_gene_id'),which(colnames(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_chr'),
                which(colnames(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_strand'),which(colnames(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_start'),which(colnames(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_end'))
    
    CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits<-unique(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits[,-indx.dep])
    
    
    # if(DEBUG == 1)
    # {
      cat("CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits_1\n")
      cat(str(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits))
      cat("\n")
    # }
    
    check<-CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits[which(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits$Symbol%in%tracking_genes),]
    
    # if(DEBUG == 1)
    # {
      cat("check_2\n")
      cat(str(check))
      cat("\n")
      cat(str(unique(check$Symbol)))
      cat("\n")
    # }
    
    
    
    
  }#dim(TSS_PoI_df_Ensembl_promoters_K562_hits)[1] >0
  
  
  #### Intersect Linked_Peaks to gr_K562_Promoters_in_TSS ----
  
  DEBUG <- 0
  
  m <- findOverlaps(gr_K562_Promoters_in_TSS,gr_Linked_PoI)
  
  if(DEBUG == 1)
  {
    cat("m\n")
    cat(str(m))
    cat("\n")
  }
  
  subjectHits_ALL_Linked_PoI<-subjectHits(m)
  
  if(DEBUG == 1)
  {
    cat("subjectHits_ALL_Linked_PoI\n")
    cat(str(subjectHits_ALL_Linked_PoI))
    cat("\n")
  }
  
  queryHits_Ensembl_promoters_K562_linked_to_genes<-queryHits(m)
  
  if(DEBUG == 1)
  {
    cat("queryHits_Ensembl_promoters_K562_linked_to_genes\n")
    cat(str(queryHits_Ensembl_promoters_K562_linked_to_genes))
    cat("\n")
  }
  
  
  
  Ensembl_promoters_K562_TSS_Peaks <- data.frame(Promoter_chr=as.character(seqnames(gr_K562_Promoters_in_TSS)),
                                                 Promoter_start=as.integer(start(gr_K562_Promoters_in_TSS)),
                                                 Promoter_end=as.integer(end(gr_K562_Promoters_in_TSS)),
                                                 Promoter_Symbol=as.character(gr_K562_Promoters_in_TSS$name2),
                                                 ensembl_gene_id=as.character(gr_K562_Promoters_in_TSS$name3),
                                                 strand=as.character(gr_K562_Promoters_in_TSS$name4),
                                                 gene_biotype=as.character(gr_K562_Promoters_in_TSS$name8),
                                                 feature=as.character(gr_K562_Promoters_in_TSS$name5),
                                                 activity=as.character(gr_K562_Promoters_in_TSS$name6),
                                                 CLASS=as.character(gr_K562_Promoters_in_TSS$name7),
                                                 value_string=names(gr_K562_Promoters_in_TSS), stringsAsFactors = F)
  
  if(DEBUG == 1)
  {
    cat("Ensembl_promoters_K562_TSS_Peaks_0\n")
    cat(str(Ensembl_promoters_K562_TSS_Peaks))
    cat("\n")
  }
  
  Ensembl_promoters_K562_TSS_Peaks_hits<-Ensembl_promoters_K562_TSS_Peaks[queryHits_Ensembl_promoters_K562_linked_to_genes,]
  
  if(DEBUG == 1)
  {
    cat("Ensembl_promoters_K562_TSS_Peaks_hits_0\n")
    cat(str(Ensembl_promoters_K562_TSS_Peaks_hits))
    cat("\n")
  }
  
  
  
  ALL_Linked_PoI_df_Ensembl_promoters_K562 <- data.frame(chr=as.character(seqnames(gr_Linked_PoI)),
                                                  start=as.integer(start(gr_Linked_PoI)),
                                                  end=as.integer(end(gr_Linked_PoI)),
                                                  Symbol=as.character(gr_Linked_PoI$name2),
                                                  Peak_ID=names(gr_Linked_PoI),
                                                  stringsAsFactors = F)
  
  
  if(DEBUG == 1)
  {
    cat("ALL_Linked_PoI_df_Ensembl_promoters_K562_0\n")
    cat(str(ALL_Linked_PoI_df_Ensembl_promoters_K562))
    cat("\n")
  }
  
  check<-ALL_Linked_PoI_df_Ensembl_promoters_K562[which(ALL_Linked_PoI_df_Ensembl_promoters_K562$Symbol%in%tracking_genes),]
  
  
  if(DEBUG == 1)
  {
    cat("check_0\n")
    cat(str(check))
    cat("\n")
    cat(str(unique(check$Symbol)))
    cat("\n")
  }
  
  check_special<-ALL_Linked_PoI_df_Ensembl_promoters_K562[which(ALL_Linked_PoI_df_Ensembl_promoters_K562$Symbol == 'RUNX1'),]
  
  if(DEBUG == 1)
  {
    cat("check_special_0\n")
    cat(str(check_special))
    cat("\n")
    cat(str(unique(check_special$Symbol)))
    cat("\n")
  }
  
  ALL_Linked_PoI_df_Ensembl_promoters_K562_hits<-ALL_Linked_PoI_df_Ensembl_promoters_K562[subjectHits_ALL_Linked_PoI,]
  
  if(dim(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits)[1] >0)
  {
    if(DEBUG == 1)
    {
      cat("ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_0\n")
      cat(str(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits))
      cat("\n")
    }
    
    ALL_Linked_PoI_df_Ensembl_promoters_K562_hits<-cbind(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits,Ensembl_promoters_K562_TSS_Peaks_hits)
    
    if(DEBUG == 1)
    {
      cat("ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_1\n")
      cat(str(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits))
      cat("\n")
    }
    
    
    
    ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_subset<-unique(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits[,c(1:5,c(
      which(colnames(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_ensembl_gene_id'),which(colnames(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_Symbol'),which(colnames(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_strand'),
      which(colnames(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'feature'),which(colnames(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'activity'),
      which(colnames(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'CLASS'),which(colnames(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'value_string')))])
    
    if(DEBUG == 1)
    {
      cat("ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_subset_0\n")
      cat(str(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_subset))
      cat("\n")
    }
    
    check<-ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_subset[which(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_subset$Symbol%in%tracking_genes),]
    
    
    if(DEBUG == 1)
    {
      cat("check_1\n")
      cat(str(check))
      cat("\n")
      cat(str(unique(check$Symbol)))
      cat("\n")
    }
    
    check_special<-ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_subset[which(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_subset$Symbol == 'RUNX1'),]
    
    if(DEBUG == 1)
    {
      cat("check_special_1\n")
      cat(str(check_special))
      cat("\n")
      cat(str(unique(check_special$Symbol)))
      cat("\n")
    }
    
    
    #### only keep promoters that match their gene with the gene regulated by the peak ----
    
    indx.match.genes<-which(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits$Symbol == ALL_Linked_PoI_df_Ensembl_promoters_K562_hits$Promoter_Symbol &
                              ALL_Linked_PoI_df_Ensembl_promoters_K562_hits$chr == ALL_Linked_PoI_df_Ensembl_promoters_K562_hits$Promoter_chr)
    
    
    if(DEBUG == 1)
    {
      cat("indx.match.genes_0\n")
      cat(str(indx.match.genes))
      cat("\n")
    }
    
    
    CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits<-ALL_Linked_PoI_df_Ensembl_promoters_K562_hits[indx.match.genes,]
    
    
    if(DEBUG == 1)
    {
      cat("CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_0\n")
      cat(str(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits))
      cat("\n")
    }
    
    indx.dep<-c(which(colnames(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_Symbol'),which(colnames(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_ensembl_gene_id'),which(colnames(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_chr'),
                which(colnames(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_strand'),which(colnames(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_start'),which(colnames(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits) == 'Promoter_end'))
    
    CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits<-unique(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits[,-indx.dep])
    
    
    # if(DEBUG == 1)
    # {
      cat("CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_1\n")
      cat(str(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits))
      cat("\n")
    # }
    
    check<-CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits[which(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits$Symbol%in%tracking_genes),]
    
    # if(DEBUG == 1)
    # {
      cat("check_2\n")
      cat(str(check))
      cat("\n")
      cat(str(unique(check$Symbol)))
      cat("\n")
    # }
    
    
    
    
  }#dim(ALL_Linked_PoI_df_Ensembl_promoters_K562_hits)[1] >0
  
  
  ###### Merge all the Peaks assigned to genes ------
  
  
  if(dim(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits)[1] >0 & dim(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits)[1] >0)
  {
    CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits$Peak_CLASS<-'TSS_Peak'
    
    cat("CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits_FINAL\n")
    cat(str(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits))
    cat("\n")
    
    
    CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_NR<-CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits[-which(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits$Peak_ID%in%CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits$Peak_ID),]
    
    if(dim(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_NR)[1] >0)
    {
      CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_NR$Peak_CLASS<-'Linked_Peak'
      
      cat("CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_NR_0\n")
      cat(str(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_NR))
      cat("\n")
      
      
      DEF_promoter_peaks<-rbind(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits,
                                CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_NR)
      
      
      
    }else{
      
      DEF_promoter_peaks<-CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits
      
    }#dim(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits_NR)[1] >0
    
    DEBUG<-1
    
    cat("DEF_promoter_peaks_0\n")
    cat(str(DEF_promoter_peaks))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(DEF_promoter_peaks$Peak_CLASS))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(DEF_promoter_peaks$Peak_CLASS)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(DEF_promoter_peaks$CLASS))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(DEF_promoter_peaks$CLASS)))))
    cat("\n")
    
    # c('Peak_ID',"chr","start",
    #   "end","Symbol",
    #   "gene_biotype",
    #   "strand","ensembl_gene_id",
    #   "Peak_CLASS")
    
    DEF_promoter_peaks.dt<-data.table(DEF_promoter_peaks, key=c('Peak_ID',"chr","start",
                                                                                "end","Peak_CLASS"))
    
    DEF_promoter_peaks_collapsed<-as.data.frame(DEF_promoter_peaks.dt[,.(variable_string=paste(unique(sort(CLASS)), collapse=';'),
                                           value_string=paste(value_string, collapse=';'),
                                           Symbol_string=paste(unique(sort(Symbol)), collapse=';')),by=key(DEF_promoter_peaks.dt)], 
                                           stringsAsFactors=F)
    
    
    DEF_promoter_peaks_collapsed$variable_string<-gsub("_",'|',DEF_promoter_peaks_collapsed$variable_string)
    
    if(DEBUG == 1)
    {
      cat("DEF_promoter_peaks_collapsed_0\n")
      cat(str(DEF_promoter_peaks_collapsed))
      cat("\n")
      cat(sprintf(paste(as.character(names(summary(as.factor(DEF_promoter_peaks_collapsed$variable_string)))), collapse="BOOM")))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(DEF_promoter_peaks_collapsed$variable_string)))))
      cat("\n")
    }
    
    ### classification of promoter activity
    
    DEF_promoter_peaks_collapsed$CLASS<-NA
    
    Promoter_ACTIVE_CLASS<-c('Promoter|ACTIVE','Promoter|ACTIVE;Promoter|INACTIVE','Promoter|ACTIVE;Promoter|POISED','Promoter|ACTIVE;Promoter|REPRESSED',
                             'Promoter|ACTIVE;Promoter|INACTIVE;Promoter|POISED','Promoter|ACTIVE;Promoter|INACTIVE;Promoter|REPRESSED',
                             'Promoter|ACTIVE;Promoter|INACTIVE;Promoter|POISED;Promoter|REPRESSED')
    Promoter_INACTIVE_CLASS<-c('Promoter|INACTIVE','Promoter|INACTIVE;Promoter|POISED','Promoter|INACTIVE;Promoter|POISED;Promoter|REPRESSED','Promoter|INACTIVE;Promoter|REPRESSED')
    Promoter_POISED_CLASS<-c('Promoter|POISED','Promoter|POISED;Promoter|REPRESSED')
    Promoter_REPRESSED_CLASS<-c('Promoter|REPRESSED')
    
    DEF_promoter_peaks_collapsed$CLASS[which(DEF_promoter_peaks_collapsed$variable_string%in%Promoter_ACTIVE_CLASS)]<-'Promoter_ACTIVE'
    DEF_promoter_peaks_collapsed$CLASS[which(DEF_promoter_peaks_collapsed$variable_string%in%Promoter_INACTIVE_CLASS)]<-'Promoter_INACTIVE'
    DEF_promoter_peaks_collapsed$CLASS[which(DEF_promoter_peaks_collapsed$variable_string%in%Promoter_POISED_CLASS)]<-'Promoter_POISED'
    DEF_promoter_peaks_collapsed$CLASS[which(DEF_promoter_peaks_collapsed$variable_string%in%Promoter_REPRESSED_CLASS)]<-'Promoter_REPRESSED'
    
    DEF_promoter_peaks_collapsed$feature<-gsub("_.+$","",DEF_promoter_peaks_collapsed$CLASS)
    DEF_promoter_peaks_collapsed$activity<-gsub("^[^_]+_","",DEF_promoter_peaks_collapsed$CLASS)
    
    if(DEBUG == 1)
    {
      cat("DEF_promoter_peaks_collapsed_1\n")
      cat(str(DEF_promoter_peaks_collapsed))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(DEF_promoter_peaks_collapsed$CLASS))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(DEF_promoter_peaks_collapsed$CLASS)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(DEF_promoter_peaks_collapsed$feature))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(DEF_promoter_peaks_collapsed$feature)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(DEF_promoter_peaks_collapsed$activity))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(DEF_promoter_peaks_collapsed$activity)))))
      cat("\n")
    }
    
    check.NA<-DEF_promoter_peaks_collapsed[is.na(DEF_promoter_peaks_collapsed$CLASS),]
    
    # if(DEBUG == 1)
    # {
      cat("check.NA_0\n")
      cat(str(check.NA))
      cat("\n")
      cat(sprintf(paste(as.character(names(summary(as.factor(check.NA$variable_string)))), collapse="BOOM")))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(check.NA$variable_string)))))
      cat("\n")
    # }
    
    
    

    DEF_promoter_peaks_collapsed$Peak_CLASS<-factor(DEF_promoter_peaks_collapsed$Peak_CLASS,
                                          levels=c('Non_TSS_Peak','TSS_Peak'),
                                          ordered=T)
    
    
    DEF_promoter_peaks_collapsed$CLASS<-factor(DEF_promoter_peaks_collapsed$CLASS,
                                          levels=c('Promoter_ACTIVE','Promoter_POISED','Promoter_REPRESSED','Promoter_INACTIVE'),
                                          ordered=T)
    
    cat("DEF_promoter_peaks_collapsed_1\n")
    cat(str(DEF_promoter_peaks_collapsed))
    cat("\n")
    cat(sprintf(as.character(names(summary(DEF_promoter_peaks_collapsed$Peak_CLASS)))))
    cat("\n")
    cat(sprintf(as.character(summary(DEF_promoter_peaks_collapsed$Peak_CLASS))))
    cat("\n")
    cat(sprintf(as.character(names(summary(DEF_promoter_peaks_collapsed$CLASS)))))
    cat("\n")
    cat(sprintf(as.character(summary(DEF_promoter_peaks_collapsed$CLASS))))
    cat("\n")
    
    
    check<-DEF_promoter_peaks_collapsed[grep(paste(tracking_genes,
                                                   collapse='|'),DEF_promoter_peaks_collapsed$Symbol_string),]
    
    # if(DEBUG == 1)
    # {
    cat("check_3\n")
    cat(str(check))
    cat("\n")
    cat(str(unique(check$Symbol)))
    cat("\n")
    # }
    
    setwd(out)
    
    saveRDS(DEF_promoter_peaks_collapsed,file='PoI_concordant_promoters.rds')
    write.table(DEF_promoter_peaks_collapsed,file='PoI_concordant_promoters.tsv',sep="\t",quote=F,row.names = F)
    
  }#dim(CONCORDANT_TSS_PoI_df_Ensembl_promoters_K562_hits)[1] >0 & dim(CONCORDANT_ALL_Linked_PoI_df_Ensembl_promoters_K562_hits)[1] >0
  
}




printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}



#### main script ----



main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--tracking_genes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--K562_Regulatory_Build"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--K562_Promoters_in_TSS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TSS_PoI"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Linked_PoI"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
 
  Assign_TSS_and_linked_peaks(opt)
 

  
}


###########################################################################

system.time( main() )