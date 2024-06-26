
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


classify_the_REST_promoters_and_non_promoters = function(option_list)
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
  

  
  gr_K562_Regulatory_Build <- GRanges(
    seqnames = as.character(K562_Regulatory_Build$seqid),
    name2=as.character(K562_Regulatory_Build$feature_type),
    name3=as.character(K562_Regulatory_Build$activity),
    ranges=IRanges(
      start=as.numeric(K562_Regulatory_Build$start),
      end=as.numeric(K562_Regulatory_Build$end),
      name=K562_Regulatory_Build$regulatory_feature_stable_id))
  
  #### Read PoI_promoters file ----
  
  PoI_promoters<-readRDS(file=opt$PoI_promoters)
  
  
  cat("PoI_promoters_0\n")
  cat(str(PoI_promoters))
  cat("\n")
  cat(sprintf(as.character(names(summary(PoI_promoters$feature)))))
  cat("\n")
  cat(sprintf(as.character(summary(PoI_promoters$feature))))
  cat("\n")
  
  PoI_promoters<-unique(PoI_promoters[,-c(which(colnames(PoI_promoters) == 'variable_string'),
                                   which(colnames(PoI_promoters) == 'CLASS'))])
  
  cat("PoI_promoters_1\n")
  cat(str(PoI_promoters))
  cat("\n")
  cat(sprintf(as.character(names(summary(PoI_promoters$feature)))))
  cat("\n")
  cat(sprintf(as.character(summary(PoI_promoters$feature))))
  cat("\n")
  
  #### Read Non_promoters_PoI file ----
  
  Non_promoters_PoI<-readRDS(file=opt$Non_promoters_PoI)
  
  cat("Non_promoters_PoI_1\n")
  cat(str(Non_promoters_PoI))
  cat("\n")
  cat(sprintf(as.character(names(summary(Non_promoters_PoI$feature)))))
  cat("\n")
  cat(sprintf(as.character(summary(Non_promoters_PoI$feature))))
  cat("\n")
  
  #### Read the ALL_PoI file ----
  
  ALL_PoI<-readRDS(file=opt$ALL_PoI)
  
  cat("ALL_PoI_0\n")
  cat(str(ALL_PoI))
  cat("\n")
  cat(str(unique(ALL_PoI$Peak_ID)))
  cat("\n")
  
 
  ALL_PoI_REST<-ALL_PoI[-which(ALL_PoI$Peak_ID%in%PoI_promoters$Peak_ID),]
  
  cat("ALL_PoI_REST_0\n")
  cat(str(ALL_PoI_REST))
  cat("\n")
  cat(str(unique(ALL_PoI_REST$Peak_ID)))
  cat("\n")
  
  ALL_PoI_REST<-ALL_PoI_REST[-which(ALL_PoI_REST$Peak_ID%in%Non_promoters_PoI$Peak_ID),]
  
  cat("ALL_PoI_REST_1\n")
  cat(str(ALL_PoI_REST))
  cat("\n")
  cat(str(unique(ALL_PoI_REST$Peak_ID)))
  cat("\n")
 
  
  
  gr_ALL_PoI_REST <- GRanges(
    seqnames = as.character(gsub("^chr","",ALL_PoI_REST$chr)),
    ranges=IRanges(
      start=as.numeric(ALL_PoI_REST$start),
      end=as.numeric(ALL_PoI_REST$end),
      name=ALL_PoI_REST$Peak_ID))
  
  # cat("gr_ALL_PoI_REST_0\n")
  # cat(str(gr_ALL_PoI_REST))
  # cat("\n")
  
 
  #### CLasses of tags ------------------------
  
  
  
  
  Promoter_ACTIVE_CLASS<-c('CTCF Binding Site|POISED;CTCF Binding Site|REPRESSED;Promoter|ACTIVE',
                           'Enhancer|INACTIVE;Promoter|ACTIVE','Promoter|ACTIVE','Promoter|ACTIVE;Promoter|INACTIVE','Promoter|ACTIVE;Promoter|POISED','Promoter|ACTIVE;Promoter|REPRESSED',
                           'Promoter|ACTIVE;Promoter|INACTIVE;Promoter|POISED','Promoter|ACTIVE;Promoter|INACTIVE;Promoter|REPRESSED',
                           'Promoter|ACTIVE;Promoter|INACTIVE;Promoter|POISED;Promoter|REPRESSED',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;Promoter|ACTIVE',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;CTCF Binding Site|REPRESSED;Promoter|ACTIVE',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Promoter|ACTIVE',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Promoter|ACTIVE;Promoter|INACTIVE',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|REPRESSED;Promoter|ACTIVE',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Promoter|ACTIVE',
                           'CTCF Binding Site|ACTIVE;Enhancer|REPRESSED;Promoter|ACTIVE',
                           'CTCF Binding Site|ACTIVE;Promoter|ACTIVE',
                           'CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;Promoter|ACTIVE',
                           'CTCF Binding Site|INACTIVE;Promoter|ACTIVE',
                           'CTCF Binding Site|INACTIVE;Promoter|ACTIVE;Promoter|INACTIVE',
                           'CTCF Binding Site|POISED;Promoter|ACTIVE',
                           'CTCF Binding Site|POISED;Promoter|ACTIVE;Promoter|INACTIVE',
                           'CTCF Binding Site|POISED;Promoter|ACTIVE;Promoter|INACTIVE;Promoter|POISED',
                           'CTCF Binding Site|REPRESSED;Promoter|ACTIVE',
                           'CTCF Binding Site|ACTIVE;Promoter|ACTIVE;Promoter|INACTIVE',
                           'CTCF Binding Site|ACTIVE;Promoter|ACTIVE;Promoter|POISED')
  Promoter_INACTIVE_CLASS<-c('Enhancer|ACTIVE;Promoter|INACTIVE',
                             'Enhancer|INACTIVE;Promoter|INACTIVE',
                             'Promoter|INACTIVE;TF binding|ACTIVE','Promoter|INACTIVE',
                             'Promoter|INACTIVE;Promoter|POISED','Promoter|INACTIVE;Promoter|POISED;Promoter|REPRESSED',
                             'Promoter|INACTIVE;Promoter|REPRESSED',
                             'CTCF Binding Site|POISED;Enhancer|ACTIVE;Promoter|INACTIVE',
                             'CTCF Binding Site|POISED;Enhancer|POISED;Open chromatin|POISED;Promoter|INACTIVE;TF binding|POISED',
                             'CTCF Binding Site|POISED;Promoter|INACTIVE',
                             'CTCF Binding Site|POISED;Promoter|INACTIVE;Promoter|REPRESSED',
                             'CTCF Binding Site|REPRESSED;Promoter|INACTIVE',
                             'CTCF Binding Site|REPRESSED;Promoter|INACTIVE;Promoter|REPRESSED',
                             'CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Promoter|INACTIVE',
                             'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Promoter|INACTIVE',
                             'CTCF Binding Site|ACTIVE;Promoter|INACTIVE',
                             'CTCF Binding Site|ACTIVE;Promoter|INACTIVE;Promoter|REPRESSED',
                             'CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;Promoter|INACTIVE',
                             'CTCF Binding Site|INACTIVE;CTCF Binding Site|REPRESSED;Promoter|INACTIVE',
                             'CTCF Binding Site|INACTIVE;Promoter|INACTIVE',
                             'CTCF Binding Site|INACTIVE;Promoter|INACTIVE;Promoter|POISED',
                             'CTCF Binding Site|INACTIVE;Promoter|INACTIVE;Promoter|REPRESSED',
                             'CTCF Binding Site|POISED;Enhancer|REPRESSED;Promoter|INACTIVE')
  Promoter_POISED_CLASS<-c('Promoter|POISED','Promoter|POISED;Promoter|REPRESSED')
  Promoter_REPRESSED_CLASS<-c('CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Promoter|POISED',
                              'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE;Promoter|POISED',
                              'CTCF Binding Site|ACTIVE;Promoter|INACTIVE;Promoter|POISED',
                              'CTCF Binding Site|ACTIVE;Promoter|POISED',
                              'CTCF Binding Site|INACTIVE;Promoter|POISED',
                              'CTCF Binding Site|REPRESSED;Promoter|POISED',
                              'CTCF Binding Site|POISED;Promoter|POISED',
                              'CTCF Binding Site|POISED;Promoter|POISED;Promoter|REPRESSED','Promoter|REPRESSED','Promoter|REPRESSED;TF binding|ACTIVE',
                              'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Promoter|REPRESSED',
                              'CTCF Binding Site|ACTIVE;Promoter|REPRESSED',
                              'CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;Promoter|REPRESSED',
                              'CTCF Binding Site|INACTIVE;Enhancer|INACTIVE;Promoter|REPRESSED',
                              'CTCF Binding Site|REPRESSED;Promoter|REPRESSED',
                              'CTCF Binding Site|POISED;Promoter|REPRESSED',
                              'CTCF Binding Site|INACTIVE;Promoter|POISED;Promoter|REPRESSED',
                              'CTCF Binding Site|INACTIVE;Promoter|REPRESSED',
                              'CTCF Binding Site|POISED;CTCF Binding Site|REPRESSED;Promoter|REPRESSED')
  
  Enhancer_ACTIVE_CLASS<-c('Enhancer|ACTIVE;Open chromatin|REPRESSED',
                           'Enhancer|ACTIVE;TF binding|REPRESSED','Enhancer|ACTIVE','Enhancer|ACTIVE;Enhancer|INACTIVE','Enhancer|ACTIVE;TF binding|ACTIVE',
                           'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE','CTCF Binding Site|INACTIVE;Enhancer|ACTIVE',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Enhancer|ACTIVE',
                           'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE;Open chromatin|ACTIVE;TF binding|ACTIVE',
                           'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE;Open chromatin|INACTIVE',
                           'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE;TF binding|ACTIVE',
                           'CTCF Binding Site|INACTIVE;Enhancer|ACTIVE;Enhancer|INACTIVE',
                           'CTCF Binding Site|INACTIVE;Enhancer|ACTIVE;Open chromatin|INACTIVE',
                           'CTCF Binding Site|INACTIVE;Enhancer|ACTIVE;Open chromatin|INACTIVE;TF binding|ACTIVE',
                           'CTCF Binding Site|INACTIVE;Enhancer|ACTIVE;TF binding|ACTIVE',
                           'CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;Enhancer|ACTIVE;Enhancer|POISED',
                           'CTCF Binding Site|POISED;Enhancer|ACTIVE',
                           'CTCF Binding Site|POISED;Enhancer|ACTIVE;Enhancer|POISED',
                           'Enhancer|ACTIVE;Open chromatin|ACTIVE;TF binding|ACTIVE',
                           'Enhancer|ACTIVE;Open chromatin|INACTIVE',
                           'CTCF Binding Site|REPRESSED;Enhancer|ACTIVE;Open chromatin|REPRESSED',
                           'CTCF Binding Site|POISED;Open chromatin|POISED;TF binding|ACTIVE;TF binding|POISED',
                           'CTCF Binding Site|POISED;Open chromatin|REPRESSED;TF binding|ACTIVE',
                           'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE;Enhancer|INACTIVE',
                           'Enhancer|ACTIVE;TF binding|INACTIVE','Enhancer|ACTIVE;Open chromatin|ACTIVE',
                           'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE;Enhancer|INACTIVE;TF binding|ACTIVE',
                           'Enhancer|ACTIVE;Enhancer|INACTIVE;TF binding|ACTIVE',
                           'Enhancer|ACTIVE;Open chromatin|INACTIVE;TF binding|ACTIVE',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Enhancer|ACTIVE;TF binding|ACTIVE',
                           'Enhancer|ACTIVE;Enhancer|INACTIVE;Open chromatin|INACTIVE',
                           'Enhancer|ACTIVE;Enhancer|POISED',
                           'Enhancer|ACTIVE;Open chromatin|ACTIVE;Open chromatin|INACTIVE',
                           'Enhancer|ACTIVE;Enhancer|REPRESSED',
                           'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE;Open chromatin|ACTIVE',
                           'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE;TF binding|INACTIVE','CTCF Binding Site|POISED;Enhancer|ACTIVE;Enhancer|REPRESSED',
                           'CTCF Binding Site|POISED;Enhancer|ACTIVE;TF binding|ACTIVE',
                           'CTCF Binding Site|INACTIVE;Enhancer|ACTIVE;Enhancer|INACTIVE;TF binding|ACTIVE',
                           'CTCF Binding Site|INACTIVE;Enhancer|ACTIVE;TF binding|INACTIVE',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|REPRESSED;Enhancer|ACTIVE;Enhancer|REPRESSED;Open chromatin|INACTIVE',
                           'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE;Open chromatin|INACTIVE;TF binding|ACTIVE',
                           'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE;Open chromatin|INACTIVE;TF binding|INACTIVE',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Enhancer|ACTIVE;Enhancer|INACTIVE;Open chromatin|INACTIVE',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Enhancer|ACTIVE',
                           'CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;Enhancer|ACTIVE')
  
  
  

  Enhancer_INACTIVE_CLASS<-c('Enhancer|INACTIVE;Enhancer|REPRESSED;Open chromatin|INACTIVE',
                             'Enhancer|INACTIVE;Enhancer|REPRESSED;Open chromatin|REPRESSED',
                             'Enhancer|INACTIVE;Open chromatin|REPRESSED',
                             'Enhancer|INACTIVE;TF binding|REPRESSED','CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;Enhancer|INACTIVE',
                             'CTCF Binding Site|ACTIVE;Enhancer|INACTIVE;Enhancer|POISED','CTCF Binding Site|POISED;Enhancer|INACTIVE;TF binding|INACTIVE','CTCF Binding Site|POISED;Enhancer|INACTIVE;TF binding|POISED',
                             'CTCF Binding Site|POISED;Enhancer|INACTIVE',
                             'CTCF Binding Site|POISED;Enhancer|INACTIVE;Open chromatin|INACTIVE',
                             'Enhancer|INACTIVE','Enhancer|INACTIVE;Enhancer|REPRESSED','Enhancer|INACTIVE;Open chromatin|ACTIVE','Enhancer|INACTIVE;Open chromatin|INACTIVE',
                             'Enhancer|INACTIVE;TF binding|ACTIVE','Enhancer|INACTIVE;TF binding|INACTIVE',
                             'CTCF Binding Site|INACTIVE;Enhancer|INACTIVE','CTCF Binding Site|ACTIVE;Enhancer|INACTIVE',
                             'CTCF Binding Site|ACTIVE;Enhancer|INACTIVE;Open chromatin|INACTIVE','CTCF Binding Site|ACTIVE;Enhancer|INACTIVE;TF binding|ACTIVE',
                             "Enhancer|INACTIVE;Enhancer|POISED",
                             'CTCF Binding Site|INACTIVE;Enhancer|INACTIVE;Enhancer|REPRESSED;Open chromatin|INACTIVE',
                             'CTCF Binding Site|INACTIVE;Enhancer|INACTIVE;Open chromatin|INACTIVE',
                             'CTCF Binding Site|INACTIVE;Enhancer|INACTIVE;TF binding|ACTIVE',
                             'CTCF Binding Site|ACTIVE;Enhancer|INACTIVE;Open chromatin|POISED',
                             'CTCF Binding Site|ACTIVE;Enhancer|INACTIVE;Open chromatin|REPRESSED',
                             'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Enhancer|INACTIVE',
                             'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Enhancer|INACTIVE;Open chromatin|INACTIVE',
                             'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Enhancer|INACTIVE;TF binding|ACTIVE',
                             'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Enhancer|INACTIVE;TF binding|INACTIVE',
                             'CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Enhancer|INACTIVE',
                             'CTCF Binding Site|ACTIVE;Enhancer|INACTIVE;TF binding|INACTIVE',
                             'CTCF Binding Site|ACTIVE;Enhancer|INACTIVE;Open chromatin|ACTIVE',
                             'CTCF Binding Site|INACTIVE;Enhancer|INACTIVE;TF binding|INACTIVE',
                             'CTCF Binding Site|ACTIVE;Enhancer|INACTIVE;Open chromatin|INACTIVE;TF binding|ACTIVE',
                             'CTCF Binding Site|ACTIVE;Enhancer|INACTIVE;Enhancer|REPRESSED',
                             'CTCF Binding Site|POISED;Enhancer|INACTIVE;Enhancer|REPRESSED',
                             'CTCF Binding Site|POISED;Enhancer|INACTIVE;TF binding|ACTIVE',
                             'Enhancer|INACTIVE;Open chromatin|INACTIVE;TF binding|ACTIVE',
                             'Enhancer|INACTIVE;TF binding|POISED',
                             'CTCF Binding Site|REPRESSED;Enhancer|INACTIVE;Enhancer|REPRESSED',
                             'CTCF Binding Site|REPRESSED;Enhancer|INACTIVE;Open chromatin|REPRESSED',
                             'CTCF Binding Site|ACTIVE;Enhancer|ACTIVE;Enhancer|INACTIVE;Open chromatin|INACTIVE',
                             'CTCF Binding Site|REPRESSED;Enhancer|INACTIVE',
                             'CTCF Binding Site|POISED;Enhancer|INACTIVE;Open chromatin|REPRESSED')
  
  
  
  
  Enhancer_POISED_CLASS<-c('Enhancer|POISED;TF binding|ACTIVE',
                           'Enhancer|POISED;TF binding|INACTIVE','CTCF Binding Site|POISED;Enhancer|POISED;TF binding|REPRESSED',
                           'CTCF Binding Site|REPRESSED;Enhancer|POISED;Enhancer|REPRESSED','Enhancer|POISED;Open chromatin|POISED','Enhancer|POISED','Enhancer|POISED;TF binding|POISED','CTCF Binding Site|POISED;Enhancer|POISED','CTCF Binding Site|INACTIVE;Enhancer|POISED',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Enhancer|POISED','CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Enhancer|POISED',
                           'CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;Enhancer|POISED',
                           'CTCF Binding Site|POISED;Enhancer|POISED;Enhancer|REPRESSED',
                           'CTCF Binding Site|POISED;Enhancer|POISED;Open chromatin|POISED',
                           'CTCF Binding Site|POISED;Enhancer|POISED;TF binding|INACTIVE',
                           'CTCF Binding Site|POISED;Enhancer|POISED;TF binding|POISED',
                           'Enhancer|POISED;Open chromatin|INACTIVE','CTCF Binding Site|POISED;Enhancer|POISED;Open chromatin|POISED;TF binding|POISED',
                           'CTCF Binding Site|ACTIVE;Enhancer|POISED',
                           'CTCF Binding Site|REPRESSED;Enhancer|POISED',
                           'Enhancer|POISED;Open chromatin|REPRESSED','CTCF Binding Site|POISED;Enhancer|INACTIVE;Enhancer|POISED',
                           'CTCF Binding Site|POISED;CTCF Binding Site|REPRESSED;Enhancer|POISED',
                           'CTCF Binding Site|INACTIVE;Enhancer|INACTIVE;Enhancer|POISED',
                           'CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;Enhancer|POISED;TF binding|INACTIVE')
  
  Enhancer_REPRESSED_CLASS<-c('Enhancer|REPRESSED;TF binding|ACTIVE','CTCF Binding Site|POISED;CTCF Binding Site|REPRESSED;Enhancer|POISED;Enhancer|REPRESSED',
                              'CTCF Binding Site|ACTIVE;Enhancer|REPRESSED;TF binding|REPRESSED',
                              'CTCF Binding Site|REPRESSED;Enhancer|REPRESSED;TF binding|POISED',
                              'Enhancer|REPRESSED;Open chromatin|POISED','Enhancer|REPRESSED','Enhancer|REPRESSED;Open chromatin|REPRESSED','CTCF Binding Site|INACTIVE;Enhancer|REPRESSED','CTCF Binding Site|POISED;Enhancer|REPRESSED',
                              'CTCF Binding Site|REPRESSED;Enhancer|REPRESSED','CTCF Binding Site|REPRESSED;Enhancer|REPRESSED;Open chromatin|REPRESSED',
                              'CTCF Binding Site|REPRESSED;Enhancer|REPRESSED;TF binding|REPRESSED','CTCF Binding Site|ACTIVE;CTCF Binding Site|REPRESSED;Enhancer|REPRESSED',
                              'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Enhancer|REPRESSED','CTCF Binding Site|ACTIVE;Enhancer|REPRESSED',
                              'CTCF Binding Site|POISED;Enhancer|REPRESSED;Open chromatin|INACTIVE',
                              'CTCF Binding Site|POISED;Enhancer|REPRESSED;TF binding|ACTIVE',
                              'Enhancer|POISED;Enhancer|REPRESSED',
                              'Enhancer|REPRESSED;Open chromatin|INACTIVE',
                              'Enhancer|REPRESSED;TF binding|INACTIVE',
                              'Enhancer|REPRESSED;TF binding|POISED',
                              'Enhancer|REPRESSED;TF binding|REPRESSED',
                              'CTCF Binding Site|INACTIVE;CTCF Binding Site|REPRESSED;Enhancer|REPRESSED;Open chromatin|INACTIVE',
                              'CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Enhancer|REPRESSED;TF binding|ACTIVE',
                              'CTCF Binding Site|ACTIVE;Enhancer|REPRESSED;Open chromatin|REPRESSED',
                              'CTCF Binding Site|POISED;CTCF Binding Site|REPRESSED;Enhancer|REPRESSED',
                              'CTCF Binding Site|REPRESSED;Enhancer|REPRESSED;Open chromatin|INACTIVE',
                              'CTCF Binding Site|ACTIVE;Enhancer|REPRESSED;Open chromatin|INACTIVE',
                              'CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;Enhancer|REPRESSED',
                              'CTCF Binding Site|INACTIVE;CTCF Binding Site|REPRESSED;Enhancer|REPRESSED',
                              'CTCF Binding Site|POISED;Enhancer|REPRESSED;Open chromatin|POISED',
                              'CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Enhancer|REPRESSED',
                              'CTCF Binding Site|POISED;Enhancer|REPRESSED;Open chromatin|REPRESSED',
                              'CTCF Binding Site|POISED;Enhancer|REPRESSED;TF binding|INACTIVE',
                              'CTCF Binding Site|POISED;Enhancer|REPRESSED;TF binding|POISED',
                              'CTCF Binding Site|ACTIVE;Enhancer|REPRESSED;TF binding|ACTIVE')
  
  TFBS_ACTIVE_CLASS<-c('TF binding|ACTIVE;TF binding|INACTIVE',
                       'TF binding|ACTIVE;TF binding|REPRESSED',
                       'CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;TF binding|ACTIVE',
                       'CTCF Binding Site|ACTIVE;Open chromatin|ACTIVE;Open chromatin|INACTIVE;TF binding|ACTIVE','Open chromatin|REPRESSED;TF binding|ACTIVE','CTCF Binding Site|ACTIVE;TF binding|ACTIVE;TF binding|INACTIVE',
                       'CTCF Binding Site|INACTIVE;Open chromatin|ACTIVE;TF binding|ACTIVE',
                       'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Open chromatin|INACTIVE;TF binding|ACTIVE',
                       'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;TF binding|ACTIVE','TF binding|ACTIVE','Open chromatin|INACTIVE;TF binding|ACTIVE',
                       'CTCF Binding Site|POISED;TF binding|ACTIVE','CTCF Binding Site|ACTIVE;TF binding|ACTIVE',
                       'CTCF Binding Site|ACTIVE;TF binding|ACTIVE;TF binding|POISED','CTCF Binding Site|INACTIVE;TF binding|ACTIVE',
                       'CTCF Binding Site|ACTIVE;Open chromatin|ACTIVE;TF binding|ACTIVE','CTCF Binding Site|ACTIVE;Open chromatin|INACTIVE;TF binding|ACTIVE',
                       'CTCF Binding Site|INACTIVE;Open chromatin|INACTIVE;TF binding|ACTIVE',
                       'Open chromatin|ACTIVE;TF binding|ACTIVE','CTCF Binding Site|REPRESSED;Open chromatin|REPRESSED;TF binding|ACTIVE')
  
  
  
  
  TFBS_INACTIVE_CLASS<-c('TF binding|INACTIVE',
                         'CTCF Binding Site|ACTIVE;TF binding|INACTIVE',
                         'CTCF Binding Site|INACTIVE;TF binding|INACTIVE',
                         'CTCF Binding Site|ACTIVE;Open chromatin|INACTIVE;TF binding|INACTIVE',
                         'Open chromatin|INACTIVE;TF binding|INACTIVE','CTCF Binding Site|POISED;TF binding|INACTIVE',
                         'CTCF Binding Site|ACTIVE;Open chromatin|ACTIVE;TF binding|INACTIVE',
                         'CTCF Binding Site|INACTIVE;Open chromatin|INACTIVE;TF binding|INACTIVE')
  TFBS_POISED_CLASS<-c('TF binding|POISED',
                       'CTCF Binding Site|REPRESSED;TF binding|POISED;TF binding|REPRESSED',
                       'CTCF Binding Site|POISED;Open chromatin|REPRESSED;TF binding|POISED',
                       'CTCF Binding Site|ACTIVE;Open chromatin|REPRESSED;TF binding|POISED',
                       'CTCF Binding Site|POISED;TF binding|POISED',
                       'Open chromatin|POISED;TF binding|POISED','Open chromatin|REPRESSED;TF binding|POISED',
                       'TF binding|ACTIVE;TF binding|POISED','CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Open chromatin|POISED;TF binding|POISED',
                       'Open chromatin|POISED','CTCF Binding Site|ACTIVE;TF binding|POISED','CTCF Binding Site|REPRESSED;Open chromatin|REPRESSED;TF binding|POISED',
                       'CTCF Binding Site|REPRESSED;TF binding|POISED',
                       'CTCF Binding Site|POISED;Open chromatin|POISED;TF binding|POISED',
                       'CTCF Binding Site|POISED;TF binding|ACTIVE;TF binding|POISED',
                       'CTCF Binding Site|POISED;CTCF Binding Site|REPRESSED;TF binding|POISED',
                       'CTCF Binding Site|POISED;Open chromatin|INACTIVE;TF binding|POISED',
                       'Open chromatin|ACTIVE;TF binding|POISED')
  
  
  
  
  
  TFBS_REPRESSED_CLASS<-c('CTCF Binding Site|REPRESSED;Open chromatin|REPRESSED;TF binding|REPRESSED',
                            'TF binding|REPRESSED','CTCF Binding Site|REPRESSED;TF binding|REPRESSED',
                          'CTCF Binding Site|POISED;TF binding|REPRESSED',
                          'Open chromatin|REPRESSED;TF binding|REPRESSED',
                          'CTCF Binding Site|INACTIVE;TF binding|REPRESSED',
                          'CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Open chromatin|INACTIVE;TF binding|REPRESSED')
  
  Open_chromatin_ACTIVE_CLASS<-c('Open chromatin|ACTIVE','CTCF Binding Site|ACTIVE;Open chromatin|ACTIVE',
                                 'CTCF Binding Site|INACTIVE;Open chromatin|ACTIVE',
                                 'CTCF Binding Site|ACTIVE;Open chromatin|ACTIVE;Open chromatin|INACTIVE',
                                 'CTCF Binding Site|POISED;Open chromatin|ACTIVE')
  Open_chromatin_INACTIVE_CLASS<-c('Open chromatin|INACTIVE',
                                   'CTCF Binding Site|ACTIVE;Open chromatin|INACTIVE',
                                   'CTCF Binding Site|INACTIVE;Open chromatin|INACTIVE',
                                   'CTCF Binding Site|POISED;Open chromatin|INACTIVE',
                                   'CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Open chromatin|INACTIVE')
  Open_chromatin_POISED_CLASS<-c('Open_chromatin|POISED',
                                 'CTCF Binding Site|POISED;Open chromatin|POISED;Open chromatin|REPRESSED',
                                 'CTCF Binding Site|REPRESSED;Open chromatin|POISED',
                                 'CTCF Binding Site|POISED;Open chromatin|POISED',
                                 'CTCF Binding Site|ACTIVE;Open chromatin|POISED')
  Open_chromatin_REPRESSED_CLASS<-c('Open chromatin|REPRESSED',
                                    'CTCF Binding Site|REPRESSED;Open chromatin|REPRESSED',
                                    'CTCF Binding Site|ACTIVE;Open chromatin|REPRESSED',
                                    'CTCF Binding Site|POISED;Open chromatin|REPRESSED',
                                    'Open chromatin|INACTIVE;Open chromatin|REPRESSED',
                                    'CTCF Binding Site|INACTIVE;Open chromatin|REPRESSED',
                                    'CTCF Binding Site|POISED;CTCF Binding Site|REPRESSED;Open chromatin|REPRESSED')
  
  CTCF_ACTIVE_CLASS<-c('CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;CTCF Binding Site|REPRESSED',
                       'CTCF Binding Site|ACTIVE','CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED',
                       'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE',
                       'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED',
                       'CTCF Binding Site|ACTIVE;CTCF Binding Site|REPRESSED')
  CTCF_INACTIVE_CLASS<-c('CTCF Binding Site|INACTIVE','CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED',
                         'CTCF Binding Site|INACTIVE;CTCF Binding Site|REPRESSED')
  CTCF_POISED_CLASS<-c('CTCF Binding Site|POISED','CTCF Binding Site|POISED;CTCF Binding Site|REPRESSED',
                       'CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;CTCF Binding Site|REPRESSED')
  CTCF_REPRESSED_CLASS<-c('CTCF Binding Site|REPRESSED')
  
  #### Intersect Regulatory features ----
  
  DEBUG <- 1
  
  m <- findOverlaps(gr_K562_Regulatory_Build,gr_ALL_PoI_REST)
  
  if(DEBUG == 1)
  {
    cat("m\n")
    cat(str(m))
    cat("\n")
  }
  
  subjectHits_ALL_PoI_REST<-subjectHits(m)
  
  if(DEBUG == 1)
  {
    cat("subjectHits_ALL_PoI_REST\n")
    cat(str(subjectHits_ALL_PoI_REST))
    cat("\n")
  }
  
  queryHits_K562_Regulatory_Build<-queryHits(m)
  
  if(DEBUG == 1)
  {
    cat("queryHits_K562_Regulatory_Build\n")
    cat(str(queryHits_K562_Regulatory_Build))
    cat("\n")
  }
  
  K562_Regulatory_Build_df <- data.frame(chr=as.character(seqnames(gr_K562_Regulatory_Build)),
                                              start=as.integer(start(gr_K562_Regulatory_Build)),
                                              end=as.integer(end(gr_K562_Regulatory_Build)),
                                              feature_type=as.character(gr_K562_Regulatory_Build$name2),
                                              activity=as.character(gr_K562_Regulatory_Build$name3),
                                              regulatory_feature_stable_id=names(gr_K562_Regulatory_Build), stringsAsFactors = F)
  
  if(DEBUG == 1)
  {
    cat("K562_Regulatory_Build_df_0\n")
    cat(str(K562_Regulatory_Build_df))
    cat("\n")
  }
  
  K562_Regulatory_Build_df_hits<-K562_Regulatory_Build_df[queryHits_K562_Regulatory_Build,]
  
  if(DEBUG == 1)
  {
    cat("K562_Regulatory_Build_df_hits_0\n")
    cat(str(K562_Regulatory_Build_df_hits))
    cat("\n")
  }
  
  ALL_PoI_REST_Ensembl_REST_K562 <- data.frame(chr=as.character(seqnames(gr_ALL_PoI_REST)),
                                             start=as.integer(start(gr_ALL_PoI_REST)),
                                             end=as.integer(end(gr_ALL_PoI_REST)),
                                             Peak_ID=names(gr_ALL_PoI_REST),
                                             stringsAsFactors = F)
  
  
  
  if(DEBUG == 1)
  {
    cat("ALL_PoI_REST_Ensembl_REST_K562_0\n")
    cat(str(ALL_PoI_REST_Ensembl_REST_K562))
    cat("\n")
  }
  
  ALL_PoI_REST_Ensembl_REST_K562_hits<-ALL_PoI_REST_Ensembl_REST_K562[subjectHits_ALL_PoI_REST,]
  
  if(dim(ALL_PoI_REST_Ensembl_REST_K562_hits)[1] >0)
  {
    if(DEBUG == 1)
    {
      cat("ALL_PoI_REST_Ensembl_REST_K562_hits_0\n")
      cat(str(ALL_PoI_REST_Ensembl_REST_K562_hits))
      cat("\n")
    }
    
    ALL_PoI_REST_Ensembl_REST_K562_hits<-cbind(ALL_PoI_REST_Ensembl_REST_K562_hits,K562_Regulatory_Build_df_hits)
    
    if(DEBUG == 1)
    {
      cat("ALL_PoI_REST_Ensembl_REST_K562_hits_1\n")
      cat(str(ALL_PoI_REST_Ensembl_REST_K562_hits))
      cat("\n")
    }
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset<-unique(ALL_PoI_REST_Ensembl_REST_K562_hits[,c(1:4,c(which(colnames(ALL_PoI_REST_Ensembl_REST_K562_hits) == 'feature_type'),
                                                                                                which(colnames(ALL_PoI_REST_Ensembl_REST_K562_hits) == 'activity'),
                                                                                                which(colnames(ALL_PoI_REST_Ensembl_REST_K562_hits) == 'regulatory_feature_stable_id')))])
    
    if(DEBUG == 1)
    {
      cat("ALL_PoI_REST_Ensembl_REST_K562_hits_subset_0\n")
      cat(str(ALL_PoI_REST_Ensembl_REST_K562_hits_subset))
      cat("\n")
    }
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset$variable<-paste(ALL_PoI_REST_Ensembl_REST_K562_hits_subset$feature_type,ALL_PoI_REST_Ensembl_REST_K562_hits_subset$activity,sep='|')
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset$value<-paste(ALL_PoI_REST_Ensembl_REST_K562_hits_subset$regulatory_feature_stable_id,
                                                          ALL_PoI_REST_Ensembl_REST_K562_hits_subset$feature_type,
                                                          ALL_PoI_REST_Ensembl_REST_K562_hits_subset$activity,sep='|')
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset<-ALL_PoI_REST_Ensembl_REST_K562_hits_subset[,-c(which(colnames(ALL_PoI_REST_Ensembl_REST_K562_hits_subset) == 'regulatory_feature_stable_id'),
                                                                                           which(colnames(ALL_PoI_REST_Ensembl_REST_K562_hits_subset) == 'feature_type'),
                                                                                           which(colnames(ALL_PoI_REST_Ensembl_REST_K562_hits_subset) == 'activity'))]
    
    if(DEBUG == 1)
    {
      cat("ALL_PoI_REST_Ensembl_REST_K562_hits_subset_1\n")
      cat(str(ALL_PoI_REST_Ensembl_REST_K562_hits_subset))
      cat("\n")
    }
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset.dt<-data.table(ALL_PoI_REST_Ensembl_REST_K562_hits_subset, key=c("chr","start","end","Peak_ID"))
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed<-as.data.frame(ALL_PoI_REST_Ensembl_REST_K562_hits_subset.dt[,.(variable_string=paste(unique(sort(variable)), collapse=';'),
                                                                                                                     value_string=paste(value, collapse=';')),by=key(ALL_PoI_REST_Ensembl_REST_K562_hits_subset.dt)], stringsAsFactors=F)
    
    
    if(DEBUG == 1)
    {
      cat("ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed_0\n")
      cat(str(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed))
      cat("\n")
      cat(sprintf(paste(as.character(names(summary(as.factor(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string)))), collapse="BOOM")))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string)))))
      cat("\n")
    }
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS<-NA
    
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Promoter_ACTIVE_CLASS)]<-'Promoter_ACTIVE'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Promoter_INACTIVE_CLASS)]<-'Promoter_INACTIVE'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Promoter_POISED_CLASS)]<-'Promoter_POISED'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Promoter_REPRESSED_CLASS)]<-'Promoter_REPRESSED'
    
    
    
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_ACTIVE_CLASS)]<-'Enhancer_ACTIVE'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_INACTIVE_CLASS)]<-'Enhancer_INACTIVE'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_POISED_CLASS)]<-'Enhancer_POISED'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_REPRESSED_CLASS)]<-'Enhancer_REPRESSED'
    
    
    
    
    
   
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_ACTIVE_CLASS)]<-'OpenChromatin_ACTIVE'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_INACTIVE_CLASS)]<-'OpenChromatin_INACTIVE'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_POISED_CLASS)]<-'OpenChromatin_POISED'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_REPRESSED_CLASS)]<-'OpenChromatin_REPRESSED'
    
    
    
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_ACTIVE_CLASS)]<-'TFBS_ACTIVE'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_INACTIVE_CLASS)]<-'TFBS_INACTIVE'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_POISED_CLASS)]<-'TFBS_POISED'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_REPRESSED_CLASS)]<-'TFBS_REPRESSED'
    
    
    
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_ACTIVE_CLASS)]<-'CTCF_ACTIVE'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_INACTIVE_CLASS)]<-'CTCF_INACTIVE'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_POISED_CLASS)]<-'CTCF_POISED'
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_REPRESSED_CLASS)]<-'CTCF_REPRESSED'
    
    
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$feature<-gsub("_.+$","",ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS)
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$activity<-gsub("^[^_]+_","",ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS)
    
    check.NA<-ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed[is.na(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS),]
    
    
    # if(DEBUG == 1)
    # {
    cat("check.NA_ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed_0\n")
    cat(str(check.NA))
    cat("\n")
    cat(sprintf(paste(as.character(names(summary(as.factor(check.NA$variable_string)))), collapse="BOOM")))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(check.NA$variable_string)))))
    cat("\n")
    # }
    
    indx.dep<-c(which(colnames(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed) == 'variable_string'),
                which(colnames(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed) == 'CLASS'))
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed<-unique(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed[,-indx.dep])
    
    
    # if(DEBUG == 1)
    # {
    cat("ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed_1\n")
    cat(str(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$feature))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$feature)))))
    cat("\n")
    # }
    
    check.NA<-ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed[is.na(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$feature),]
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$Peak_CLASS<-'Non_gene_associated_peak'
    
    ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$Symbol_string<-NA
    
    cat("ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed_2\n")
    cat(str(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed))
    cat("\n")
   
    #### SAVE RDS ----
    
    
    
    setwd(out)
    
    saveRDS(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed,file='Rest_of_PoI_overlapping_regulatory.rds')
    write.table(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed,file='Rest_of_PoI_overlapping_regulatory.tsv',sep="\t",quote=F,row.names = F)
  }#dim(ALL_PoI_REST_Ensembl_REST_K562_hits)[1] >0
  
}

Consolidate_ALL_PoI = function(option_list)
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
  
  
  
  #### Read PoI_promoters file ----
  
  PoI_promoters<-readRDS(file=opt$PoI_promoters)
  
  
  cat("PoI_promoters_0\n")
  cat(str(PoI_promoters))
  cat("\n")
  cat(sprintf(as.character(names(summary(PoI_promoters$feature)))))
  cat("\n")
  cat(sprintf(as.character(summary(PoI_promoters$feature))))
  cat("\n")
  
  PoI_promoters<-unique(PoI_promoters[,-c(which(colnames(PoI_promoters) == 'variable_string'),
                                          which(colnames(PoI_promoters) == 'CLASS'))])
  
  cat("PoI_promoters_1\n")
  cat(str(PoI_promoters))
  cat("\n")
  cat(sprintf(as.character(names(summary(PoI_promoters$feature)))))
  cat("\n")
  cat(sprintf(as.character(summary(PoI_promoters$feature))))
  cat("\n")
  
  #### Read Non_promoters_PoI file ----
  
  Non_promoters_PoI<-readRDS(file=opt$Non_promoters_PoI)
  
  cat("Non_promoters_PoI_1\n")
  cat(str(Non_promoters_PoI))
  cat("\n")
  cat(sprintf(as.character(names(summary(Non_promoters_PoI$feature)))))
  cat("\n")
  cat(sprintf(as.character(summary(Non_promoters_PoI$feature))))
  cat("\n")
  
  #### Rest of PoI classified ------
  
  setwd(out)
  
  ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed<-readRDS(file='Rest_of_PoI_overlapping_regulatory.rds')
  
  cat("ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed_1\n")
  cat(str(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$feature)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$feature))))
  cat("\n")
  
  #### Read the ALL_PoI file ----
  
  ALL_PoI<-readRDS(file=opt$ALL_PoI)
  
  cat("ALL_PoI_0\n")
  cat(str(ALL_PoI))
  cat("\n")
  cat(str(unique(ALL_PoI$Peak_ID)))
  cat("\n")
  
  
  ALL_PoI_REST<-ALL_PoI[-which(ALL_PoI$Peak_ID%in%PoI_promoters$Peak_ID),]
  
  cat("ALL_PoI_REST_0\n")
  cat(str(ALL_PoI_REST))
  cat("\n")
  cat(str(unique(ALL_PoI_REST$Peak_ID)))
  cat("\n")
  
  ALL_PoI_REST<-ALL_PoI_REST[-which(ALL_PoI_REST$Peak_ID%in%Non_promoters_PoI$Peak_ID),]
  
  cat("ALL_PoI_REST_1\n")
  cat(str(ALL_PoI_REST))
  cat("\n")
  cat(str(unique(ALL_PoI_REST$Peak_ID)))
  cat("\n")
  
  
  ALL_PoI_REST<-ALL_PoI_REST[-which(ALL_PoI_REST$Peak_ID%in%ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$Peak_ID),]
  
  ALL_PoI_REST$chr<-gsub("^chr","",ALL_PoI_REST$chr)
  
  ALL_PoI_REST$value_string<-NA
  ALL_PoI_REST$feature<-NA
  ALL_PoI_REST$activity<-NA
  ALL_PoI_REST$Symbol_string<-NA
  
  ALL_PoI_REST$Peak_CLASS<-'Non_gene_associated_peak'
  
  
  cat("ALL_PoI_REST_1\n")
  cat(str(ALL_PoI_REST))
  cat("\n")
  cat(str(unique(ALL_PoI_REST$Peak_ID)))
  cat("\n")
  
  #### Rbind Master peak file -------
  
    Master_peak_file<-rbind(PoI_promoters,Non_promoters_PoI,
                            ALL_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed,ALL_PoI_REST)
  
  cat("Master_peak_file_0\n")
  cat(str(Master_peak_file))
  cat("\n")
  cat(str(unique(Master_peak_file$Peak_ID)))
  cat("\n")
    #### SAVE RDS ----
    
    
    
    setwd(out)
    
    saveRDS(Master_peak_file,file='Master_peak_file.rds')
    write.table(Master_peak_file,file='Master_peak_file.tsv',sep="\t",quote=F,row.names = F)

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
    make_option(c("--PoI_promoters"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Non_promoters_PoI"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_PoI"), type="character", default=NULL, 
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
  
  
  classify_the_REST_promoters_and_non_promoters(opt)
  Consolidate_ALL_PoI(opt)
  
  
}


###########################################################################

system.time( main() )