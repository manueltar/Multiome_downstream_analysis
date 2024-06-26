
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


classify_NON_PROMOTERS = function(option_list)
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
 
  # #### Read the ensembl_gtf ----
  # 
  # ensembl_gtf<-readGFF(file=opt$ensembl_gtf)
  # 
  # 
  # cat("ensembl_gtf_0\n")
  # cat(str(ensembl_gtf))
  # cat("\n")
  # cat(sprintf(as.character(names(summary((as.factor(ensembl_gtf$type)))))))
  # cat("\n")
  # cat(sprintf(as.character(summary((as.factor(ensembl_gtf$type))))))
  # cat("\n")
  # 
  # 
  # ensembl_gtf_gene<-unique(ensembl_gtf[which(ensembl_gtf$type == 'gene'),c(which(colnames(ensembl_gtf) == 'gene_id'),
  #                                                                          which(colnames(ensembl_gtf) == 'gene_name'),
  #                                                                          which(colnames(ensembl_gtf) == 'gene_biotype'),
  #                                                                          which(colnames(ensembl_gtf) == 'strand'))])
  # 
  # colnames(ensembl_gtf_gene)[which(colnames(ensembl_gtf_gene) == 'gene_id')]<-'ensembl_gene_id'
  # colnames(ensembl_gtf_gene)[which(colnames(ensembl_gtf_gene) == 'gene_name')]<-'Symbol'
  # 
  # 
  # cat("ensembl_gtf_gene_0\n")
  # cat(str(ensembl_gtf_gene))
  # cat("\n")
  # cat(str(unique(ensembl_gtf_gene$gene_id)))
  # cat("\n")
  # 
  # rm(ensembl_gtf)
  
  #### Read the K562_Regulatory_Build ----
  
  K562_Regulatory_Build<-readGFF(file=opt$K562_Regulatory_Build)
 

  cat("K562_Regulatory_Build_0\n")
  cat(str(K562_Regulatory_Build))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_Regulatory_Build$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_Regulatory_Build$type))))
  cat("\n")
  
  K562_Regulatory_Build_REST<-K562_Regulatory_Build[which(K562_Regulatory_Build$type != 'promoter'),]
  
  
  cat("K562_Regulatory_Build_REST_0\n")
  cat(str(K562_Regulatory_Build_REST))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_Regulatory_Build_REST$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_Regulatory_Build_REST$type))))
  cat("\n")
  
  gr_K562_Regulatory_Build_REST <- GRanges(
    seqnames = as.character(K562_Regulatory_Build_REST$seqid),
    name2=as.character(K562_Regulatory_Build_REST$feature_type),
    name3=as.character(K562_Regulatory_Build_REST$activity),
    ranges=IRanges(
      start=as.numeric(K562_Regulatory_Build_REST$start),
      end=as.numeric(K562_Regulatory_Build_REST$end),
      name=K562_Regulatory_Build_REST$regulatory_feature_stable_id))
  
  #### Read promoters file ----
  
  PoI_promoters<-readRDS(file=opt$PoI_promoters)
  
  
  cat("PoI_promoters_0\n")
  cat(str(PoI_promoters))
  cat("\n")
  cat(sprintf(as.character(names(summary(PoI_promoters$Peak_ID)))))
  cat("\n")
  cat(sprintf(as.character(summary(PoI_promoters$Peak_ID))))
  cat("\n")
  
  #### Read the Peak to genes results file ----
  
  TSS_PoI<-readRDS(file=opt$TSS_PoI)
  
  
  cat("TSS_PoI_0\n")
  cat(str(TSS_PoI))
  cat("\n")
  
  
  TSS_PoI_REST<-TSS_PoI[-which(TSS_PoI$Peak_ID%in%PoI_promoters$Peak_ID),]
  
  cat("TSS_PoI_REST_0\n")
  cat(str(TSS_PoI_REST))
  cat("\n")
 
  
  gr_TSS_PoI_REST <- GRanges(
    seqnames = as.character(TSS_PoI_REST$chr),
    name2=as.character(TSS_PoI_REST$TSS_Symbol),
    name3=as.character(TSS_PoI_REST$ensembl_gene_id),
    name4=as.character(TSS_PoI_REST$gene_biotype),
    name5=as.character(TSS_PoI_REST$strand),
    ranges=IRanges(
      start=as.numeric(TSS_PoI_REST$start),
      end=as.numeric(TSS_PoI_REST$end),
      name=TSS_PoI_REST$Peak_ID))
  
  # cat("gr_TSS_PoI_REST_0\n")
  # cat(str(gr_TSS_PoI_REST))
  # cat("\n")
  
  #### Read the Peak to genes results file ----
  
  Linked_PoI<-readRDS(file=opt$Linked_PoI)
  
  
  cat("Linked_PoI_0\n")
  cat(str(Linked_PoI))
  cat("\n")
 
  Linked_PoI_REST<-Linked_PoI[-which(Linked_PoI$Peak_ID%in%PoI_promoters$Peak_ID),]
  
  cat("Linked_PoI_REST_0\n")
  cat(str(Linked_PoI_REST))
  cat("\n")
  
  Linked_PoI_REST<-Linked_PoI_REST[-which(Linked_PoI_REST$Peak_ID%in%TSS_PoI_REST$Peak_ID),]
  
  cat("Linked_PoI_REST_1\n")
  cat(str(Linked_PoI_REST))
  cat("\n")
  
  
  gr_Linked_PoI_REST <- GRanges(
    seqnames = as.character(gsub("^chr","",gsub("-.+$","",Linked_PoI_REST$Peak_ID))),
    name2=as.character(Linked_PoI_REST$Symbol),
    ranges=IRanges(
      start=as.integer(gsub("-.+$","",gsub("^[^-]+-","",Linked_PoI_REST$Peak_ID))), 
      end=as.integer(gsub("^[^-]+-[^-]+-","",Linked_PoI_REST$Peak_ID)),
      name=Linked_PoI_REST$Peak_ID))
  
  cat("gr_Linked_PoI_REST_0\n")
  cat(str(gr_Linked_PoI_REST))
  cat("\n")
  
  #### CLasses of tags ------------------------
  
  Enhancer_ACTIVE_CLASS<-c('Enhancer|ACTIVE','Enhancer|ACTIVE;Enhancer|INACTIVE','Enhancer|ACTIVE;TF binding|ACTIVE',
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
                           'CTCF Binding Site|POISED;Open chromatin|REPRESSED;TF binding|ACTIVE')
  
  Enhancer_INACTIVE_CLASS<-c('CTCF Binding Site|POISED;Enhancer|INACTIVE',
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
                             'CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Enhancer|INACTIVE')
  Enhancer_POISED_CLASS<-c('Enhancer|POISED','Enhancer|POISED;TF binding|POISED','CTCF Binding Site|POISED;Enhancer|POISED','CTCF Binding Site|INACTIVE;Enhancer|POISED',
                           'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;Enhancer|POISED','CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Enhancer|POISED',
                           'CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED;Enhancer|POISED',
                           'CTCF Binding Site|POISED;Enhancer|POISED;Enhancer|REPRESSED',
                           'CTCF Binding Site|POISED;Enhancer|POISED;Open chromatin|POISED',
                           'CTCF Binding Site|POISED;Enhancer|POISED;TF binding|INACTIVE',
                           'CTCF Binding Site|POISED;Enhancer|POISED;TF binding|POISED',
                           'Enhancer|POISED;Open chromatin|INACTIVE',
                           'CTCF Binding Site|ACTIVE;Enhancer|POISED',
                           'CTCF Binding Site|REPRESSED;Enhancer|POISED')
  
  Enhancer_REPRESSED_CLASS<-c('Enhancer|REPRESSED','Enhancer|REPRESSED;Open chromatin|REPRESSED','CTCF Binding Site|INACTIVE;Enhancer|REPRESSED','CTCF Binding Site|POISED;Enhancer|REPRESSED',
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
                              'CTCF Binding Site|REPRESSED;Enhancer|REPRESSED;Open chromatin|INACTIVE')
  
  TFBS_ACTIVE_CLASS<-c('TF binding|ACTIVE','Open chromatin|INACTIVE;TF binding|ACTIVE',
                       'CTCF Binding Site|POISED;TF binding|ACTIVE','CTCF Binding Site|ACTIVE;TF binding|ACTIVE',
                       'CTCF Binding Site|ACTIVE;TF binding|ACTIVE;TF binding|POISED','CTCF Binding Site|INACTIVE;TF binding|ACTIVE',
                       'CTCF Binding Site|ACTIVE;Open chromatin|ACTIVE;TF binding|ACTIVE','CTCF Binding Site|ACTIVE;Open chromatin|INACTIVE;TF binding|ACTIVE',
                       'CTCF Binding Site|INACTIVE;Open chromatin|INACTIVE;TF binding|ACTIVE',
                       'Open chromatin|ACTIVE;TF binding|ACTIVE')
  TFBS_INACTIVE_CLASS<-c('TF binding|INACTIVE',
                         'CTCF Binding Site|ACTIVE;TF binding|INACTIVE',
                         'CTCF Binding Site|INACTIVE;TF binding|INACTIVE',
                         'CTCF Binding Site|ACTIVE;Open chromatin|INACTIVE;TF binding|INACTIVE')
  TFBS_POISED_CLASS<-c('TF binding|POISED',
                       'CTCF Binding Site|POISED;TF binding|POISED',
                       'Open chromatin|POISED;TF binding|POISED','Open chromatin|REPRESSED;TF binding|POISED',
                       'TF binding|ACTIVE;TF binding|POISED','CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED;Open chromatin|POISED;TF binding|POISED',
                       'Open chromatin|POISED','CTCF Binding Site|ACTIVE;TF binding|POISED','CTCF Binding Site|REPRESSED;Open chromatin|REPRESSED;TF binding|POISED',
                       'CTCF Binding Site|REPRESSED;TF binding|POISED')
  TFBS_REPRESSED_CLASS<-c('TF binding|REPRESSED','CTCF Binding Site|REPRESSED;TF binding|REPRESSED')
  
  Open_chromatin_ACTIVE_CLASS<-c('Open chromatin|ACTIVE','CTCF Binding Site|ACTIVE;Open chromatin|ACTIVE',
                                 'CTCF Binding Site|INACTIVE;Open chromatin|ACTIVE')
  Open_chromatin_INACTIVE_CLASS<-c('Open chromatin|INACTIVE',
                                   'CTCF Binding Site|ACTIVE;Open chromatin|INACTIVE',
                                   'CTCF Binding Site|INACTIVE;Open chromatin|INACTIVE',
                                   'CTCF Binding Site|POISED;Open chromatin|INACTIVE')
  Open_chromatin_POISED_CLASS<-c('Open_chromatin|POISED',
                                 'CTCF Binding Site|POISED;Open chromatin|POISED;Open chromatin|REPRESSED',
                                 'CTCF Binding Site|REPRESSED;Open chromatin|POISED')
  Open_chromatin_REPRESSED_CLASS<-c('Open chromatin|REPRESSED',
                                    'CTCF Binding Site|REPRESSED;Open chromatin|REPRESSED',
                                    'CTCF Binding Site|ACTIVE;Open chromatin|REPRESSED')
  
  CTCF_ACTIVE_CLASS<-c('CTCF Binding Site|ACTIVE','CTCF Binding Site|ACTIVE;CTCF Binding Site|POISED',
                       'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE',
                       'CTCF Binding Site|ACTIVE;CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED',
                       'CTCF Binding Site|ACTIVE;CTCF Binding Site|REPRESSED')
  CTCF_INACTIVE_CLASS<-c('CTCF Binding Site|INACTIVE','CTCF Binding Site|INACTIVE;CTCF Binding Site|POISED')
  CTCF_POISED_CLASS<-c('CTCF Binding Site|POISED','CTCF Binding Site|POISED;CTCF Binding Site|REPRESSED')
  CTCF_REPRESSED_CLASS<-c('CTCF Binding Site|REPRESSED')
  
  #### Intersect REST TSS Peak to genes with REST Regulatory features ----
  
  DEBUG <- 0
  
  m <- findOverlaps(gr_K562_Regulatory_Build_REST,gr_Linked_PoI_REST)
  
  if(DEBUG == 1)
  {
    cat("m\n")
    cat(str(m))
    cat("\n")
  }
  
  subjectHits_Linked_PoI_REST<-subjectHits(m)
  
  if(DEBUG == 1)
  {
    cat("subjectHits_Linked_PoI_REST\n")
    cat(str(subjectHits_Linked_PoI_REST))
    cat("\n")
  }
  
  queryHits_K562_Regulatory_Build_REST<-queryHits(m)
  
  if(DEBUG == 1)
  {
    cat("queryHits_K562_Regulatory_Build_REST\n")
    cat(str(queryHits_K562_Regulatory_Build_REST))
    cat("\n")
  }
  
  K562_Regulatory_Build_REST_df <- data.frame(chr=as.character(seqnames(gr_K562_Regulatory_Build_REST)),
                                              start=as.integer(start(gr_K562_Regulatory_Build_REST)),
                                              end=as.integer(end(gr_K562_Regulatory_Build_REST)),
                                              feature_type=as.character(gr_K562_Regulatory_Build_REST$name2),
                                              activity=as.character(gr_K562_Regulatory_Build_REST$name3),
                                              regulatory_feature_stable_id=names(gr_K562_Regulatory_Build_REST), stringsAsFactors = F)
  
  if(DEBUG == 1)
  {
    cat("K562_Regulatory_Build_REST_df_0\n")
    cat(str(K562_Regulatory_Build_REST_df))
    cat("\n")
  }
  
  K562_Regulatory_Build_REST_df_hits<-K562_Regulatory_Build_REST_df[queryHits_K562_Regulatory_Build_REST,]
  
  if(DEBUG == 1)
  {
    cat("K562_Regulatory_Build_REST_df_hits_0\n")
    cat(str(K562_Regulatory_Build_REST_df_hits))
    cat("\n")
  }
  
  Linked_PoI_REST_Ensembl_REST_K562 <- data.frame(chr=as.character(seqnames(gr_Linked_PoI_REST)),
                                             start=as.integer(start(gr_Linked_PoI_REST)),
                                             end=as.integer(end(gr_Linked_PoI_REST)),
                                             Symbol=as.character(gr_Linked_PoI_REST$name2),
                                             Peak_ID=names(gr_Linked_PoI_REST),
                                             stringsAsFactors = F)
  
  
  
  if(DEBUG == 1)
  {
    cat("Linked_PoI_REST_Ensembl_REST_K562_0\n")
    cat(str(Linked_PoI_REST_Ensembl_REST_K562))
    cat("\n")
  }
  
  Linked_PoI_REST_Ensembl_REST_K562_hits<-Linked_PoI_REST_Ensembl_REST_K562[subjectHits_Linked_PoI_REST,]
  
  if(dim(Linked_PoI_REST_Ensembl_REST_K562_hits)[1] >0)
  {
    if(DEBUG == 1)
    {
      cat("Linked_PoI_REST_Ensembl_REST_K562_hits_0\n")
      cat(str(Linked_PoI_REST_Ensembl_REST_K562_hits))
      cat("\n")
    }
    
    Linked_PoI_REST_Ensembl_REST_K562_hits<-cbind(Linked_PoI_REST_Ensembl_REST_K562_hits,K562_Regulatory_Build_REST_df_hits)
    
    if(DEBUG == 1)
    {
      cat("Linked_PoI_REST_Ensembl_REST_K562_hits_1\n")
      cat(str(Linked_PoI_REST_Ensembl_REST_K562_hits))
      cat("\n")
    }
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset<-unique(Linked_PoI_REST_Ensembl_REST_K562_hits[,c(1:5,c(which(colnames(Linked_PoI_REST_Ensembl_REST_K562_hits) == 'feature_type'),
                                                                                                which(colnames(Linked_PoI_REST_Ensembl_REST_K562_hits) == 'activity'),
                                                                                                which(colnames(Linked_PoI_REST_Ensembl_REST_K562_hits) == 'regulatory_feature_stable_id')))])
    
    if(DEBUG == 1)
    {
      cat("Linked_PoI_REST_Ensembl_REST_K562_hits_subset_0\n")
      cat(str(Linked_PoI_REST_Ensembl_REST_K562_hits_subset))
      cat("\n")
    }
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset$variable<-paste(Linked_PoI_REST_Ensembl_REST_K562_hits_subset$feature_type,Linked_PoI_REST_Ensembl_REST_K562_hits_subset$activity,sep='|')
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset$value<-paste(Linked_PoI_REST_Ensembl_REST_K562_hits_subset$regulatory_feature_stable_id,
                                                          Linked_PoI_REST_Ensembl_REST_K562_hits_subset$feature_type,
                                                          Linked_PoI_REST_Ensembl_REST_K562_hits_subset$activity,sep='|')
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset<-Linked_PoI_REST_Ensembl_REST_K562_hits_subset[,-c(which(colnames(Linked_PoI_REST_Ensembl_REST_K562_hits_subset) == 'regulatory_feature_stable_id'),
                                                                                           which(colnames(Linked_PoI_REST_Ensembl_REST_K562_hits_subset) == 'feature_type'),
                                                                                           which(colnames(Linked_PoI_REST_Ensembl_REST_K562_hits_subset) == 'activity'))]
    
    if(DEBUG == 1)
    {
      cat("Linked_PoI_REST_Ensembl_REST_K562_hits_subset_1\n")
      cat(str(Linked_PoI_REST_Ensembl_REST_K562_hits_subset))
      cat("\n")
    }
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset.dt<-data.table(Linked_PoI_REST_Ensembl_REST_K562_hits_subset, key=c("chr","start","end","Peak_ID"))
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed<-as.data.frame(Linked_PoI_REST_Ensembl_REST_K562_hits_subset.dt[,.(variable_string=paste(unique(sort(variable)), collapse=';'),
                                                                                                                     value_string=paste(value, collapse=';'),
                                                                                                                     Symbol_string=paste(unique(sort(Symbol)), collapse=';')),by=key(Linked_PoI_REST_Ensembl_REST_K562_hits_subset.dt)], stringsAsFactors=F)
    
    
    if(DEBUG == 1)
    {
      cat("Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed_0\n")
      cat(str(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed))
      cat("\n")
      cat(sprintf(paste(as.character(names(summary(as.factor(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string)))), collapse="BOOM")))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string)))))
      cat("\n")
    }
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS<-NA
    
    
    
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_ACTIVE_CLASS)]<-'Enhancer_ACTIVE'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_INACTIVE_CLASS)]<-'Enhancer_INACTIVE'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_POISED_CLASS)]<-'Enhancer_POISED'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_REPRESSED_CLASS)]<-'Enhancer_REPRESSED'
    
    
    
    
    
   
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_ACTIVE_CLASS)]<-'OpenChromatin_ACTIVE'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_INACTIVE_CLASS)]<-'OpenChromatin_INACTIVE'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_POISED_CLASS)]<-'OpenChromatin_POISED'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_REPRESSED_CLASS)]<-'OpenChromatin_REPRESSED'
    
    
    
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_ACTIVE_CLASS)]<-'TFBS_ACTIVE'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_INACTIVE_CLASS)]<-'TFBS_INACTIVE'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_POISED_CLASS)]<-'TFBS_POISED'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_REPRESSED_CLASS)]<-'TFBS_REPRESSED'
    
    
    
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_ACTIVE_CLASS)]<-'CTCF_ACTIVE'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_INACTIVE_CLASS)]<-'CTCF_INACTIVE'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_POISED_CLASS)]<-'CTCF_POISED'
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_REPRESSED_CLASS)]<-'CTCF_REPRESSED'
    
    
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$feature<-gsub("_.+$","",Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS)
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$activity<-gsub("^[^_]+_","",Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS)
    
    check.NA<-Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed[is.na(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS),]
    
    
    # if(DEBUG == 1)
    # {
    cat("check.NA_Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed_0\n")
    cat(str(check.NA))
    cat("\n")
    cat(sprintf(paste(as.character(names(summary(as.factor(check.NA$variable_string)))), collapse="BOOM")))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(check.NA$variable_string)))))
    cat("\n")
    # }
    
    indx.dep<-c(which(colnames(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed) == 'variable_string'),
                which(colnames(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed) == 'CLASS'))
    
    Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed<-unique(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed[,-indx.dep])
    
    
    # if(DEBUG == 1)
    # {
    cat("Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed_1\n")
    cat(str(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$CLASS)))))
    cat("\n")
    # }
    
    check.NA<-Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed[is.na(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$feature),]
    
   
    
  }#dim(Linked_PoI_REST_Ensembl_REST_K562_hits)[1] >0
  
  
 
  
 
  #### Intersect REST TSS Peak to genes with REST Regulatory features ----
  
  DEBUG <- 0
  
  m <- findOverlaps(gr_K562_Regulatory_Build_REST,gr_TSS_PoI_REST)
  
  if(DEBUG == 1)
  {
    cat("m\n")
    cat(str(m))
    cat("\n")
  }
  
  subjectHits_TSS_PoI_REST<-subjectHits(m)
  
  if(DEBUG == 1)
  {
    cat("subjectHits_TSS_PoI_REST\n")
    cat(str(subjectHits_TSS_PoI_REST))
    cat("\n")
  }
  
  queryHits_K562_Regulatory_Build_REST<-queryHits(m)
  
  if(DEBUG == 1)
  {
    cat("queryHits_K562_Regulatory_Build_REST\n")
    cat(str(queryHits_K562_Regulatory_Build_REST))
    cat("\n")
  }
  
  K562_Regulatory_Build_REST_df <- data.frame(chr=as.character(seqnames(gr_K562_Regulatory_Build_REST)),
                        start=as.integer(start(gr_K562_Regulatory_Build_REST)),
                        end=as.integer(end(gr_K562_Regulatory_Build_REST)),
                        feature_type=as.character(gr_K562_Regulatory_Build_REST$name2),
                        activity=as.character(gr_K562_Regulatory_Build_REST$name3),
                        regulatory_feature_stable_id=names(gr_K562_Regulatory_Build_REST), stringsAsFactors = F)
  
  if(DEBUG == 1)
  {
    cat("K562_Regulatory_Build_REST_df_0\n")
    cat(str(K562_Regulatory_Build_REST_df))
    cat("\n")
  }
  
  K562_Regulatory_Build_REST_df_hits<-K562_Regulatory_Build_REST_df[queryHits_K562_Regulatory_Build_REST,]
  
  if(DEBUG == 1)
  {
    cat("K562_Regulatory_Build_REST_df_hits_0\n")
    cat(str(K562_Regulatory_Build_REST_df_hits))
    cat("\n")
  }
  
  TSS_PoI_df_Ensembl_REST_K562 <- data.frame(chr=as.character(seqnames(gr_TSS_PoI_REST)),
                                                  start=as.integer(start(gr_TSS_PoI_REST)),
                                                  end=as.integer(end(gr_TSS_PoI_REST)),
                                                  Symbol=as.character(gr_TSS_PoI_REST$name2),
                                                  ensembl_gene_id=as.character(gr_TSS_PoI_REST$name3),
                                                  gene_biotype=as.character(gr_TSS_PoI_REST$name4),
                                                  strand=as.character(gr_TSS_PoI_REST$name5),
                                                  Peak_ID=names(gr_TSS_PoI_REST),
                                                  stringsAsFactors = F)
  
  
  
  if(DEBUG == 1)
  {
    cat("TSS_PoI_df_Ensembl_REST_K562_0\n")
    cat(str(TSS_PoI_df_Ensembl_REST_K562))
    cat("\n")
  }
  
  TSS_PoI_df_Ensembl_REST_K562_hits<-TSS_PoI_df_Ensembl_REST_K562[subjectHits_TSS_PoI_REST,]
  
  if(dim(TSS_PoI_df_Ensembl_REST_K562_hits)[1] >0)
  {
    if(DEBUG == 1)
    {
      cat("TSS_PoI_df_Ensembl_REST_K562_hits_0\n")
      cat(str(TSS_PoI_df_Ensembl_REST_K562_hits))
      cat("\n")
    }
    
    TSS_PoI_df_Ensembl_REST_K562_hits<-cbind(TSS_PoI_df_Ensembl_REST_K562_hits,K562_Regulatory_Build_REST_df_hits)
    
    if(DEBUG == 1)
    {
      cat("TSS_PoI_df_Ensembl_REST_K562_hits_1\n")
      cat(str(TSS_PoI_df_Ensembl_REST_K562_hits))
      cat("\n")
    }
    
    TSS_PoI_df_Ensembl_REST_K562_hits_subset<-unique(TSS_PoI_df_Ensembl_REST_K562_hits[,c(1:8,c(which(colnames(TSS_PoI_df_Ensembl_REST_K562_hits) == 'feature_type'),
                                                                                which(colnames(TSS_PoI_df_Ensembl_REST_K562_hits) == 'activity'),
                                                                                which(colnames(TSS_PoI_df_Ensembl_REST_K562_hits) == 'regulatory_feature_stable_id')))])

    if(DEBUG == 1)
    {
      cat("TSS_PoI_df_Ensembl_REST_K562_hits_subset_0\n")
      cat(str(TSS_PoI_df_Ensembl_REST_K562_hits_subset))
      cat("\n")
    }

    TSS_PoI_df_Ensembl_REST_K562_hits_subset$variable<-paste(TSS_PoI_df_Ensembl_REST_K562_hits_subset$feature_type,TSS_PoI_df_Ensembl_REST_K562_hits_subset$activity,sep='|')
    TSS_PoI_df_Ensembl_REST_K562_hits_subset$value<-paste(TSS_PoI_df_Ensembl_REST_K562_hits_subset$regulatory_feature_stable_id,
                                                                    TSS_PoI_df_Ensembl_REST_K562_hits_subset$feature_type,
                                                                    TSS_PoI_df_Ensembl_REST_K562_hits_subset$activity,sep='|')

    TSS_PoI_df_Ensembl_REST_K562_hits_subset<-TSS_PoI_df_Ensembl_REST_K562_hits_subset[,-c(which(colnames(TSS_PoI_df_Ensembl_REST_K562_hits_subset) == 'regulatory_feature_stable_id'),
                                                                           which(colnames(TSS_PoI_df_Ensembl_REST_K562_hits_subset) == 'feature_type'),
                                                                           which(colnames(TSS_PoI_df_Ensembl_REST_K562_hits_subset) == 'activity'))]

    if(DEBUG == 1)
    {
      cat("TSS_PoI_df_Ensembl_REST_K562_hits_subset_1\n")
      cat(str(TSS_PoI_df_Ensembl_REST_K562_hits_subset))
      cat("\n")
    }
    
    TSS_PoI_df_Ensembl_REST_K562_hits_subset.dt<-data.table(TSS_PoI_df_Ensembl_REST_K562_hits_subset, key=c("chr","start","end","Peak_ID"))
    
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed<-as.data.frame(TSS_PoI_df_Ensembl_REST_K562_hits_subset.dt[,.(variable_string=paste(unique(sort(variable)), collapse=';'),
                                                                                                                                                   value_string=paste(value, collapse=';'),
                                                                                                                     Symbol_string=paste(unique(sort(Symbol)), collapse=';')),by=key(TSS_PoI_df_Ensembl_REST_K562_hits_subset.dt)], stringsAsFactors=F)
    
    
    if(DEBUG == 1)
    {
      cat("TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed_0\n")
      cat(str(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed))
      cat("\n")
      cat(sprintf(paste(as.character(names(summary(as.factor(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string)))), collapse="BOOM")))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string)))))
      cat("\n")
    }
    
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS<-NA
    
   
    
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_ACTIVE_CLASS)]<-'Enhancer_ACTIVE'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_INACTIVE_CLASS)]<-'Enhancer_INACTIVE'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_POISED_CLASS)]<-'Enhancer_POISED'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Enhancer_REPRESSED_CLASS)]<-'Enhancer_REPRESSED'
    
    
    
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_ACTIVE_CLASS)]<-'OpenChromatin_ACTIVE'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_INACTIVE_CLASS)]<-'OpenChromatin_INACTIVE'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_POISED_CLASS)]<-'OpenChromatin_POISED'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%Open_chromatin_REPRESSED_CLASS)]<-'OpenChromatin_REPRESSED'
    
    
  
    
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_ACTIVE_CLASS)]<-'TFBS_ACTIVE'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_INACTIVE_CLASS)]<-'TFBS_INACTIVE'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_POISED_CLASS)]<-'TFBS_POISED'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%TFBS_REPRESSED_CLASS)]<-'TFBS_REPRESSED'
    
    
    
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_ACTIVE_CLASS)]<-'CTCF_ACTIVE'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_INACTIVE_CLASS)]<-'CTCF_INACTIVE'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_POISED_CLASS)]<-'CTCF_POISED'
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS[which(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$variable_string%in%CTCF_REPRESSED_CLASS)]<-'CTCF_REPRESSED'
    
    
    
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$feature<-gsub("_.+$","",TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS)
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$activity<-gsub("^[^_]+_","",TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS)
    
    check.NA<-TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed[is.na(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS),]
    
    
    # if(DEBUG == 1)
    # {
    cat("check.NA_TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed_0\n")
    cat(str(check.NA))
    cat("\n")
    cat(sprintf(paste(as.character(names(summary(as.factor(check.NA$variable_string)))), collapse="BOOM")))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(check.NA$variable_string)))))
    cat("\n")
    # }
    # 
    
    indx.dep<-c(which(colnames(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed) == 'variable_string'),
                which(colnames(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed) == 'CLASS'))
    
    TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed<-unique(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed[,-indx.dep])
    
    # if(DEBUG == 1)
    # {
      cat("TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed_1\n")
      cat(str(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$CLASS)))))
      cat("\n")
    # }
    
    
    
   
    
  }#dim(TSS_PoI_df_Ensembl_REST_K562_hits)[1] >0
  
  ################ Merge ---------------------------------
  
  Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed$Peak_CLASS<-'Linked_Peak'
  
  TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed$Peak_CLASS<-'TSS_Peak'
  
  
  cat("TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed_REMEMBER\n")
  cat(str(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed))
  cat("\n")
  
  cat("Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed_REMEMBER\n")
  cat(str(Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed))
  cat("\n")
  
  DEF_collapsed<-rbind(TSS_PoI_df_Ensembl_REST_K562_hits_subset_collapsed,
                       Linked_PoI_REST_Ensembl_REST_K562_hits_subset_collapsed)

  # if(DEBUG == 1)
  # {
    cat("DEF_collapsed_0\n")
    cat(str(DEF_collapsed))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(DEF_collapsed$feature))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(DEF_collapsed$feature)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(DEF_collapsed$activity))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(DEF_collapsed$activity)))))
    cat("\n")
  # }
  
  DEF_collapsed$Peak_CLASS<-factor(DEF_collapsed$Peak_CLASS,
                                        levels=c('Linked_Peak','TSS_Peak'),
                                        ordered=T)


  # DEF_collapsed$feature<-factor(DEF_collapsed$feature,
  #                                    levels=c('Promoter','Enhancer','TFBS','OpenChromatin','CTCF'),
  #                                    ordered=T)
  # 
  # DEF_collapsed$activity<-factor(DEF_collapsed$activity,
  #                                      levels=c("ACTIVE",'INACTIVE','POISED',"REPRESSED"),
  #                                      ordered=T)


  if(DEBUG == 1)
  {
    cat("DEF_collapsed_0\n")
    cat(str(DEF_collapsed))
    cat("\n")
    cat(sprintf(as.character(names(summary(DEF_collapsed$feature)))))
    cat("\n")
    cat(sprintf(as.character(summary(DEF_collapsed$feature))))
    cat("\n")
    cat(sprintf(as.character(names(summary(DEF_collapsed$activity)))))
    cat("\n")
    cat(sprintf(as.character(summary(DEF_collapsed$activity))))
    cat("\n")
  }

  check<-DEF_collapsed[grep(paste(tracking_genes,
                                                 collapse='|'),DEF_collapsed$Symbol_string),]
  

  # if(DEBUG == 1)
  # {
  cat("check_3\n")
  cat(str(check))
  cat("\n")
  cat(str(unique(check$Symbol)))
  cat("\n")
  # }
  
  check.NA<-DEF_collapsed[is.na(DEF_collapsed$feature),]
  
  cat("check.NA_0\n")
  cat(str(check.NA))
  cat("\n")
  cat(str(unique(check.NA$Symbol)))
  cat("\n")

  #### SAVE RDS ----



  setwd(out)

  saveRDS(DEF_collapsed,file='PoI_non_promoters.rds')
  write.table(DEF_collapsed,file='PoI_non_promoters.tsv',sep="\t",quote=F,row.names = F)
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
    make_option(c("--TSS_PoI"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Linked_PoI"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ensembl_gtf"), type="character", default=NULL, 
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
  
  
  classify_NON_PROMOTERS(opt)
  
  
}


###########################################################################

system.time( main() )