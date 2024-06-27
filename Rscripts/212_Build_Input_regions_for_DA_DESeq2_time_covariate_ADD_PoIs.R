
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
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)


Build_input_regions = function(option_list)
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
  
  #### READ and transform Peaks_to_display ----
  
  Peaks_to_display = unlist(strsplit(opt$Peaks_to_display, split=','))
  
  cat("Peaks_to_display_\n")
  cat(sprintf(as.character(Peaks_to_display)))
  cat("\n")
  
  Peaks_to_display_df <- data.frame(matrix(vector(), length(Peaks_to_display), 2,
                                        dimnames=list(c(),
                                                      c('Symbol','Peaks_to_display'))),stringsAsFactors=F)
  
  
  
  Peaks_to_display_df$Symbol<-gsub("\\>.+","",Peaks_to_display)
  Peaks_to_display_df$Peaks_to_display<-Peaks_to_display
  
  
  cat("Peaks_to_display_df_0\n")
  cat(str(Peaks_to_display_df))
  cat("\n")
  
  Peaks_to_display_df$Peaks_to_display<-gsub(paste(paste(Peaks_to_display_df$Symbol,'>', sep=''), collapse='|'),"",Peaks_to_display_df$Peaks_to_display)
  
  
  cat("Peaks_to_display_df_1\n")
  cat(str(Peaks_to_display_df))
  cat("\n")
  
  Peaks_to_display_df$Peaks_to_display<-gsub("\\>","",Peaks_to_display_df$Peaks_to_display)
  
  cat("Peaks_to_display_df_2\n")
  cat(str(Peaks_to_display_df))
  cat("\n")
  
  
  
  Peaks_to_display_df_long<-unique(as.data.frame(cSplit(Peaks_to_display_df,sep = '|', direction = "long",
                                                                             splitCols = "Peaks_to_display"),stringsAsFactors=F))
  

  cat("Peaks_to_display_df_long\n")
  cat(str(Peaks_to_display_df_long))
  cat("\n")
  
  #### READ Input_regions ----
  
  Input_regions<-readRDS(file=opt$Input_regions)
  
  
  cat("Input_regions\n")
  cat(str(Input_regions))
  cat("\n")
  
  check_CUX1<-Input_regions[which(Input_regions$Symbol == 'CUX1'),]
  
  cat("check_CUX1_0\n")
  cat(str(check_CUX1))
  cat("\n")
  
  Display_genes_array<-unique(unlist(strsplit(Input_regions$Display_genes, split=";")))
  
  cat("Display_genes_array\n")
  cat(str(Display_genes_array))
  cat("\n")
  
  
  check_CUX1<-Display_genes_array[which(Display_genes_array == 'CUX1')]
  
  cat("check_CUX1_1\n")
  cat(str(check_CUX1))
  cat("\n")
  
  #### READ Linked_peak_to_selected_genes ----
  
  Linked_peak_to_selected_genes<-readRDS(file=opt$Linked_peak_to_selected_genes)
  
  
  cat("Linked_peak_to_selected_genes\n")
  cat(str(Linked_peak_to_selected_genes))
  cat("\n")
  
  
  colnames(Linked_peak_to_selected_genes)[which(colnames(Linked_peak_to_selected_genes) == 'peak')]<-'Peak_ID'
  colnames(Linked_peak_to_selected_genes)[which(colnames(Linked_peak_to_selected_genes) == 'gene')]<-'Symbol'
  Linked_peak_to_selected_genes$Minus_logpval<-round(-1*log10(Linked_peak_to_selected_genes$pvalue),2)
  
  cat("Linked_peak_to_selected_genes_0\n")
  cat(str(Linked_peak_to_selected_genes))
  cat("\n")
  cat(str(unique(Linked_peak_to_selected_genes$Symbol)))
  cat("\n")
  cat(str(unique(Linked_peak_to_selected_genes$Peak_ID)))
  cat("\n")
  
  gr_Links <- GRanges(
    seqnames = as.character(gsub("^chr","",Linked_peak_to_selected_genes$seqnames)),
    name2=as.character(Linked_peak_to_selected_genes$Peak_ID),
    name3=as.numeric(Linked_peak_to_selected_genes$zscore),
    name4=as.numeric(Linked_peak_to_selected_genes$Minus_logpval),
    ranges=IRanges(
      start=as.numeric(Linked_peak_to_selected_genes$start),
      end=as.numeric(Linked_peak_to_selected_genes$end),
      name=paste('Link',as.character(Linked_peak_to_selected_genes$seqnames),
                 as.character(Linked_peak_to_selected_genes$start),
                 as.character(Linked_peak_to_selected_genes$end),                                   
                 sep='_')))
  
 
  
  #### READ Master_peak_file_with_SNP_numbered_peaks ----
  
  Master_peak_file_with_SNP_numbered_peaks<-readRDS(file=opt$Master_peak_file_with_SNP_numbered_peaks)
  
  
  cat("Master_peak_file_with_SNP_numbered_peaks\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks))
  cat("\n")
  
  # cat("gr_Links_0\n")
  # cat(str(gr_Links))
  # cat("\n")
  
  ALL_PoIs<-unique(Master_peak_file_with_SNP_numbered_peaks$Peak_ID)
  
  gr_ALL_PoIs <- GRanges(
    seqnames = as.character(gsub("^chr","",gsub("-.+$","",ALL_PoIs))),  
    ranges=IRanges(
      start=as.integer(gsub("-.+$","",gsub("^[^-]+-","",ALL_PoIs))), 
      end=as.integer(gsub("^[^-]+-[^-]+-","",ALL_PoIs)),
      name=ALL_PoIs))
  
  cat("gr_ALL_PoIs_0\n")
  cat(str(gr_ALL_PoIs))
  cat("\n")
  
  
  Master_peak_file_with_SNP_numbered_peaks_long<-unique(as.data.frame(cSplit(Master_peak_file_with_SNP_numbered_peaks,sep = ';', direction = "long",
                                                                             splitCols = "Symbol_string"),stringsAsFactors=F))
  
  
  colnames(Master_peak_file_with_SNP_numbered_peaks_long)[which(colnames(Master_peak_file_with_SNP_numbered_peaks_long) == 'Symbol_string')]<-'Symbol'
  
  cat("Master_peak_file_with_SNP_numbered_peaks_long_0\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks_long))
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
  
  check_chr7<-K562_Regulatory_Build[which(K562_Regulatory_Build$seqid == '7'),]
  
  
  cat("check_chr7_0\n")
  cat(str(check_chr7))
  cat("\n")
  
  
  
  #### Read the ensembl_gtf ----
  
  ensembl_gtf<-readGFF(file=opt$ensembl_gtf)
  
  
  cat("ensembl_gtf_0\n")
  cat(str(ensembl_gtf))
  cat("\n")
  cat(sprintf(as.character(names(summary((as.factor(ensembl_gtf$type)))))))
  cat("\n")
  cat(sprintf(as.character(summary((as.factor(ensembl_gtf$type))))))
  cat("\n")
  
  
  ensembl_gtf_gene<-unique(ensembl_gtf[which(ensembl_gtf$type == 'gene'),])
  
  cat("ensembl_gtf_gene_0\n")
  cat(str(ensembl_gtf_gene))
  cat("\n")
  cat(str(unique(ensembl_gtf_gene$gene_id)))
  cat("\n")
  
 
  gr_gene <- GRanges(
    seqnames = as.character(ensembl_gtf_gene$seqid),
    name2=as.character(ensembl_gtf_gene$gene_name),
    name3=as.character(ensembl_gtf_gene$gene_biotype),
    strand=ensembl_gtf_gene$strand,
    ranges=IRanges(
      start=as.numeric(ensembl_gtf_gene$start),
      end=as.numeric(ensembl_gtf_gene$end),
      name=ensembl_gtf_gene$gene_id))
  
  gene_df <- data.frame(chr=as.character(seqnames(gr_gene)),
                        start=as.integer(start(gr_gene)),
                        end=as.integer(end(gr_gene)),
                        gene_name=as.character(gr_gene$name2),
                        gene_biotype=as.character(gr_gene$name3),
                        strand=strand(gr_gene),
                        feature_type='gene',
                        activity=NA,
                        Peak_name=NA,                                          
                        zscore=NA,
                        Minus_logpval=NA,
                        stable_id=names(gr_gene), stringsAsFactors = F)
  
  
  cat("gene_df_0\n")
  cat(str(gene_df))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor((gene_df$strand)))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor((gene_df$strand))))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor((gene_df$gene_biotype)))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor((gene_df$gene_biotype))))))
  cat("\n")
 
  
  ##### LOOP to populate regions --------------
  
  DEBUG<-0
  
  
  
  IR_df<-data.frame()
  Results_Genomic_Features<-data.frame()
  
  START<-1
  
  for(i in START:dim(Input_regions)[1])
  {
    Input_regions_sel<-Input_regions[i,]
    
    
    Symbol_sel<-Input_regions_sel$Symbol
    
    
    if(Symbol_sel == 'CUX1')
    {
      DEBUG<-1
    }else{
      DEBUG<-0
    }
    
    Display_genes_sel<-unique(unlist(strsplit(Input_regions_sel$Display_genes, split=";")))
    
    Region_sel<-Input_regions_sel$region_name
    chr_sel<-Input_regions_sel$chr
    start_sel<-Input_regions_sel$start
    end_sel<-Input_regions_sel$end
    FLAG_INFINITE<-sum(is.infinite(start_sel))
    
    cat("-------------------------------------->\t")
    cat(sprintf(as.character(paste(i,Region_sel,chr_sel,start_sel,end_sel,Symbol_sel,paste(Display_genes_sel, collapse=';'),FLAG_INFINITE,collapse="  "))))
    cat("\n")
    
    if(FLAG_INFINITE == 0){
      
      if(Symbol_sel == 'ABCA1')
      {
        DEBUG<-1

      }else{

        DEBUG<-0
      }
      
      if(DEBUG == 1)
      {
        cat("Input_regions_sel_0\n")
        cat(str(Input_regions_sel))
        cat("\n")
      }
      
      
      Peaks_to_display_df_long_sel<-unique(Peaks_to_display_df_long[which(Peaks_to_display_df_long$Symbol == Symbol_sel),])
      
      if(DEBUG == 1)
      {
        cat("Peaks_to_display_df_long_sel_0\n")
        cat(str(Peaks_to_display_df_long_sel))
        cat("\n")
      }
      
      if(dim(Peaks_to_display_df_long_sel)[1] >0){
        
        string_of_Peaks_to_display<-paste(unique(Peaks_to_display_df_long_sel$Peaks_to_display),collapse='|')
        
        
      }else{
        
        string_of_Peaks_to_display<-NA
        
      }#dim(Peaks_to_display_df_long_sel)[1] >0
      
      
      if(DEBUG == 1)
      {
        cat("string_of_Peaks_to_display_0\n")
        cat(str(string_of_Peaks_to_display))
        cat("\n")
      }
      
      ### Sel TSS Peaks in my Input regions ----
      
      Master_peak_file_with_SNP_numbered_peaks_long_sel<-Master_peak_file_with_SNP_numbered_peaks_long[which(Master_peak_file_with_SNP_numbered_peaks_long$Peak_CLASS == 'TSS_Peak' &
                                                                                                               Master_peak_file_with_SNP_numbered_peaks_long$Symbol%in%Display_genes_sel),]
      
      if(DEBUG == 1)
      {
        cat("Master_peak_file_with_SNP_numbered_peaks_long_sel_0\n")
        cat(str(Master_peak_file_with_SNP_numbered_peaks_long_sel))
        cat("\n")
        cat(str(unique(Master_peak_file_with_SNP_numbered_peaks_long_sel$Peak_ID)))
        cat("\n")
      }
      
      if(dim(Master_peak_file_with_SNP_numbered_peaks_long_sel)[1] >0){
        
        string_of_TSS_peaks_to_display<-unique(Master_peak_file_with_SNP_numbered_peaks_long_sel$Peak_Number)
        
      }else{
        
        string_of_TSS_peaks_to_display<-NA
        
      }#dim(Master_peak_file_with_SNP_numbered_peaks_long_sel)[1] >0)
      
      
      if(DEBUG == 1)
      {
        cat("string_of_TSS_peaks_to_display_0\n")
        cat(str(string_of_TSS_peaks_to_display))
        cat("\n")
      }
      
      
      bundle_of_selected_peaks<-paste(unique(c(string_of_TSS_peaks_to_display,string_of_Peaks_to_display)), collapse ='|')
      
      if(DEBUG == 1)
      {
        cat("bundle_of_selected_peaks_0\n")
        cat(str(bundle_of_selected_peaks))
        cat("\n")
      }
      
      
      Input_regions_sel$Display_peaks<-bundle_of_selected_peaks
      
      
      if(DEBUG == 1)
      {
        cat("Input_regions_sel_0\n")
        cat(str(Input_regions_sel))
        cat("\n")
      }
      
      IR_df<-rbind(Input_regions_sel,IR_df)
      
      ##### GR object for the region ------------
      
      gr_region <- GRanges(
        seqnames = as.character(gsub("^chr","",chr_sel)),    
        ranges=IRanges(
          start=as.numeric(start_sel),
          end=as.numeric(end_sel),
          name=Region_sel))
      
      if(DEBUG == 1)
      {
        cat("gr_region_0\n")
        cat(str(gr_region))
        cat("\n")
        
      }
      
      
      ##### Find overlaps with Peaks ------------
      
      
      m <- findOverlaps(gr_ALL_PoIs,gr_region,
                        ignore.strand = TRUE)
      
      if(DEBUG == 1)
      {
        cat("m\n")
        cat(str(m))
        cat("\n")
      }
      
      subjectHits_region<-subjectHits(m)
      
      if(DEBUG == 1)
      {
        cat("subjectHits_region\n")
        cat(str(subjectHits_region))
        cat("\n")
      }
      
      queryHits_Peaks<-queryHits(m)
      
      if(DEBUG == 1)
      {
        cat("queryHits_Peaks\n")
        cat(str(queryHits_Peaks))
        cat("\n")
      }
      
      Peaks_df <- data.frame(chr=as.character(seqnames(gr_ALL_PoIs)),
                             start=as.integer(start(gr_ALL_PoIs)),
                             end=as.integer(end(gr_ALL_PoIs)),
                             feature_type='Peak',
                             strand='+',
                             activity=NA,
                             gene_name=NA,
                             gene_biotype=NA,
                             Peak_name=NA,                                          
                             zscore=NA,
                             Minus_logpval=NA,
                             stable_id=names(gr_ALL_PoIs), stringsAsFactors = F)
      
      if(DEBUG == 1)
      {
        cat("Peaks_df_0\n")
        cat(str(Peaks_df))
        cat("\n")
      }
      
      Peaks_df_hits<-Peaks_df[queryHits_Peaks,]
      
      if(DEBUG == 1)
      {
        cat("Peaks_df_hits_0\n")
        cat(str(Peaks_df_hits))
        cat("\n")
      }
      
      
      
      if(dim(Peaks_df_hits)[1] >0)
      {
        if(DEBUG == 1)
        {
          cat("Peaks_df_hits_0\n")
          cat(str(Peaks_df_hits))
          cat("\n")
        }
        
        Peaks_df_hits<-unique(Peaks_df_hits)
        
        Peaks_df_hits$region_string<-Region_sel
        Peaks_df_hits$Symbol<-Symbol_sel
        
        if(DEBUG == 1)
        {
          cat("Peaks_df_hits_1\n")
          cat(str(Peaks_df_hits))
          cat("\n")
        }
        
      }else{
        
        Peaks_df_hits<-data.frame()
      }#dim(Peaks_df_hits)[1] >0,
      
      
      ##### Find overlaps with Links ------------
      
      
      m <- findOverlaps(gr_Links,gr_region,
                        ignore.strand = TRUE)
      
      if(DEBUG == 1)
      {
        cat("m\n")
        cat(str(m))
        cat("\n")
      }
      
      subjectHits_region<-subjectHits(m)
      
      if(DEBUG == 1)
      {
        cat("subjectHits_region\n")
        cat(str(subjectHits_region))
        cat("\n")
      }
      
      queryHits_gene<-queryHits(m)
      
      if(DEBUG == 1)
      {
        cat("queryHits_gene\n")
        cat(str(queryHits_gene))
        cat("\n")
      }
      
      Links_df <- data.frame(chr=as.character(seqnames(gr_Links)),
                             start=as.integer(start(gr_Links)),
                             end=as.integer(end(gr_Links)),
                             Peak_name=as.character(gr_Links$name2),
                             strand='+',
                             feature_type='Link',
                             gene_name=NA,
                             gene_biotype=NA,
                             activity=NA,
                             zscore=as.numeric(gr_Links$name3),
                             Minus_logpval=as.numeric(gr_Links$name4),
                             stable_id=names(gr_Links), stringsAsFactors = F)
      
      if(DEBUG == 1)
      {
        cat("Links_df_0\n")
        cat(str(Links_df))
        cat("\n")
      }
      
      Links_df_hits<-Links_df[queryHits_gene,]
      
      if(DEBUG == 1)
      {
        cat("Links_df_hits_0\n")
        cat(str(Links_df_hits))
        cat("\n")
      }
      
      
      
      if(dim(Links_df_hits)[1] >0)
      {
        if(DEBUG == 1)
        {
          cat("Links_df_hits_0\n")
          cat(str(Links_df_hits))
          cat("\n")
        }
        
        Links_df_hits<-unique(Links_df_hits)
        
        Links_df_hits$region_string<-Region_sel
        Links_df_hits$Symbol<-Symbol_sel
        
        if(DEBUG == 1)
        {
          cat("Links_df_hits_1\n")
          cat(str(Links_df_hits))
          cat("\n")
        }
        
      }else{
        
        Links_df_hits<-data.frame()
      }#dim(Links_df_hits)[1] >0
      
      
      ##### Find overlaps with genes ------------
      
      m <- findOverlaps(gr_gene,gr_region,
                        ignore.strand = TRUE)
      
      if(DEBUG == 1)
      {
        cat("m\n")
        cat(str(m))
        cat("\n")
      }
      
      subjectHits_region<-subjectHits(m)
      
      if(DEBUG == 1)
      {
        cat("subjectHits_region\n")
        cat(str(subjectHits_region))
        cat("\n")
      }
      
      queryHits_gene<-queryHits(m)
      
      if(DEBUG == 1)
      {
        cat("queryHits_gene\n")
        cat(str(queryHits_gene))
        cat("\n")
      }
      
      gene_df <- data.frame(chr=as.character(seqnames(gr_gene)),
                            start=as.integer(start(gr_gene)),
                            end=as.integer(end(gr_gene)),
                            gene_name=as.character(gr_gene$name2),
                            gene_biotype=as.character(gr_gene$name3),
                            strand=strand(gr_gene),
                            feature_type='gene',
                            activity=NA,
                            Peak_name=NA,                                          
                            zscore=NA,
                            Minus_logpval=NA,
                            stable_id=names(gr_gene), stringsAsFactors = F)
      
      if(DEBUG == 1)
      {
        cat("gene_df_0\n")
        cat(str(gene_df))
        cat("\n")
        cat(sprintf(as.character(names(summary(as.factor((gene_df$strand)))))))
        cat("\n")
        cat(sprintf(as.character(summary(as.factor((gene_df$strand))))))
        cat("\n")
        cat(sprintf(as.character(names(summary(as.factor((gene_df$gene_biotype)))))))
        cat("\n")
        cat(sprintf(as.character(summary(as.factor((gene_df$gene_biotype))))))
        cat("\n")
      }
      
      gene_df_hits<-gene_df[queryHits_gene,]
      
      if(DEBUG == 1)
      {
        cat("gene_df_hits_0\n")
        cat(str(gene_df_hits))
        cat("\n")
      }
      
      
      
      if(dim(gene_df_hits)[1] >0)
      {
        if(DEBUG == 1)
        {
          cat("gene_df_hits_0\n")
          cat(str(gene_df_hits))
          cat("\n")
        }
        
        gene_df_hits<-unique(gene_df_hits)
        
        gene_df_hits$region_string<-Region_sel
        gene_df_hits$Symbol<-Symbol_sel
        
        
        if(DEBUG == 1)
        {
          cat("gene_df_hits_1\n")
          cat(str(gene_df_hits))
          cat("\n")
          cat(sprintf(as.character(names(summary(as.factor(gene_df_hits$gene_name))))))
          cat("\n")
          cat(sprintf(as.character(summary(as.factor(gene_df_hits$gene_name)))))
          cat("\n")
        }
        
      }else{
        
        gene_df_hits<-data.frame()
      }#dim(gene_df_hits)[1] >0
      
      
      
      ##### Find overlaps with K-562 regulatory build ------------
      
      m <- findOverlaps(gr_K562_Regulatory_Build,gr_region,
                        ignore.strand = TRUE)
      
      if(DEBUG == 1)
      {
        cat("m\n")
        cat(str(m))
        cat("\n")
      }
      
      subjectHits_region<-subjectHits(m)
      
      if(DEBUG == 1)
      {
        cat("subjectHits_region\n")
        cat(str(subjectHits_region))
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
                                             strand='+',
                                             gene_biotype=NA,
                                             activity=as.character(gr_K562_Regulatory_Build$name3),
                                             gene_name=NA,
                                             Peak_name=NA,                                          
                                             zscore=NA,
                                             Minus_logpval=NA,
                                             stable_id=names(gr_K562_Regulatory_Build), stringsAsFactors = F)
      
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
      
      
      
      if(dim(K562_Regulatory_Build_df_hits)[1] >0)
      {
        if(DEBUG == 1)
        {
          cat("K562_Regulatory_Build_df_hits_0\n")
          cat(str(K562_Regulatory_Build_df_hits))
          cat("\n")
        }
        
        K562_Regulatory_Build_df_hits<-unique(K562_Regulatory_Build_df_hits)
        
        K562_Regulatory_Build_df_hits$region_string<-Region_sel
        K562_Regulatory_Build_df_hits$Symbol<-Symbol_sel
        
        if(DEBUG == 1)
        {
          cat("K562_Regulatory_Build_df_hits_1\n")
          cat(str(K562_Regulatory_Build_df_hits))
          cat("\n")
        }
        
      }else{
        
        K562_Regulatory_Build_df_hits<-data.frame()
      }#dim(K562_Regulatory_Build_df_hits)[1] >0
      
      
      
      
      
      
      
      DEF<-unique(rbind(gene_df_hits,K562_Regulatory_Build_df_hits))
      
      if(DEBUG == 1)
      {
        cat("DEF_0\n")
        cat(str(DEF))
        cat("\n")
      }
      
      if(DEBUG == 1)
      {
        cat("Peaks_df_hits_REMEMBER\n")
        cat(str(Peaks_df_hits))
        cat("\n")
      }
      
      DEF<-rbind(DEF,Peaks_df_hits)
      
      if(DEBUG == 1)
      {
        cat("DEF_1\n")
        cat(str(DEF))
        cat("\n")
      }
      
      DEF<-rbind(DEF,Links_df_hits)
      
      if(DEBUG == 1)
      {
        cat("DEF_2\n")
        cat(str(DEF))
        cat("\n")
      }
      
      DEF$region_name<-Region_sel
      
      
      if(dim(DEF)[1] >0){
        
        Results_Genomic_Features<-unique(rbind(DEF,Results_Genomic_Features))
        
        if(DEBUG == 1)
        {
          cat("Results_Genomic_Features_0\n")
          cat(str(Results_Genomic_Features))
          cat("\n")
        }
      }#dim(DEF)[1] >0
      
      
    }else{
      
      cat("-->Hello_world\n")
      
    }#FLAG_INFINITE == 0
    

  }# i in 1:dim(Input_regions)[1]
  
  
  ############# SAVE Input_Regions_DEF ------------------
  
  if(dim(IR_df)[1])
  {
    cat("IR_df_FINAL\n")
    cat(str(IR_df))
    cat("\n")
    
    check_CUX1<-IR_df[which(IR_df$Symbol == 'CUX1'),]
    
    cat("check_CUX1_FINAL\n")
    cat(str(check_CUX1))
    cat("\n")
    
    setwd(out)
    
    saveRDS(IR_df, file="Input_Regions_DEF.rds")
    write.table(IR_df, file="Input_Regions_DEF.tsv", sep="\t",quote = F, row.names = F)
  }
  
  ############# SAVE Results_Genomic_Features ------------------
  
  if(dim(Results_Genomic_Features)[1])
  {
    cat("Results_Genomic_Features_FINAL\n")
    cat(str(Results_Genomic_Features))
    cat("\n")
    
    Results_Genomic_Features$feature_type[which(Results_Genomic_Features$feature_type == 'CTCF Binding Site')]<-'CTCF'
    Results_Genomic_Features$feature_type[which(Results_Genomic_Features$feature_type == 'Open chromatin')]<-'OpenChromatin'
    Results_Genomic_Features$feature_type[which(Results_Genomic_Features$feature_type == 'TF binding')]<-'TFBS'
    
    cat(sprintf(as.character(names(summary(as.factor(Results_Genomic_Features$feature))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Results_Genomic_Features$feature)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(Results_Genomic_Features$activity))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(Results_Genomic_Features$activity)))))
    cat("\n")
    
    
    Results_Genomic_Features$feature_type<-factor(Results_Genomic_Features$feature_type,
                                                  levels=c('gene','Promoter','Enhancer','TFBS','OpenChromatin','CTCF','Peak','Link'),
                                                  ordered=T)
    
    Results_Genomic_Features$activity<-factor(Results_Genomic_Features$activity,
                                              levels=c("ACTIVE",'INACTIVE','POISED',"REPRESSED"),
                                              ordered=T)
    
    cat(sprintf(as.character(names(summary(Results_Genomic_Features$feature)))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_Genomic_Features$feature))))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results_Genomic_Features$activity)))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_Genomic_Features$activity))))
    cat("\n")
    
    indx.features<-which(Results_Genomic_Features$feature_type%in%c('Promoter','Enhancer','TFBS','OpenChromatin','CTCF'))
    indx.Peaks<-which(Results_Genomic_Features$feature_type%in%c('Peak'))
    indx.Links<-which(Results_Genomic_Features$feature_type%in%c('Link'))
    
    Results_Genomic_Features$interaction<-NA
    
    Results_Genomic_Features$interaction[indx.features]<-paste(Results_Genomic_Features$feature_type[indx.features],
                                                           Results_Genomic_Features$activity[indx.features],
                                                           sep='|')
    
    
    Results_Genomic_Features$interaction<-factor(Results_Genomic_Features$interaction,
                                             levels=c('Promoter|ACTIVE','Promoter|POISED','Promoter|REPRESSED','Promoter|INACTIVE',
                                                      'Enhancer|ACTIVE','Enhancer|POISED','Enhancer|REPRESSED','Enhancer|INACTIVE',
                                                      'TFBS|ACTIVE','TFBS|POISED','TFBS|REPRESSED','TFBS|INACTIVE',
                                                      'OpenChromatin|ACTIVE','OpenChromatin|POISED','OpenChromatin|REPRESSED','OpenChromatin|INACTIVE',
                                                      'CTCF|ACTIVE','CTCF|POISED','CTCF|REPRESSED','CTCF|INACTIVE'),
                                             ordered=T)
    
    
    cat("----------------------------->INTERACTION\n")
    cat(sprintf(as.character(names(summary(Results_Genomic_Features$interaction)))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_Genomic_Features$interaction))))
    cat("\n")
    
    
    
    
    Results_Genomic_Features$ymax<-NA
    Results_Genomic_Features$ymin<-NA
    
    Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'gene' & Results_Genomic_Features$strand == '+')]<-4
    Results_Genomic_Features$ymax[which(Results_Genomic_Features$feature_type == 'gene' & Results_Genomic_Features$strand == '+')]<-Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'gene' & Results_Genomic_Features$strand == '+')]+0.2
    Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'gene' & Results_Genomic_Features$strand == '-')]<-3
    Results_Genomic_Features$ymax[which(Results_Genomic_Features$feature_type == 'gene' & Results_Genomic_Features$strand == '-')]<-Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'gene' & Results_Genomic_Features$strand == '-')]+0.2
    
    
    
    indx.Reg<-which(Results_Genomic_Features$feature_type%in%c('Promoter','Enhancer','TFBS','OpenChromatin','CTCF'))
    
    cat("indx.Reg_0\n")
    cat(str(indx.Reg))
    cat("\n")
    
    random_vec <- round(runif(n=length(indx.Reg), min=-0.25, max=0.25),2)
    
    cat("random_vec_0\n")
    cat(str(random_vec))
    cat("\n")
    
    Results_Genomic_Features$ymin[indx.Reg]<-2
    Results_Genomic_Features$ymax[indx.Reg]<-Results_Genomic_Features$ymin[indx.Reg]+0.2
    
    Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'Peak')]<-1
    Results_Genomic_Features$ymax[which(Results_Genomic_Features$feature_type == 'Peak')]<-1.2
    Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'Link')]<-0.5
    Results_Genomic_Features$ymax[which(Results_Genomic_Features$feature_type == 'Link')]<-1
    
    cat(sprintf(as.character(names(summary(Results_Genomic_Features$ymin[indx.Reg])))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_Genomic_Features$ymin[indx.Reg]))))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results_Genomic_Features$ymax[indx.Reg])))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_Genomic_Features$ymax[indx.Reg]))))
    cat("\n")
    
    cat(sprintf(as.character(names(summary(Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'gene')])))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'gene')]))))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results_Genomic_Features$ymax[which(Results_Genomic_Features$feature_type == 'gene')])))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_Genomic_Features$ymax[which(Results_Genomic_Features$feature_type == 'gene')]))))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'Peak')])))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'Peak')]))))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results_Genomic_Features$ymax[which(Results_Genomic_Features$feature_type == 'Peak')])))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_Genomic_Features$ymax[which(Results_Genomic_Features$feature_type == 'Peak')]))))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'Link')])))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_Genomic_Features$ymin[which(Results_Genomic_Features$feature_type == 'Link')]))))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results_Genomic_Features$ymax[which(Results_Genomic_Features$feature_type == 'Link')])))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_Genomic_Features$ymax[which(Results_Genomic_Features$feature_type == 'Link')]))))
    cat("\n")
    
    
    setwd(out)
    
    saveRDS(Results_Genomic_Features, file="ALL_regions_Genomic_features.rds")
    write.table(Results_Genomic_Features, file="ALL_regions_Genomic_features.tsv", sep="\t",quote = F, row.names = F)
    
    
  }#dim(Results_Genomic_Features)[1]
  
  
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
    make_option(c("--Peaks_to_display"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Input_regions"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Master_peak_file_with_SNP_numbered_peaks"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Linked_peak_to_selected_genes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ensembl_gtf"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--K562_Regulatory_Build"), type="character", default=NULL, 
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
  
  
  Build_input_regions(opt)
  

  
}


###########################################################################

system.time( main() )