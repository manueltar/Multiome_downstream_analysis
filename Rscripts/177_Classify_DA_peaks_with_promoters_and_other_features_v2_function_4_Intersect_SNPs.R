
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


Intersect_SNPS = function(option_list)
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
  
  #### READ and transform selected_variants ----
  
  selected_variants = unlist(strsplit(opt$selected_variants, split=","))
  
  cat("selected_variants_0\n")
  cat(sprintf(as.character(selected_variants)))
  cat("\n")
  
  rs_vector<-gsub("__.+$","",selected_variants)
  VAR_vector<-gsub("[^__]+__","",selected_variants)
  
  cat("rs_vector_0\n")
  cat(str(rs_vector))
  cat("\n")
  
  cat("VAR_vector_0\n")
  cat(str(VAR_vector))
  cat("\n")
  
  
  selected_variants_df<-as.data.frame(cbind(rs_vector,VAR_vector))
  
  colnames(selected_variants_df)<-c('rs','VAR')
  
  selected_variants_df$chr<-gsub("_.+$","",selected_variants_df$VAR)
  selected_variants_df$pos38<-gsub("^[^_]+_","",selected_variants_df$VAR)
  selected_variants_df$pos38<-as.integer(gsub("_.+$","",selected_variants_df$pos38))
  selected_variants_df$ref<-gsub("^[^_]+_[^_]+_","",selected_variants_df$VAR)
  selected_variants_df$ref<-gsub("_.+$","",selected_variants_df$ref)
  selected_variants_df$alt<-gsub("^[^_]+_[^_]+_[^_]+_","",selected_variants_df$VAR)
  
  cat("selected_variants_df_0\n")
  cat(str(selected_variants_df))
  cat("\n")
  
  gr_VARS <- GRanges(
    seqnames = as.character(gsub("chr","",selected_variants_df$chr)),
    name2=as.character(selected_variants_df$rs),
    ranges=IRanges(
      start=as.numeric(selected_variants_df$pos38),
      end=as.numeric(selected_variants_df$pos38),
      name=selected_variants_df$VAR))
  
  cat("gr_VARS_0\n")
  cat(str(gr_VARS))
  cat("\n")
  
  #### Read linked and classified peaks----
  
  
  Master_peak_file<-readRDS(file=opt$Master_peak_file)
  
  Master_peak_file$chr<-paste('chr',Master_peak_file$chr,sep='')
  
  cat("Master_peak_file_0\n")
  cat(str(Master_peak_file))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Master_peak_file$chr))))))
  cat("\n")
  cat(sprintf(as.character(names(as.factor(Master_peak_file$chr)))))
  cat("\n")
  
  gr_peaks <- GRanges(
    seqnames = as.character(gsub("chr","",Master_peak_file$chr)),
    ranges=IRanges(
      start=as.numeric(Master_peak_file$start),
      end=as.numeric(Master_peak_file$end),
      name=Master_peak_file$Peak_ID))
  
  cat("gr_peaks_0\n")
  cat(str(gr_peaks))
  cat("\n")
  
  
  #### Intersect Peak to genes with SNP ----
  
  DEBUG <- 1
  
  m <- findOverlaps(gr_VARS,gr_peaks)
  
  if(DEBUG == 1)
  {
    cat("m\n")
    cat(str(m))
    cat("\n")
  }
  
  subjectHits_peaks<-subjectHits(m)
  
  if(DEBUG == 1)
  {
    cat("subjectHits_peaks\n")
    cat(str(subjectHits_peaks))
    cat("\n")
  }
  
  queryHits_VARS<-queryHits(m)
  
  if(DEBUG == 1)
  {
    cat("queryHits_VARS\n")
    cat(str(queryHits_VARS))
    cat("\n")
  }
  
  VARS_df <- data.frame(rs=as.character(gr_VARS$name2),
                        VAR=names(gr_VARS), stringsAsFactors = F)
  
  if(DEBUG == 1)
  {
    cat("VARS_df_0\n")
    cat(str(VARS_df))
    cat("\n")
  }
  
  VARS_df_hits<-VARS_df[queryHits_VARS,]
  
  if(DEBUG == 1)
  {
    cat("VARS_df_hits_0\n")
    cat(str(VARS_df_hits))
    cat("\n")
  }
  
  peaks_df_SNP <- data.frame(chr=as.character(seqnames(gr_peaks)),
                                     start=as.integer(start(gr_peaks)),
                                     end=as.integer(end(gr_peaks)),
                                     Peak_ID=names(gr_peaks),
                                     stringsAsFactors = F)
  
  
  if(DEBUG == 1)
  {
    cat("peaks_df_SNP_0\n")
    cat(str(peaks_df_SNP))
    cat("\n")
  }
  
  peaks_df_SNP_hits<-peaks_df_SNP[subjectHits_peaks,]
  
  if(dim(peaks_df_SNP_hits)[1] >0)
  {
    if(DEBUG == 1)
    {
      cat("peaks_df_SNP_hits_0\n")
      cat(str(peaks_df_SNP_hits))
      cat("\n")
    }
    
    peaks_df_SNP_hits<-cbind(peaks_df_SNP_hits,VARS_df_hits)
    
    peaks_df_SNP_hits$chr<-paste('chr',peaks_df_SNP_hits$chr,sep='')
    
    if(DEBUG == 1)
    {
      cat("peaks_df_SNP_hits_1\n")
      cat(str(peaks_df_SNP_hits))
      cat("\n")
    }
    
 
    Master_peak_file<-merge(Master_peak_file,
                                       peaks_df_SNP_hits,
                                       by=c('chr','start','end','Peak_ID'),
                                       all.x=T)
    
    if(DEBUG == 1)
    {
      cat("Master_peak_file_POST_MERGE\n")
      cat(str(Master_peak_file))
      cat("\n")
      cat(str(unique(Master_peak_file$rs)))
      cat("\n")
    }
    
    
    #### SAVE RDS ----
    
    
    
    setwd(out)
    
    saveRDS(Master_peak_file,file='Master_peak_file_with_SNP.rds')
    write.table(Master_peak_file,file='Master_peak_file_with_SNP.tsv',sep="\t",quote=F,row.names = F)
    
  }#dim(peaks_df_SNP_hits)[1] >0
  
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
    make_option(c("--selected_variants"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Master_peak_file"), type="character", default=NULL, 
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
  
  
  Intersect_SNPS(opt)
  

  
}


###########################################################################

system.time( main() )