
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


Number_the_peaks = function(option_list)
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
  
  #### Read linked and classified peaks----
  
  
  Master_peak_file<-readRDS(file=opt$Master_peak_file)
  
  cat("Master_peak_file_0\n")
  cat(str(Master_peak_file))
  cat("\n")
  cat(str(Master_peak_file$Peak_ID))
  cat("\n")
  
  
  #### Unassigned peaks to genes -----
  
  Master_peak_file$Peak_Number<-NA
  
  indx.unassigned<-which(Master_peak_file$Peak_CLASS == 'Non_gene_associated_peak')
  
  cat("indx.unassigned_0\n")
  cat(str(indx.unassigned))
  cat("\n")
  
  indx.no_gene_found<-which(Master_peak_file$Symbol_string == '')
  
  
  cat("indx.no_gene_found_0\n")
  cat(str(indx.no_gene_found))
  cat("\n")
  
  Master_peak_file_identified<-Master_peak_file[-c(indx.unassigned,indx.no_gene_found),]
  
  cat("Master_peak_file_identified_0\n")
  cat(str(Master_peak_file_identified))
  cat("\n")
  cat(str(Master_peak_file_identified$Peak_ID))
  cat("\n")
  
  Master_peak_file_not_identified<-Master_peak_file[c(indx.unassigned,indx.no_gene_found),]
  
  Master_peak_file_not_identified$Total_peaks<-NA
  
  cat("Master_peak_file_not_identified_0\n")
  cat(str(Master_peak_file_not_identified))
  cat("\n")
  cat(str(Master_peak_file_not_identified$Peak_ID))
  cat("\n")
 
  
  #### Number the peaks by Peak_CLASS and Symbol_string ----
  
  Master_peak_file_identified.dt<-data.table(Master_peak_file_identified,
                                  key=c("Symbol_string","Peak_CLASS"))
  
  
  Master_peak_file_identified_summary<-as.data.frame(Master_peak_file_identified.dt[,.(Total_peaks=.N), by=key(Master_peak_file_identified.dt)], stringsAsFactors=F)
  
  
  cat("Master_peak_file_identified_summary\n")
  cat(str(Master_peak_file_identified_summary))
  cat("\n")
  cat(sprintf(as.character(names(summary(Master_peak_file_identified_summary$Total_peaks)))))
  cat("\n")
  cat(sprintf(as.character(summary(Master_peak_file_identified_summary$Total_peaks))))
  cat("\n")
  
  
  Master_peak_file_identified<-merge(Master_peak_file_identified,
                                     Master_peak_file_identified_summary,
                                     by=c("Symbol_string","Peak_CLASS"))
  
  cat("Master_peak_file_identified_1\n")
  cat(str(Master_peak_file_identified))
  cat("\n")
  cat(str(Master_peak_file_identified$Peak_ID))
  cat("\n")
  
  
  #### LOOP to assign Peaks ----
  
  
  array_Peak_CLASS<-levels(droplevels(Master_peak_file_identified$Peak_CLASS))
  
  cat("array_Peak_CLASS_0\n")
  cat(str(array_Peak_CLASS))
  cat("\n")
  
  DEBUG<-0
  
  
  Numbered_df<-data.frame()
  
  for(i in 1:length(array_Peak_CLASS))
  {
    array_Peak_CLASS_sel<-array_Peak_CLASS[i]
    
    cat("-------------------------------------------->\t")
    cat(sprintf(as.character(array_Peak_CLASS_sel)))
    cat("\n")
    
    
    Master_peak_file_identified_sel<-Master_peak_file_identified[which(Master_peak_file_identified$Peak_CLASS == array_Peak_CLASS_sel),]
    
    
    if(DEBUG == 1)
    {
      cat("Master_peak_file_identified_sel_0\n")
      cat(str(Master_peak_file_identified_sel))
      cat("\n")
    }
    
    array_Symbol_string<-unique(Master_peak_file_identified_sel$Symbol_string)
    
    if(DEBUG == 1)
    {
      cat("array_Symbol_string_0\n")
      cat(str(array_Symbol_string))
      cat("\n")
    }
    
    for(k in 1:length(array_Symbol_string))
    {
      
      array_Symbol_string_sel<-array_Symbol_string[k]
      
      cat("-->\t")
      cat(sprintf(as.character(array_Symbol_string_sel)))
      cat("\n")
      
      
      Master_peak_file_identified_sel_Symbol_string_sel<-unique(Master_peak_file_identified_sel[which(Master_peak_file_identified_sel$Symbol_string == array_Symbol_string_sel),])
      
      Master_peak_file_identified_sel_Symbol_string_sel<-Master_peak_file_identified_sel_Symbol_string_sel[order(Master_peak_file_identified_sel_Symbol_string_sel$start, decreasing = FALSE),]
      
      if(DEBUG == 1)
      {
        cat("Master_peak_file_identified_sel_Symbol_string_sel_0\n")
        cat(str(Master_peak_file_identified_sel_Symbol_string_sel))
        cat("\n")
      }
      
      total_peaks<-unique(Master_peak_file_identified_sel_Symbol_string_sel$Total_peaks)
      
      if(DEBUG == 1)
      {
        cat("total_peaks_0\n")
        cat(str(total_peaks))
        cat("\n")
      }
      
      vector_of_Peak_numbers<-seq(1,total_peaks, by=1)
      
      if(DEBUG == 1)
      {
        cat("vector_of_Peak_numbers_0\n")
        cat(str(vector_of_Peak_numbers))
        cat("\n")
      }
      
      Master_peak_file_identified_sel_Symbol_string_sel$Peak_Number<-paste(Master_peak_file_identified_sel_Symbol_string_sel$Peak_CLASS,vector_of_Peak_numbers,
                                                                           Master_peak_file_identified_sel_Symbol_string_sel$Symbol_string, sep='__')
      
      
      if(DEBUG == 1)
      {
        cat("Master_peak_file_identified_sel_Symbol_string_sel_1\n")
        cat(str(Master_peak_file_identified_sel_Symbol_string_sel))
        cat("\n")
      }
      
      Numbered_df<-rbind(Master_peak_file_identified_sel_Symbol_string_sel,Numbered_df)
      
      
      
    }#k in 1:length(array_Symbol_string)
    
  }# i in 1:length(array_Peak_CLASS)
  
  #### Add together ----
  
  cat("Master_peak_file_not_identified_REMEMBER\n")
  cat(str(Master_peak_file_not_identified))
  cat("\n")
  cat(str(Master_peak_file_not_identified$Peak_ID))
  cat("\n")
  
  cat("Numbered_df_REMEMBER\n")
  cat(str(Numbered_df))
  cat("\n")
  cat(str(Numbered_df$Peak_ID))
  cat("\n")
  
  
  Master_peak_file_FINAL<-rbind(Numbered_df,Master_peak_file_not_identified)
  
  cat("Master_peak_file_FINAL_0\n")
  cat(str(Master_peak_file_FINAL))
  cat("\n")
  cat(str(Master_peak_file_FINAL$Peak_ID))
  cat("\n")
  
  #### SAVE RDS ----
  
  
  
  setwd(out)
  
  saveRDS(Master_peak_file_FINAL,file='Master_peak_file_with_SNP_numbered_peaks.rds')
  write.table(Master_peak_file_FINAL,file='Master_peak_file_with_SNP_numbered_peaks.tsv',sep="\t",quote=F,row.names = F)
  
  
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
  
  
  Number_the_peaks(opt)
  

  
}


###########################################################################

system.time( main() )