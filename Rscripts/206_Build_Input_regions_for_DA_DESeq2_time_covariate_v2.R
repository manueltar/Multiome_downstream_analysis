
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
  
  #### READ and transform ery_genes_array ----
  
  ery_genes_array = unlist(strsplit(opt$ery_genes_array, split=','))
  
  cat("ery_genes_array_\n")
  cat(sprintf(as.character(ery_genes_array)))
  cat("\n")
  
  #### READ and transform out ----
  
  myeloid_genes_array = unique(unlist(strsplit(opt$myeloid_genes_array, split=',')))
  
  cat("myeloid_genes_array_\n")
  cat(sprintf(as.character(myeloid_genes_array)))
  cat("\n")
  
  #### READ and transform megak_genes_array ----
  
  megak_genes_array = unique(unlist(strsplit(opt$megak_genes_array, split=',')))
  
  cat("megak_genes_array_\n")
  cat(sprintf(as.character(megak_genes_array)))
  cat("\n")
  
  #### READ and transform marker_genes_array ----
  
  marker_genes_array = unique(unlist(strsplit(opt$marker_genes_array, split=',')))
  
  cat("marker_genes_array_\n")
  cat(sprintf(as.character(marker_genes_array)))
  cat("\n")
  
  #### READ and transform ANAPC_genes ----
  
  ANAPC_genes = unique(unlist(strsplit(opt$ANAPC_genes, split=',')))
  
  cat("ANAPC_genes_\n")
  cat(sprintf(as.character(ANAPC_genes)))
  cat("\n")
  
  #### READ and transform DE_genes ----
  
  DE_genes<-readRDS(file=opt$DE_genes)
  
  cat("DE_genes_0\n")
  cat(str(DE_genes))
  cat("\n")
  
  DE_genes_SIG<-unique(DE_genes[which(DE_genes$Minus_logpval >= 1.3),])
  
  cat("DE_genes_SIG_0\n")
  cat(str(DE_genes_SIG))
  cat("\n")
  
  #### READ and transform Selected_genes_classified ----
  
  Selected_genes_classified<-readRDS(file=opt$Selected_genes_classified)
  
  cat("Selected_genes_classified_0\n")
  cat(str(Selected_genes_classified))
  cat("\n")
  cat(sprintf(as.character(names(summary(Selected_genes_classified$GENE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Selected_genes_classified$GENE_CLASS))))
  cat("\n")
  
  orginal_levels_of_Selected_genes_classified<-levels(Selected_genes_classified$GENE_CLASS)
  
  cat("orginal_levels_of_Selected_genes_classified_0\n")
  cat(sprintf(as.character(orginal_levels_of_Selected_genes_classified)))
  cat("\n")
  
  
  Selected_genes_classified$GENE_CLASS<-as.character(Selected_genes_classified$GENE_CLASS)
  
  Selected_genes_classified$GENE_CLASS[which(Selected_genes_classified$Symbol == 'CUX1')]<-'CUX1'
  
  Selected_genes_classified$GENE_CLASS[which(Selected_genes_classified$Symbol%in%marker_genes_array)]<-'Marker_genes'
  
  
  selection<-unique(c('CUX1','Marker_genes',orginal_levels_of_Selected_genes_classified))
  Selected_genes_classified_subset<-Selected_genes_classified[which(Selected_genes_classified$GENE_CLASS%in%selection),]
  
  cat("Selected_genes_classified_subset_0\n")
  cat(str(Selected_genes_classified_subset))
  cat("\n")
  
  levels_selected_genes<-selection
  
  cat("levels_selected_genes_0\n")
  cat(str(levels_selected_genes))
  cat("\n")
  
  #### Add missing genes ----
  
  
  new_genes_to_add_string<-unique(c(ery_genes_array,megak_genes_array,marker_genes_array,myeloid_genes_array,ANAPC_genes,EZH2_signature,CUX1,unique(DE_genes_SIG$Symbol)))
  
  cat("new_genes_to_add_string_0\n")
  cat(str(new_genes_to_add_string))
  cat("\n")
  
  
  
  new_genes_to_add <- data.frame(matrix(vector(), length(new_genes_to_add_string), dim(Selected_genes_classified_subset)[2],
                                        dimnames=list(c(),
                                                      colnames(Selected_genes_classified_subset))),stringsAsFactors=F)
  
  
  
  new_genes_to_add$Symbol<-new_genes_to_add_string
  new_genes_to_add$GENE_CLASS<-NA
  
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%ery_genes_array)]<-'Erythrocyte'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%megak_genes_array)]<-'Megakaryocyte'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%myeloid_genes_array)]<-'AML_myeloid'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%marker_genes_array)]<-'Marker_genes'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%ANAPC_genes)]<-'ANAPC_genes'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%EZH2_signature)]<-'EZH2_signature'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%CUX1)]<-'CUX1'
  
  new_genes_to_add$GENE_CLASS[is.na(new_genes_to_add$GENE_CLASS)]<-'Unclassified'
  
  
  
  cat("new_genes_to_add_0\n")
  cat(str(new_genes_to_add))
  cat("\n")
  cat(str(unique(new_genes_to_add$Symbol)))
  cat("\n")
  
  
  new_genes_to_add<-new_genes_to_add[-which(new_genes_to_add$Symbol%in%Selected_genes_classified_subset$Symbol),]
  
  cat("new_genes_to_add_After eliminating already present genes\n")
  cat(str(new_genes_to_add))
  cat("\n")
  cat(str(sprintf(as.character(new_genes_to_add$Symbol))))
  cat("\n")
  
  #### rbind ALL genes ----
  
  
  ALL_genes<-rbind(Selected_genes_classified_subset,
                   new_genes_to_add)
  
  cat("ALL_genes_0\n")
  cat(str(ALL_genes))
  cat("\n")
  cat(str(unique(ALL_genes$Symbol)))
  cat("\n")
  
  ALL_genes$GENE_CLASS<-droplevels(factor(ALL_genes$GENE_CLASS,
                                          levels=unique(c(levels_selected_genes,'CUX1','EZH2_signature','Marker_genes','CUX1_target_genes',
                                                          'MEP_progenitor','Megakaryocyte','Platelet','Platelet_volume',
                                                          'Erythrocyte','Mitosis','ANAPC_genes',
                                                          'PU1_target_genes','RUNX1_target_genes','GATA_target_genes',
                                                          'AML_myeloid','Unclassified')),
                                          ordered=T))
  
  cat("ALL_genes_1\n")
  cat(str(ALL_genes))
  cat("\n")
  cat(str(unique(ALL_genes$Symbol)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_genes$GENE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_genes$GENE_CLASS))))
  cat("\n")
  
  check.NA<-ALL_genes[is.na(ALL_genes$GENE_CLASS),]
  
  cat("-------------------------------------------------->check.NA_0\n")
  cat(str(check.NA))
  cat("\n")
  cat(str(unique(check.NA$Symbol)))
  cat("\n")
  cat(sprintf(as.character(names(summary(check.NA$GENE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(check.NA$GENE_CLASS))))
  cat("\n")
  
  ALL_genes$Display_genes<-ALL_genes$Symbol
  
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'BACH2')]<-'BACH2;MIR4464;ENSG00000260271'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'ITGA2B')]<-'ITGA2B;GPATCH8'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'RUNX1')]<-'RUNX1;LINC01436'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'TBXAS1')]<-'TBXAS1;HIPK2'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'CD109')]<-'CD109;CD109-AS1'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'WIPF1')]<-'WIPF1;CHN1'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'HBA1')]<-'HBA1;HBA2;HBQ1;HBZ;HBM'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'HBA2')]<-'HBA1;HBA2;HBQ1;HBZ;HBM'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'WIPF1')]<-'WIPF1;CHN1'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'KMT2E')]<-'KMT2E;KMT2E-AS1'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'LIFR')]<-'LIFR;MIR3650'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'SOX5')]<-'SOX5;SOX5-AS1'
  ALL_genes$Display_genes[which(ALL_genes$Symbol == 'MOCOS')]<-'MOCOS;COSMOC'
  
  
  cat("ALL_genes_1\n")
  cat(str(ALL_genes))
  cat("\n")
  cat(str(unique(ALL_genes$Symbol)))
  cat("\n")
  cat(sprintf(as.character(names(summary(ALL_genes$GENE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(ALL_genes$GENE_CLASS))))
  cat("\n")
  
  check.NA<-ALL_genes[is.na(ALL_genes$GENE_CLASS),]
  
  cat("-------------------------------------------------->check.NA_0\n")
  cat(str(check.NA))
  cat("\n")
  cat(str(unique(check.NA$Symbol)))
  cat("\n")
  cat(sprintf(as.character(names(summary(check.NA$GENE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(check.NA$GENE_CLASS))))
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
  
  ensembl_gtf_gene$region_START<-NA
  ensembl_gtf_gene$region_END<-NA
  
  
  indx.positive_strand<-which(ensembl_gtf_gene$strand == '+')
  
  cat("indx.positive_strand_0\n")
  cat(str(indx.positive_strand))
  cat("\n")
  indx.negative_strand<-which(ensembl_gtf_gene$strand == '-')
  
  cat("indx.negative_strand_0\n")
  cat(str(indx.negative_strand))
  cat("\n")
  
  
  #### Read linked and classified peaks----
  
  
  Master_peak_file<-readRDS(file=opt$Master_peak_file)
  

  cat("Master_peak_file_0\n")
  cat(str(Master_peak_file))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Master_peak_file$chr))))))
  cat("\n")
  cat(sprintf(as.character(names(as.factor(Master_peak_file$chr)))))
  cat("\n")
  
  Master_peak_file<-unique(as.data.frame(cSplit(Master_peak_file,sep = ';', direction = "long",
                                               splitCols = "Symbol_string"),stringsAsFactors=F))
  
  cat("Master_peak_file_1\n")
  cat(str(Master_peak_file))
  cat("\n")
 
  ##### LOOP to populate regions --------------
  
  DEBUG<-0
  
  IR_df<-data.frame()
  
  
  for(i in 1:dim(ALL_genes)[1])
  {
    ALL_genes_sel<-ALL_genes[i,]
    
    
    Symbol_sel<-ALL_genes_sel$Symbol
    Display_genes_sel<-unique(unlist(strsplit(ALL_genes_sel$Display_genes, split=";")))
    
    cat("-------------------------------------->\t")
    cat(sprintf(as.character(paste(Symbol_sel,Display_genes_sel,collapse=" "))))
    cat("\n")
    
    if(Symbol_sel == 'HBA1')
    {
      DEBUG<-1
      
    }else{
      
      DEBUG<-0
    }
    
    if(DEBUG == 1)
    {
      cat("ALL_genes_sel_0\n")
      cat(str(ALL_genes_sel))
      cat("\n")
    }
    
    ensembl_gtf_gene_sel<-ensembl_gtf_gene[which(ensembl_gtf_gene$gene_name%in%Display_genes_sel |
                                                   ensembl_gtf_gene$gene_id%in%Display_genes_sel),]
    
    if(DEBUG == 1)
    {
      cat("ensembl_gtf_gene_sel_0\n")
      cat(str(ensembl_gtf_gene_sel))
      cat("\n")
    }
    
    Master_peak_file_sel<-Master_peak_file[which(Master_peak_file$Symbol_string%in%Display_genes_sel),]
    
    if(DEBUG == 1)
    {
      cat("Master_peak_file_sel_0\n")
      cat(str(Master_peak_file_sel))
      cat("\n")
    }
    
    chr_sel<-unique(c(paste('chr',ensembl_gtf_gene_sel$seqid,sep=''),Master_peak_file_sel$chr))
    
    if(DEBUG == 1)
    {
      cat("chr_sel_0\n")
      cat(str(chr_sel))
      cat("\n")
    }
    
    # Master_peak_file_sel$start,
    
    START_vector<-c(ensembl_gtf_gene_sel$start,Master_peak_file_sel$start)
    
    min_START<-min(START_vector)-1000
    
    if(min_START < 0){
      
      min_START<-0
      
    }#min_START
    
    if(DEBUG == 1)
    {
      cat("START_vector_0\n")
      cat(str(START_vector))
      cat("\n")
      cat(str(min_START))
      cat("\n")
    }
    
    # Master_peak_file_sel$end,
    
    END_vector<-c(ensembl_gtf_gene_sel$end,Master_peak_file_sel$end)
    
    max_END<-max(END_vector)+1000
    
    if(DEBUG == 1)
    {
      cat("END_vector_0\n")
      cat(str(END_vector))
      cat("\n")
      cat(str(max_END))
      cat("\n")
    }
    
    ALL_genes_sel$chr<-chr_sel
    ALL_genes_sel$start<-min_START
    ALL_genes_sel$end<-max_END
    ALL_genes_sel$region<-paste(ALL_genes_sel$chr,
                                     ALL_genes_sel$start,
                                     ALL_genes_sel$end,
                                     sep='-')
    
    ALL_genes_sel$region_name<-paste('Region',Symbol_sel,ALL_genes_sel$region,sep='__')
    
    
    IR_df<-rbind(ALL_genes_sel,IR_df)
    
    if(DEBUG == 1)
    {
      cat("IR_df_0\n")
      cat(str(IR_df))
      cat("\n")
      
    }
   
  }# i in 1:dim(ALL_genes)[1]
  
  
  #### SAVE RDS ----

  cat("IR_df_FINAL\n")
  cat(str(IR_df))
  cat("\n")
  cat(str(unique(IR_df$region_name)))
  cat("\n")

  setwd(out)

  saveRDS(IR_df,file='Input_regions.rds')
  write.table(IR_df,file='Input_regions.tsv',sep="\t",quote=F,row.names = F)

 
  
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
    make_option(c("--Selected_genes_classified"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--DE_genes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ANAPC_genes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ery_genes_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--myeloid_genes_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--megak_genes_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--marker_genes_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ensembl_gtf"), type="character", default=NULL, 
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
  
  
  Build_input_regions(opt)
  

  
}


###########################################################################

system.time( main() )