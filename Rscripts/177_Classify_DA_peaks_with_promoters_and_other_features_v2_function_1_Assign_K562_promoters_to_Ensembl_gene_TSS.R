
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

classify_promoters = function(option_list)
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
  
  #### READ and transform Promoter_distance_to_TSS ----
  
  Promoter_distance_to_TSS = opt$Promoter_distance_to_TSS
  
  cat("Promoter_distance_to_TSS_\n")
  cat(sprintf(as.character(Promoter_distance_to_TSS)))
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
  
  K562_Regulatory_Build_promoter<-K562_Regulatory_Build[which(K562_Regulatory_Build$type == 'promoter'),]
  
  
  cat("K562_Regulatory_Build_promoter_0\n")
  cat(str(K562_Regulatory_Build_promoter))
  cat("\n")
  cat(sprintf(as.character(names(summary(K562_Regulatory_Build_promoter$type)))))
  cat("\n")
  cat(sprintf(as.character(summary(K562_Regulatory_Build_promoter$type))))
  cat("\n")
  
  gr_K562_Regulatory_Build_promoter <- GRanges(
    seqnames = as.character(K562_Regulatory_Build_promoter$seqid),
    name2=as.character(K562_Regulatory_Build_promoter$feature_type),
    name3=as.character(K562_Regulatory_Build_promoter$activity),
    ranges=IRanges(
      start=as.numeric(K562_Regulatory_Build_promoter$start),
      end=as.numeric(K562_Regulatory_Build_promoter$end),
      name=K562_Regulatory_Build_promoter$regulatory_feature_stable_id))
  
  check_chr7<-K562_Regulatory_Build_promoter[which(K562_Regulatory_Build_promoter$seqid == '7'),]
  
  
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
  
  ensembl_gtf_gene$region_START[indx.positive_strand]<-ensembl_gtf_gene$start[indx.positive_strand] - Promoter_distance_to_TSS
  ensembl_gtf_gene$region_END[indx.positive_strand]<-ensembl_gtf_gene$start[indx.positive_strand] + Promoter_distance_to_TSS
  
  
  ensembl_gtf_gene$region_START[indx.negative_strand]<-ensembl_gtf_gene$end[indx.negative_strand] - Promoter_distance_to_TSS
  ensembl_gtf_gene$region_END[indx.negative_strand]<-ensembl_gtf_gene$end[indx.negative_strand] + Promoter_distance_to_TSS
  
  ensembl_gtf_gene$region_START[which(ensembl_gtf_gene$region_START < Promoter_distance_to_TSS)]<-1
  ensembl_gtf_gene$region_END[which(ensembl_gtf_gene$region_END < Promoter_distance_to_TSS)]<-Promoter_distance_to_TSS
  
  ensembl_gtf_gene$region_START<-as.integer(ensembl_gtf_gene$region_START)
  ensembl_gtf_gene$region_END<-as.integer(ensembl_gtf_gene$region_END)
  
  
  indx.discordance<-which(ensembl_gtf_gene$region_START >= ensembl_gtf_gene$region_END)
  
  check<-ensembl_gtf_gene[indx.discordance,]
  
  cat("ensembl_gtf_gene_0\n")
  cat(str(ensembl_gtf_gene))
  cat("\n")
  cat(sprintf(as.character(names(summary(ensembl_gtf_gene$region_START)))))
  cat("\n")
  cat(sprintf(as.character(summary(ensembl_gtf_gene$region_START))))
  cat("\n")
  cat(sprintf(as.character(names(summary(ensembl_gtf_gene$region_END)))))
  cat("\n")
  cat(sprintf(as.character(summary(ensembl_gtf_gene$region_END))))
  cat("\n")
  
  cat("check_0\n")
  cat(str(check))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$region_START)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$region_START))))
  cat("\n")
  cat(sprintf(as.character(names(summary(check$region_END)))))
  cat("\n")
  cat(sprintf(as.character(summary(check$region_END))))
  cat("\n")
  
  check_tracking_genes<-ensembl_gtf_gene[which(ensembl_gtf_gene$gene_name%in%tracking_genes),]
  
  
  cat("check_tracking_genes_0\n")
  cat(str(check_tracking_genes))
  cat("\n")
  
 
  
  rm(ensembl_gtf)
  
  gr_TSS <- GRanges(
    seqnames = as.character(ensembl_gtf_gene$seqid),
    name2=as.character(ensembl_gtf_gene$gene_name),
    name3=as.character(ensembl_gtf_gene$gene_biotype),
    strand=ensembl_gtf_gene$strand,
    ranges=IRanges(
      start=as.numeric(ensembl_gtf_gene$region_START),
      end=as.numeric(ensembl_gtf_gene$region_END),
      name=ensembl_gtf_gene$gene_id))
  
  
  
  #### Intersect TSS to genes with Regulatory features ----
  
  DEBUG <- 1
  
  m <- findOverlaps(gr_K562_Regulatory_Build_promoter,gr_TSS,
                    ignore.strand = TRUE)
  
  if(DEBUG == 1)
  {
    cat("m\n")
    cat(str(m))
    cat("\n")
  }
  
  subjectHits_TSS<-subjectHits(m)
  
  if(DEBUG == 1)
  {
    cat("subjectHits_TSS\n")
    cat(str(subjectHits_TSS))
    cat("\n")
  }
  
  queryHits_K562_Regulatory_Build_promoter<-queryHits(m)
  
  if(DEBUG == 1)
  {
    cat("queryHits_K562_Regulatory_Build_promoter\n")
    cat(str(queryHits_K562_Regulatory_Build_promoter))
    cat("\n")
  }
  
  K562_Regulatory_Build_promoter_df <- data.frame(chr=as.character(seqnames(gr_K562_Regulatory_Build_promoter)),
                                         start=as.integer(start(gr_K562_Regulatory_Build_promoter)),
                                         end=as.integer(end(gr_K562_Regulatory_Build_promoter)),
                                         feature_type=as.character(gr_K562_Regulatory_Build_promoter$name2),
                                         activity=as.character(gr_K562_Regulatory_Build_promoter$name3),
                                         regulatory_feature_stable_id=names(gr_K562_Regulatory_Build_promoter), stringsAsFactors = F)
  
  if(DEBUG == 1)
  {
    cat("K562_Regulatory_Build_promoter_df_0\n")
    cat(str(K562_Regulatory_Build_promoter_df))
    cat("\n")
  }
  
  K562_Regulatory_Build_promoter_df_hits<-K562_Regulatory_Build_promoter_df[queryHits_K562_Regulatory_Build_promoter,]
  
  if(DEBUG == 1)
  {
    cat("K562_Regulatory_Build_promoter_df_hits_0\n")
    cat(str(K562_Regulatory_Build_promoter_df_hits))
    cat("\n")
  }
  
  TSS_df_K562_Regulatory_Build_promoter <- data.frame(chr=as.character(seqnames(gr_TSS)),
                                                       start=as.integer(start(gr_TSS)),
                                                       end=as.integer(end(gr_TSS)),
                                                       Symbol=as.character(gr_TSS$name2),
                                                       gene_biotype=as.character(gr_TSS$name3),
                                                       strand=strand(gr_TSS),
                                                       ensembl_gene_id=names(gr_TSS),
                                                       stringsAsFactors = F)
  
  
  if(DEBUG == 1)
  {
    cat("TSS_df_K562_Regulatory_Build_promoter_0\n")
    cat(str(TSS_df_K562_Regulatory_Build_promoter))
    cat("\n")
  }
  
  TSS_df_K562_Regulatory_Build_promoter_hits<-TSS_df_K562_Regulatory_Build_promoter[subjectHits_TSS,]
  
  if(dim(TSS_df_K562_Regulatory_Build_promoter_hits)[1] >0)
  {
    if(DEBUG == 1)
    {
      cat("TSS_df_K562_Regulatory_Build_promoter_hits_0\n")
      cat(str(TSS_df_K562_Regulatory_Build_promoter_hits))
      cat("\n")
    }
    
    TSS_df_K562_Regulatory_Build_promoter_hits<-cbind(TSS_df_K562_Regulatory_Build_promoter_hits,K562_Regulatory_Build_promoter_df_hits)
    
    if(DEBUG == 1)
    {
      cat("TSS_df_K562_Regulatory_Build_promoter_hits_1\n")
      cat(str(TSS_df_K562_Regulatory_Build_promoter_hits))
      cat("\n")
    }
    
   
    
    TSS_df_K562_Regulatory_Build_promoter_hits_subset<-unique(TSS_df_K562_Regulatory_Build_promoter_hits[,c(1:7,c(which(colnames(TSS_df_K562_Regulatory_Build_promoter_hits) == 'feature_type'),
                                                                                                                    which(colnames(TSS_df_K562_Regulatory_Build_promoter_hits) == 'activity'),
                                                                                                                    which(colnames(TSS_df_K562_Regulatory_Build_promoter_hits) == 'regulatory_feature_stable_id')))])
    
    if(DEBUG == 1)
    {
      cat("TSS_df_K562_Regulatory_Build_promoter_hits_subset_0\n")
      cat(str(TSS_df_K562_Regulatory_Build_promoter_hits_subset))
      cat("\n")
    }
    
    TSS_df_K562_Regulatory_Build_promoter_hits_subset$variable<-paste(TSS_df_K562_Regulatory_Build_promoter_hits_subset$feature_type,TSS_df_K562_Regulatory_Build_promoter_hits_subset$activity,sep='|')
    TSS_df_K562_Regulatory_Build_promoter_hits_subset$value<-paste(TSS_df_K562_Regulatory_Build_promoter_hits_subset$regulatory_feature_stable_id,
                                                                    TSS_df_K562_Regulatory_Build_promoter_hits_subset$feature_type,
                                                                    TSS_df_K562_Regulatory_Build_promoter_hits_subset$activity,sep='|')
    
    TSS_df_K562_Regulatory_Build_promoter_hits_subset<-TSS_df_K562_Regulatory_Build_promoter_hits_subset[,-c(which(colnames(TSS_df_K562_Regulatory_Build_promoter_hits_subset) == 'regulatory_feature_stable_id'),
                                                                                                               which(colnames(TSS_df_K562_Regulatory_Build_promoter_hits_subset) == 'feature_type'),
                                                                                                               which(colnames(TSS_df_K562_Regulatory_Build_promoter_hits_subset) == 'activity'))]
    
    if(DEBUG == 1)
    {
      cat("TSS_df_K562_Regulatory_Build_promoter_hits_subset_1\n")
      cat(str(TSS_df_K562_Regulatory_Build_promoter_hits_subset))
      cat("\n")
    }
    
   
    
    DEF.dt<-data.table(TSS_df_K562_Regulatory_Build_promoter_hits_subset, key=c("chr","start","end","Symbol","gene_biotype","strand","ensembl_gene_id"))
    
    DEF_collapsed<-as.data.frame(DEF.dt[,.(variable_string=paste(unique(sort(variable)), collapse=';'),
                                           value_string=paste(value, collapse=';')),by=key(DEF.dt)], stringsAsFactors=F)
    
    
    
    
    if(DEBUG == 1)
    {
      cat("DEF_collapsed_0\n")
      cat(str(DEF_collapsed))
      cat("\n")
      cat(sprintf(paste(as.character(names(summary(as.factor(DEF_collapsed$variable_string)))), collapse="BOOM")))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(DEF_collapsed$variable_string)))))
      cat("\n")
    }
    
    ### classification of promoter activity
    
    DEF_collapsed$CLASS<-NA
    
    Promoter_ACTIVE_CLASS<-c('Promoter|ACTIVE','Promoter|ACTIVE;Promoter|INACTIVE','Promoter|ACTIVE;Promoter|POISED','Promoter|ACTIVE;Promoter|REPRESSED',
                             'Promoter|ACTIVE;Promoter|INACTIVE;Promoter|POISED','Promoter|ACTIVE;Promoter|INACTIVE;Promoter|REPRESSED',
                             'Promoter|ACTIVE;Promoter|INACTIVE;Promoter|POISED;Promoter|REPRESSED')
    Promoter_INACTIVE_CLASS<-c('Promoter|INACTIVE','Promoter|INACTIVE;Promoter|POISED','Promoter|INACTIVE;Promoter|POISED;Promoter|REPRESSED','Promoter|INACTIVE;Promoter|REPRESSED')
    Promoter_POISED_CLASS<-c('Promoter|POISED','Promoter|POISED;Promoter|REPRESSED')
    Promoter_REPRESSED_CLASS<-c('Promoter|REPRESSED')
    
    DEF_collapsed$CLASS[which(DEF_collapsed$variable_string%in%Promoter_ACTIVE_CLASS)]<-'Promoter_ACTIVE'
    DEF_collapsed$CLASS[which(DEF_collapsed$variable_string%in%Promoter_INACTIVE_CLASS)]<-'Promoter_INACTIVE'
    DEF_collapsed$CLASS[which(DEF_collapsed$variable_string%in%Promoter_POISED_CLASS)]<-'Promoter_POISED'
    DEF_collapsed$CLASS[which(DEF_collapsed$variable_string%in%Promoter_REPRESSED_CLASS)]<-'Promoter_REPRESSED'
    
    DEF_collapsed$feature<-gsub("_.+$","",DEF_collapsed$CLASS)
    DEF_collapsed$activity<-gsub("^[^_]+_","",DEF_collapsed$CLASS)
    
    if(DEBUG == 1)
    {
      cat("DEF_collapsed_1\n")
      cat(str(DEF_collapsed))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(DEF_collapsed$CLASS))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(DEF_collapsed$CLASS)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(DEF_collapsed$feature))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(DEF_collapsed$feature)))))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(DEF_collapsed$activity))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(DEF_collapsed$activity)))))
      cat("\n")
    }
    
    check.NA<-DEF_collapsed[is.na(DEF_collapsed$CLASS),]
    
    if(DEBUG == 1)
    {
      cat("check.NA_0\n")
      cat(str(check.NA))
      cat("\n")
      cat(sprintf(paste(as.character(names(summary(as.factor(check.NA$variable_string)))), collapse="BOOM")))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(check.NA$variable_string)))))
      cat("\n")
    }
    
    check.dt<-data.table(TSS_df_K562_Regulatory_Build_promoter_hits_subset, key=c("variable","value"))
    
    check_collapsed<-as.data.frame(check.dt[,.(ensembl_gene_id_string=paste(unique(sort(ensembl_gene_id)), collapse=';'),
                                               Symbol_string=paste(unique(sort(Symbol)), collapse=';')),by=key(check.dt)], stringsAsFactors=F)
    
    
    if(DEBUG == 1)
    {
      cat("check_collapsed_0\n")
      cat(str(check_collapsed))
      cat("\n")
      cat(sprintf(paste(as.character(names(summary(as.factor(check_collapsed$Symbol_string)))), collapse="BOOM")))
      cat("\n")
      cat(sprintf(paste(as.character(summary(as.factor(check_collapsed$Symbol_string)))), collapse="BOOM"))
      cat("\n")
    }
    
    # cat("check_tracking_genes_1\n")
    # cat(str(check_tracking_genes))
    # cat("\n")
    
    #### SAVE RDS ----
    
    reference_output_dir = paste(out,'reference_files','/',sep='')
    
    if (file.exists(reference_output_dir)){
      
      
    }else{
      
      dir.create(file.path(reference_output_dir))
      
    }#reference_output_dir
    
    setwd(reference_output_dir)
    
    saveRDS(DEF_collapsed,file='Ensembl_promoters_K562_linked_to_genes.rds')
    write.table(DEF_collapsed,file='Ensembl_promoters_K562_linked_to_genes.tsv',sep="\t",quote=F,row.names = F)
    
    
    
  }#dim(TSS_df_K562_Regulatory_Build_promoter_hits)[1] >0
  
  
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
    make_option(c("--Promoter_distance_to_TSS"), type="numeric", default=NULL, 
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
  
  classify_promoters(opt)
 

  
}


###########################################################################

system.time( main() )