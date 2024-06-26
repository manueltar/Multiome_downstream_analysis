
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


classify_peaks_by_TSS_proximity = function(option_list)
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
  
  #### READ and transform distance_to_TSS ----
  
  distance_to_TSS = opt$distance_to_TSS
  
  cat("distance_to_TSS_\n")
  cat(sprintf(as.character(distance_to_TSS)))
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
  
  ensembl_gtf_gene$region_START[indx.positive_strand]<-ensembl_gtf_gene$start[indx.positive_strand] - distance_to_TSS
  ensembl_gtf_gene$region_END[indx.positive_strand]<-ensembl_gtf_gene$start[indx.positive_strand] + distance_to_TSS
  
  
  indx.negative_strand<-which(ensembl_gtf_gene$strand == '-')
  
  cat("indx.negative_strand_0\n")
  cat(str(indx.negative_strand))
  cat("\n")
  
  ensembl_gtf_gene$region_START[indx.negative_strand]<-ensembl_gtf_gene$end[indx.negative_strand] - distance_to_TSS
  ensembl_gtf_gene$region_END[indx.negative_strand]<-ensembl_gtf_gene$end[indx.negative_strand] + distance_to_TSS

  ensembl_gtf_gene$region_START[which(ensembl_gtf_gene$region_START < distance_to_TSS)]<-1
  ensembl_gtf_gene$region_END[which(ensembl_gtf_gene$region_END < distance_to_TSS)]<-distance_to_TSS
  
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
  
  check_negative<-ensembl_gtf_gene[which(ensembl_gtf_gene$gene_name == 'HSD11B1-AS1'),]
  
  
  cat("check_negative_0\n")
  cat(str(check_negative))
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
  
 
  
 
  #### Read the Peak to genes results file ----
  
  ALL_PoI<-readRDS(file=opt$ALL_PoI)
  
  # colnames(ALL_PoI)[which(colnames(ALL_PoI) == 'peak')]<-'Peak_ID'
  
  
  cat("ALL_PoI_0\n")
  cat(str(ALL_PoI))
  cat("\n")
 
  
  gr_ALL_PoI <- GRanges(
    seqnames = as.character(gsub("chr","",ALL_PoI$chr)),
    ranges=IRanges(
      start=as.numeric(ALL_PoI$start),
      end=as.numeric(ALL_PoI$end),
      name=ALL_PoI$Peak_ID))
  
  cat("gr_ALL_PoI_0\n")
  cat(str(gr_ALL_PoI))
  cat("\n")
  
  
 
  #### Intersect Peak to genes with TSS ----
  
  DEBUG <- 1
  
  m <- findOverlaps(gr_TSS,gr_ALL_PoI,
                    ignore.strand	= TRUE)
  
  if(DEBUG == 1)
  {
    cat("m\n")
    cat(str(m))
    cat("\n")
  }
  
  subjectHits_ALL_PoI<-subjectHits(m)
  
  if(DEBUG == 1)
  {
    cat("subjectHits_ALL_PoI\n")
    cat(str(subjectHits_ALL_PoI))
    cat("\n")
  }
  
  queryHits_TSS<-queryHits(m)
  
  if(DEBUG == 1)
  {
    cat("queryHits_TSS\n")
    cat(str(queryHits_TSS))
    cat("\n")
  }
  
  TSS_df <- data.frame(chr=as.character(seqnames(gr_TSS)),
                                                          start=as.integer(start(gr_TSS)),
                                                          end=as.integer(end(gr_TSS)),
                                                          strand=as.character(strand(gr_TSS)),
                                                          TSS_Symbol=as.character(gr_TSS$name2),
                                                          gene_biotype=as.character(gr_TSS$name3),
                                                          ensembl_gene_id=names(gr_TSS), stringsAsFactors = F)
  
  if(DEBUG == 1)
  {
    cat("TSS_df_0\n")
    cat(str(TSS_df))
    cat("\n")
  }
  
  TSS_df_hits<-TSS_df[queryHits_TSS,]
  
  if(DEBUG == 1)
  {
    cat("TSS_df_hits_0\n")
    cat(str(TSS_df_hits))
    cat("\n")
  }
  
  ALL_PoI_df_TSS <- data.frame(chr=as.character(seqnames(gr_ALL_PoI)),
                                                                        start=as.integer(start(gr_ALL_PoI)),
                                                                        end=as.integer(end(gr_ALL_PoI)),
                                                                        Peak_ID=names(gr_ALL_PoI),
                                                                        stringsAsFactors = F)
  
  
  if(DEBUG == 1)
  {
    cat("ALL_PoI_df_TSS_0\n")
    cat(str(ALL_PoI_df_TSS))
    cat("\n")
  }
  
  ALL_PoI_df_TSS_hits<-ALL_PoI_df_TSS[subjectHits_ALL_PoI,]
  
  if(dim(ALL_PoI_df_TSS_hits)[1] >0)
  {
    if(DEBUG == 1)
    {
      cat("ALL_PoI_df_TSS_hits_0\n")
      cat(str(ALL_PoI_df_TSS_hits))
      cat("\n")
    }
    
    ALL_PoI_df_TSS_hits<-cbind(ALL_PoI_df_TSS_hits,TSS_df_hits)
    
    if(DEBUG == 1)
    {
      cat("ALL_PoI_df_TSS_hits_1\n")
      cat(str(ALL_PoI_df_TSS_hits))
      cat("\n")
    }
    
    
    
    ALL_PoI_df_TSS_hits_subset<-unique(ALL_PoI_df_TSS_hits[,c(1:4,c(
      which(colnames(ALL_PoI_df_TSS_hits) == 'TSS_Symbol'),which(colnames(ALL_PoI_df_TSS_hits) == 'gene_biotype'),
      which(colnames(ALL_PoI_df_TSS_hits) == 'strand'),which(colnames(ALL_PoI_df_TSS_hits) == 'ensembl_gene_id')))])
    
    if(DEBUG == 1)
    {
      cat("ALL_PoI_df_TSS_hits_subset_0\n")
      cat(str(ALL_PoI_df_TSS_hits_subset))
      cat("\n")
    }
    
    #### only keep promoters that match their gene with the gene regulated by the peak ----
    
    
    
    setwd(out)

    saveRDS(file='ALL_PoI_in_TSS.rds', ALL_PoI_df_TSS_hits_subset)
    write.table(ALL_PoI_df_TSS_hits_subset, file="ALL_PoI_in_TSS.tsv",sep="\t",quote=F, row.names = F)
    
  }#dim(ALL_PoI_df_TSS_hits)[1] >0
  
  
  

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
    make_option(c("--distance_to_TSS"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ensembl_gtf"), type="character", default=NULL, 
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
  
  classify_peaks_by_TSS_proximity(opt)
 

  
}


###########################################################################

system.time( main() )