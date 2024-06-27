
.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_downstream/lib/R/library"))
.libPaths()

Sys.setenv(RETICULATE_PYTHON="/home/manuel.tardaguila/conda_envs/multiome_downstream/bin/python")
library(reticulate)
reticulate::use_python("/home/manuel.tardaguila/conda_envs/multiome_downstream/bin/python")
reticulate::use_condaenv("/home/manuel.tardaguila/conda_envs/multiome_downstream")
reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg')
suppressMessages(library("optparse"))
suppressMessages(library(hdf5r))
suppressMessages(library(Seurat)) 
suppressMessages(library(Signac)) 
suppressMessages(library(EnsDb.Hsapiens.v86)) 
suppressMessages(library(dplyr)) 
suppressMessages(library(ggplot2)) 
suppressMessages(library(Matrix)) 
suppressMessages(library(data.table)) 
suppressMessages(library(ggpubr)) 
suppressMessages(library(ggplot2))
suppressMessages(library(chromVAR))
suppressMessages(library(enrichR))
suppressMessages(library(JASPAR2020))
suppressMessages(library(TFBSTools))
suppressMessages(library(motifmatchr))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(pheatmap))
suppressMessages(library(presto))
suppressMessages(library("qlcMatrix"))
suppressMessages(library("cowplot"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("plyr"))
suppressMessages(library("forcats"))
suppressMessages(library('ggeasy'))
suppressMessages(library('dplyr'))
suppressMessages(library("svglite"))
suppressMessages(library("ape"))
suppressMessages(library("ggforce"))
suppressMessages(library("ggrepel"))
suppressMessages(library("ggnewscale"))
suppressMessages(library('ggtext'))
suppressMessages(library("tidyr"))
suppressMessages(library("edgeR"))
suppressMessages(library("apeglm"))
suppressMessages(library("DESeq2"))
suppressMessages(library("tibble"))


opt = NULL

options(warn = 0)


data_wrangling = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform cluster_group ----
  
  cluster_group = opt$cluster_group
  
  cat("cluster_group_\n")
  cat(sprintf(as.character(cluster_group)))
  cat("\n")
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  path_graphs = paste(out,'graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
 
  #### READ Seurat_object ----
  
  Seurat_object<-readRDS(file=opt$Seurat_object)
  

  # cat("Seurat_object_0\n")
  # cat(str(Seurat_object))
  # cat("\n")
  
  #### Define sample_id concept ----
  
  Seurat_object$sample_id<-droplevels(interaction(Seurat_object$time_point, Seurat_object$Assigned_GFPbc, sep="__", lex.order=TRUE))
  
  cat(sprintf(as.character(names(summary(Seurat_object$sample_id)))))
  cat("\n")
  cat(sprintf(as.character(summary(Seurat_object$sample_id))))
  cat("\n")
  
  
  #####extract metadata -----
  
  metadata<-Seurat_object[[]]
  
  # cat("metadata\n")
  # cat(str(metadata))
  # cat("\n")
  
  cat(sprintf(as.character(names(summary(metadata$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(metadata$time_point))))
  cat("\n")
  
  
  cat(sprintf(as.character(names(summary(metadata$diff_groups)))))
  cat("\n")
  cat(sprintf(as.character(summary(metadata$diff_groups))))
  cat("\n")
  
  
  cat(sprintf(as.character(names(summary(as.factor(metadata$Assigned_GFPbc))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(metadata$Assigned_GFPbc)))))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(metadata$sample_id))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(metadata$sample_id)))))
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
  
  #### READ and transform EZH2_signature ----
  
  EZH2_signature = unique(unlist(strsplit(opt$EZH2_signature, split=',')))
  
  cat("EZH2_signature_\n")
  cat(sprintf(as.character(EZH2_signature)))
  cat("\n")
  
  #### READ and transform CUX1 ----
  
  CUX1 = unique(unlist(strsplit(opt$CUX1, split=',')))
  
  cat("CUX1_\n")
  cat(sprintf(as.character(CUX1)))
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
  cat(str(orginal_levels_of_Selected_genes_classified))
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
  
  
  # ############################################################################## LINK PEAK TO GENES #####################################################################################################################################################################
  
  Seurat_object <- LinkPeaks(
    object = Seurat_object,
    peak.assay = "ATAC",
    expression.assay = "SCT",
    genes.use = unique(ALL_genes$Symbol))
  
  
  ##################### Use the result -------------
  
  links = Links(Seurat_object)
  lf = as.data.frame(links)
  colnames(lf)[which(colnames(lf) == 'peak')]<-'Peak_ID'
  colnames(lf)[which(colnames(lf) == 'gene')]<-'Symbol'
  
  cat("lf_0\n")
  cat(str(lf))
  cat("\n")
  
  
  # readRDS(file=paste("Peak_to_selected_genes_assigned",".rds", sep=''))
  
  Links_Peaks_to_assigned_genes<-lf
  
  Links_Peaks_to_assigned_genes$Minus_logpval<-round(-1*log10(Links_Peaks_to_assigned_genes$pvalue),2)
  
  cat("Links_Peaks_to_assigned_genes_0\n")
  cat(str(Links_Peaks_to_assigned_genes))
  cat("\n")
  cat(str(unique(Links_Peaks_to_assigned_genes$gene)))
  cat("\n")
  cat(str(unique(Links_Peaks_to_assigned_genes$Peak_ID)))
  cat("\n")
  
  gr_Links <- GRanges(
    seqnames = as.character(gsub("^chr","",Links_Peaks_to_assigned_genes$seqnames)),
    name2=as.character(Links_Peaks_to_assigned_genes$Peak_ID),
    name3=as.numeric(Links_Peaks_to_assigned_genes$zscore),
    name4=as.numeric(Links_Peaks_to_assigned_genes$Minus_logpval),
    ranges=IRanges(
      start=as.numeric(Links_Peaks_to_assigned_genes$start),
      end=as.numeric(Links_Peaks_to_assigned_genes$end),
      name=paste('Link',as.character(Links_Peaks_to_assigned_genes$seqnames),
                 as.character(Links_Peaks_to_assigned_genes$start),
                 as.character(Links_Peaks_to_assigned_genes$end),                                   
                 sep='_')))
  
  # cat("gr_Links_0\n")
  # cat(str(gr_Links))
  # cat("\n")
  
  ALL_Linked_Peaks<-Links_Peaks_to_assigned_genes$Peak_ID
  
  gr_Linked_Peaks <- GRanges(
    seqnames = as.character(gsub("^chr","",gsub("-.+$","",ALL_Linked_Peaks))),  
    ranges=IRanges(
      start=as.integer(gsub("-.+$","",gsub("^[^-]+-","",ALL_Linked_Peaks))), 
      end=as.integer(gsub("^[^-]+-[^-]+-","",ALL_Linked_Peaks)),
      name=paste("Region", Links_Peaks_to_assigned_genes$Symbol, sep="__")))
  
  # cat("gr_Linked_Peaks_0\n")
  # cat(str(gr_Linked_Peaks))
  # cat("\n")
  
  ######### Append all sets of peaks of interest --------------
  
  gr_Peaks_DEF<-gr_Linked_Peaks
  
  Peaks_DEF_df <- unique(data.frame(chr=as.character(seqnames(gr_Peaks_DEF)),
                                    start=as.integer(start(gr_Peaks_DEF)),
                                    end=as.integer(end(gr_Peaks_DEF)), stringsAsFactors = F))
  
  Peaks_DEF_df$chr<-paste('chr',Peaks_DEF_df$chr,sep='')
  
  Peaks_DEF_df$Peak_ID<-paste(Peaks_DEF_df$chr,Peaks_DEF_df$start,Peaks_DEF_df$end,sep='-')
  
  cat("Peaks_DEF_df_0\n")
  cat(str(Peaks_DEF_df))
  cat("\n")
  
  ########### save and classify the new peaks with the next Rscript -------------
  
  setwd(out)
  
  saveRDS(Peaks_DEF_df, file="ALL_PoI.rds")
  
  saveRDS(Links_Peaks_to_assigned_genes, file ="Linked_peak_to_selected_genes.rds")
  write.table(Links_Peaks_to_assigned_genes, file ="Linked_peak_to_selected_genes.tsv",sep="\t",quote=F, row.names = F)
  
  
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
    make_option(c("--Seurat_object"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Selected_genes_classified"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--DE_genes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ANAPC_genes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--EZH2_signature"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CUX1"), type="character", default=NULL, 
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
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
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
  
  data_wrangling(opt)
 

}


###########################################################################

system.time( main() )
  