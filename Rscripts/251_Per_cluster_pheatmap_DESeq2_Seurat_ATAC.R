
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
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Mm.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("pheatmap", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("edgeR", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))



opt = NULL

options(warn = 1)

heatmap_WT = function(option_list)
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
  
  path_graphs = paste(out,'ATAC_heatmaps','/', sep='')
  
  if (file.exists(path_graphs)){
    
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  
  #### READ and transform CUX1 ----
  
  CUX1 = unlist(strsplit(opt$CUX1, split=','))
  
  cat("CUX1_\n")
  cat(sprintf(as.character(CUX1)))
  cat("\n")
  
  #### READ and transform ery_genes_array ----
  
  ery_genes_array = unlist(strsplit(opt$ery_genes_array, split=','))
  
  cat("ery_genes_array_\n")
  cat(sprintf(as.character(ery_genes_array)))
  cat("\n")
  
  #### READ and transform ery_genes_array ----
  
  ery_genes_array = unlist(strsplit(opt$ery_genes_array, split=','))
  
  cat("ery_genes_array_\n")
  cat(sprintf(as.character(ery_genes_array)))
  cat("\n")
  
  #### READ Master_peak_file_with_SNP_numbered_peaks ----
  
  Master_peak_file_with_SNP_numbered_peaks<-readRDS(file=opt$Master_peak_file_with_SNP_numbered_peaks)
  
  
  cat("Master_peak_file_with_SNP_numbered_peaks\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks))
  cat("\n")
  
  Master_peak_file_with_SNP_numbered_peaks_long<-unique(as.data.frame(cSplit(Master_peak_file_with_SNP_numbered_peaks,sep = ';', direction = "long",
                                                                             splitCols = "Symbol_string"),stringsAsFactors=F))
  
  
  colnames(Master_peak_file_with_SNP_numbered_peaks_long)[which(colnames(Master_peak_file_with_SNP_numbered_peaks_long) == 'Symbol_string')]<-'Symbol'
  
  cat("Master_peak_file_with_SNP_numbered_peaks_long_0\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks_long))
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
  
  new_genes_to_add_string<-unique(c(CUX1,ery_genes_array,megak_genes_array,marker_genes_array,myeloid_genes_array,ANAPC_genes))
  
  cat("new_genes_to_add_string_0\n")
  cat(str(new_genes_to_add_string))
  cat("\n")
  
  
  
  new_genes_to_add <- data.frame(matrix(vector(), length(new_genes_to_add_string), dim(Selected_genes_classified_subset)[2],
                                        dimnames=list(c(),
                                                      colnames(Selected_genes_classified_subset))),stringsAsFactors=F)
  
  
  
  new_genes_to_add$Symbol<-new_genes_to_add_string
  
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%ery_genes_array)]<-'Erythrocyte'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%megak_genes_array)]<-'Megakaryocyte'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%myeloid_genes_array)]<-'AML_myeloid'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%marker_genes_array)]<-'Marker_genes'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%ANAPC_genes)]<-'ANAPC_genes'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%CUX1)]<-'CUX1'
  
  
  
  
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
                                          levels=unique(c(levels_selected_genes,'CUX1','Marker_genes','CUX1_target_genes',
                                                          'MEP_progenitor','Megakaryocyte','Platelet','Platelet_volume',
                                                          'Erythrocyte','Mitosis','ANAPC_genes',
                                                          'PU1_target_genes','RUNX1_target_genes','GATA_target_genes',
                                                          'AML_myeloid','REST_of_genes')),
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
  
  #### Merge ALL_genes with Master_peak_file_with_SNP_numbered_peaks_long ----
  
  Master_peak_file_with_SNP_numbered_peaks_long<-merge(Master_peak_file_with_SNP_numbered_peaks_long,
                                                       ALL_genes,
                                                       by='Symbol')
  
  cat("Master_peak_file_with_SNP_numbered_peaks_long_1\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks_long))
  cat("\n")
  cat(str(unique(Master_peak_file_with_SNP_numbered_peaks_long$Peak_ID)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Master_peak_file_with_SNP_numbered_peaks_long$GENE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Master_peak_file_with_SNP_numbered_peaks_long$GENE_CLASS))))
  cat("\n")
  
  #### order the genes by the WT_cluster of genes ---
  
  
  Master_peak_file_with_SNP_numbered_peaks_long.dt<-data.table(Master_peak_file_with_SNP_numbered_peaks_long, key=c('Peak_CLASS','chr','start','end',
                                                                                                                    'Peak_ID','value_string','feature','activity','rs','VAR','Peak_Number',
                                                                                                                    'Total_peaks'))
  
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed<-as.data.frame(Master_peak_file_with_SNP_numbered_peaks_long.dt[,
                                                                                                                     .(Symbol_string = paste(unique(sort(Symbol)), collapse=';'),
                                                                                                                       GENE_CLASS_string = paste(unique(sort(GENE_CLASS)), collapse=';')
                                                                                                                     ),
                                                                                                                     by=key(Master_peak_file_with_SNP_numbered_peaks_long.dt)])
  
  
  
  
  
  
  
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('CUX1_target_genes;AML_myeloid'))]<-'CUX1_target_genes'
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('MEP_progenitor;Platelet',
                                                                                                                                                       'MEP_progenitor;AML_myeloid',
                                                                                                                                                       'MEP_progenitor;GATA_target_genes',
                                                                                                                                                       'MEP_progenitor;Megakaryocyte'))]<-'MEP_progenitor'
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('Marker_genes;GATA_target_genes'))]<-'Marker_genes'
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('Platelet;AML_myeloid'))]<-'Platelet'
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('GATA_target_genes;AML_myeloid'))]<-'GATA_target_genes'
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('Mitosis;GATA_target_genes'))]<-'Mitosis'
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('Megakaryocyte;AML_myeloid'))]<-'Megakaryocyte'
  
  
  
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string<-factor(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string,
                                                                               levels=c('CUX1',
                                                                                        'CUX1_target_genes',
                                                                                        'Erythrocyte',
                                                                                        'GATA_target_genes',
                                                                                        'Marker_genes',
                                                                                        'Megakaryocyte',
                                                                                        'MEP_progenitor',
                                                                                        'Mitosis',
                                                                                        'Platelet',
                                                                                        'RUNX1_target_genes',
                                                                                        'AML_myeloid',
                                                                                        'ANAPC_genes'),
                                                                               ordered=T)
  
  
  
  cat("Master_peak_file_with_SNP_numbered_peaks_collapsed_0\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks_collapsed))
  cat("\n")
  cat(str(unique(Master_peak_file_with_SNP_numbered_peaks_collapsed$Peak_ID)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string)))))
  cat("\n")
  cat(sprintf(as.character(summary(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string))))
  cat("\n")
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed<-Master_peak_file_with_SNP_numbered_peaks_collapsed[order(Master_peak_file_with_SNP_numbered_peaks_collapsed$Peak_ID),]
  
  cat("Master_peak_file_with_SNP_numbered_peaks_collapsed_1\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks_collapsed))
  cat("\n")
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$Peak_ID<-as.character(Master_peak_file_with_SNP_numbered_peaks_collapsed$Peak_ID)
  
  
  Master_peak_file_FEATURE<-Master_peak_file_with_SNP_numbered_peaks_collapsed[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$Peak_CLASS == c('TSS_Peak','Linked_Peak')),]
  
  cat("Master_peak_file_FEATURE_0\n")
  cat(str(Master_peak_file_FEATURE))
  cat("\n")
  cat(str(unique(Master_peak_file_FEATURE$Peak_ID)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Master_peak_file_FEATURE$GENE_CLASS_string)))))
  cat("\n")
  cat(sprintf(as.character(summary(Master_peak_file_FEATURE$GENE_CLASS_string))))
  cat("\n")
  
  #### READ and transform metadata ----
  
  metadata<-readRDS(file=opt$metadata)
  
  cat("metadata_0\n")
  cat(str(metadata))
  cat("\n")
  
  metadata$clone_line<-as.character(metadata$clone_line)
  
  metadata$clone_line<-gsub("^chrGFP_","",metadata$clone_line)
  
  metadata$clone_line<-factor(metadata$clone_line,
                              levels=c("WTA","WTB","WTC","HET","KI_13","KI_27","KI_29","Del_16bp","Del_233","Del_235","Del_287"),
                              ordered=T)
  
  cat("metadata_1\n")
  cat(str(metadata))
  cat("\n")
  
  metadata_WT<-droplevels(metadata[which(metadata$Genotype == 'G.G'),])
  
  cat("metadata_WT_0\n")
  cat(str(metadata_WT))
  cat("\n")
  
  metadata_WT<-metadata_WT[order(metadata_WT$seurat_cluster),]
  
  cat("metadata_WT_ordered\n")
  cat(str(metadata_WT))
  cat("\n")
  
  column_order<-metadata_WT$columns_for_matrix
  
  cat("column_order_0\n")
  cat(sprintf(as.character(column_order)))
  cat("\n")
  
  metadata_WT_subset<-droplevels(unique(metadata_WT[,c(which(colnames(metadata_WT) == 'clone_line'),
                                                       which(colnames(metadata_WT) == 'Genotype'))]))
  
  cat("metadata_WT_subset_0\n")
  cat(str(metadata_WT_subset))
  cat("\n")
  
  
  
  
  
  #### READ and transform ATAC_matrix ----
  
  df_for_matrix<-as.data.frame(fread(file=opt$ATAC_matrix, sep="\t", header=T), stringsAsFactors=F)
  
  cat("df_for_matrix_0\n")
  cat(str(df_for_matrix))
  cat("\n")
  
  df_for_matrix<-merge(df_for_matrix,
                       Master_peak_file_FEATURE,
                       by='Peak_ID')
  
  cat("df_for_matrix_1\n")
  cat(str(df_for_matrix))
  cat("\n")
  
  df_for_matrix_subset<-df_for_matrix[!is.na(df_for_matrix$Peak_Number),]
  
  cat("df_for_matrix_subset_0\n")
  cat(str(df_for_matrix_subset))
  cat("\n")
  
  indx.int<-c(which(colnames(df_for_matrix_subset) == 'Peak_Number'),which(colnames(df_for_matrix_subset)%in%column_order))
  
  
  cat("indx.int\n")
  cat(str(indx.int))
  cat("\n")
  
  df_for_matrix_subset_subset<-unique(df_for_matrix_subset[,indx.int])
  
  cat("df_for_matrix_subset_subset_0\n")
  cat(str(df_for_matrix_subset_subset))
  cat("\n")
  
  indx.reorder<-c(which(colnames(df_for_matrix_subset_subset) == 'Peak_Number'),which(colnames(df_for_matrix_subset_subset)%in%column_order))
  
  
  cat("indx.reorder\n")
  cat(str(indx.reorder))
  cat("\n")
  
  
  # df_for_matrix_subset_reordered[indx.reorder]
  # 
  
  
  df_for_matrix_subset_subset.dt<-setDT(df_for_matrix_subset_subset)
  
  setcolorder(df_for_matrix_subset_subset.dt, c('Peak_Number',column_order))
  
  df_for_matrix_subset_subset_reordered<-as.data.frame(df_for_matrix_subset_subset.dt, stringsAsFactors=F)
  
  cat("df_for_matrix_subset_subset_reordered_0\n")
  cat(str(df_for_matrix_subset_subset_reordered))
  cat("\n")
  
  
  ATAC_matrix<-as.matrix(df_for_matrix_subset_subset_reordered[,-c(which(colnames(df_for_matrix_subset_subset_reordered) == 'Peak_Number'))])
  row.names(ATAC_matrix)<-df_for_matrix_subset_subset_reordered$Peak_Number
  
  cat("ATAC_matrix_0\n")
  cat(str(ATAC_matrix))
  cat("\n")
  
  
  #### Prepare the counts ----
  
  
  lcpm <- cpm(ATAC_matrix, log=TRUE)
  
  cat("lcpm_0\n")
  cat(str(lcpm))
  cat("\n")
  
  
  
  #### subset for genes of interest & columns ----
  
  
  lcpm_subset<-lcpm[which(row.names(lcpm)%in%Master_peak_file_FEATURE$Peak_Number),]
  
  
  cat("lcpm_subset_0\n")
  cat(str(lcpm_subset))
  cat("\n")
  
  cluster<-metadata_WT$seurat_cluster
  
  cat("cluster_0\n")
  cat(str(cluster))
  cat("\n")
  
  
  annot <- data.frame(cluster=paste0("cluster ", cluster))
  
  cat("annot_0\n")
  cat(str(annot))
  cat("\n")
  
  rownames(annot) <- colnames(lcpm_subset)
  
  cat("annot_1\n")
  cat(str(annot))
  cat("\n")
  
  ann_colors <- list(cluster)
  
  cat("ann_colors_0\n")
  cat(str(ann_colors))
  cat("\n")
  
  #### Annotation of rows of heatmap ----
  
  Master_peak_file_FEATURE_subset<-Master_peak_file_FEATURE[!is.na(Master_peak_file_FEATURE$Peak_Number),]
  
  cat("Master_peak_file_FEATURE_subset_0\n")
  cat(str(Master_peak_file_FEATURE_subset))
  cat("\n")
  
  annotation_row = data.frame(GeneClass = Master_peak_file_FEATURE_subset$GENE_CLASS_string)
  
  # ,
  # Symbol = Master_peak_file_FEATURE$Symbol_string
  
  rownames(annotation_row) = Master_peak_file_FEATURE_subset$Peak_Number
  
  
  cat("annotation_row_0\n")
  cat(str(annotation_row))
  cat("\n")
  cat(sprintf(as.character(names(summary(annotation_row$GeneClass)))))
  cat("\n")
  cat(sprintf(as.character(summary(annotation_row$GeneClass))))
  cat("\n")
  
  
  
  
  
  #### Annotation of columns ----
  
  annotation_col<- data.frame(matrix(vector(), length(colnames(lcpm_subset)), 3,
                                     dimnames=list(c(),
                                                   c("time_point","seurat_cluster","clone_line"))),stringsAsFactors=F)
  
  row.names(annotation_col)<-colnames(lcpm_subset)
  
  
  cat("annotation_col_0\n")
  cat(str(annotation_col))
  cat("\n")
  cat(sprintf(as.character(row.names(annotation_col))))
  cat("\n")
  
  
  annotation_col$time_point<-gsub("__.+$","",row.names(annotation_col))
  
  cat("time_point\n")
  cat(sprintf(as.character(names(summary(as.factor(annotation_col$time_point))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(annotation_col$time_point)))))
  cat("\n")
  
  annotation_col$seurat_cluster<-gsub("^[^__]+__","",row.names(annotation_col))
  
  annotation_col$seurat_cluster<-gsub("__.+$","",annotation_col$seurat_cluster)
  
  cat("seurat_cluster\n")
  cat(sprintf(as.character(names(summary(as.factor(annotation_col$seurat_cluster))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(annotation_col$seurat_cluster)))))
  cat("\n")
  
  
  annotation_col$clone_line<-gsub(paste(unique(annotation_col$time_point), collapse='|'),"",row.names(annotation_col))
  annotation_col$clone_line<-gsub(paste(unique(annotation_col$seurat_cluster), collapse='|'),"",annotation_col$clone_line)
  
  annotation_col$clone_line<-gsub("^[^__]+__[^__]+__","",row.names(annotation_col))
  annotation_col$clone_line<-gsub("__.+$","",annotation_col$clone_line)
  
  annotation_col$clone_line<-gsub("^chrGFP_","",annotation_col$clone_line)
  
  
  
  
  cat(sprintf(as.character(names(summary(as.factor(annotation_col$clone_line))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(annotation_col$clone_line)))))
  cat("\n")
  
  annotation_col$Genotype<-NA
  
  annotation_col$Genotype[which(annotation_col$clone_line%in%c('WTA','WTB','WTC'))]<-'G.G'
  annotation_col$Genotype[which(annotation_col$clone_line%in%c('HET'))]<-'A.G'
  annotation_col$Genotype[which(annotation_col$clone_line%in%c('KI_13','KI_27','KI_29'))]<-'A.A'
  annotation_col$Genotype[which(annotation_col$clone_line%in%c('Del_16bp'))]<-'Del16'
  annotation_col$Genotype[which(annotation_col$clone_line%in%c('Del_233','Del_235','Del_287'))]<-'Del80'
  
  cat(sprintf(as.character(names(summary(as.factor(annotation_col$Genotype))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(annotation_col$Genotype)))))
  cat("\n")
  
  
  
  
  annotation_col$time_point<-factor(annotation_col$time_point,
                                    levels=levels(metadata_WT$time_point),
                                    ordered=T)
  
  annotation_col$Genotype<-factor(annotation_col$Genotype,
                                  levels=levels(metadata_WT$Genotype),
                                  ordered=T)
  
  annotation_col$seurat_cluster<-factor(annotation_col$seurat_cluster,
                                        levels=levels(metadata_WT$seurat_cluster),
                                        ordered=T)
  
  annotation_col$clone_line<-factor(annotation_col$clone_line,
                                    levels=levels(metadata_WT$clone_line),
                                    ordered=T)
  
  cat("annotation_col_2\n")
  cat(str(annotation_col))
  cat("\n")
  cat(sprintf(as.character(names(summary(annotation_col$clone_line)))))
  cat("\n")
  cat(sprintf(as.character(summary(annotation_col$clone_line))))
  cat("\n")
  cat(sprintf(as.character(names(summary(annotation_col$seurat_cluster)))))
  cat("\n")
  cat(sprintf(as.character(summary(annotation_col$seurat_cluster))))
  cat("\n")
  cat(sprintf(as.character(names(summary(annotation_col$Genotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(annotation_col$Genotype))))
  cat("\n")
  cat(sprintf(as.character(names(summary(annotation_col$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(annotation_col$time_point))))
  cat("\n")
  
  #### Annotation colors ----
  
  cat("annotation_row_0\n")
  cat(str(annotation_row))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(annotation_row$color))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(annotation_row$color)))))
  cat("\n")
  
  
  ann_colors = list(
    clone_line = c(WTA = brewer.pal(8, "Oranges")[c(3)],
                   WTB = brewer.pal(8, "Oranges")[c(5)],
                   WTC = brewer.pal(8, "Oranges")[c(7)]),
    time_point = c(h24 = brewer.pal(9, "Greys")[c(3)],
                   h48 = brewer.pal(9, "Greys")[c(5)],
                   h72 = brewer.pal(9, "Greys")[c(7)],
                   h96 = brewer.pal(9, "Greys")[c(9)]),
    Genotype = c('G.G' = 'black',
                 'A.G' = brewer.pal(8, "Paired")[c(3)],
                 'A.A' =  brewer.pal(8, "Paired")[c(4)],
                 Del16 = brewer.pal(8, "Paired")[c(5)],
                 Del80 = brewer.pal(8, "Paired")[c(6)]),
    seurat_cluster = c('1' = brewer.pal(9, "YlOrRd")[c(7)],
                       '2' = brewer.pal(9, "Blues")[c(5)],
                       '3' = brewer.pal(9, "RdPu")[c(6)],
                       '4' = brewer.pal(9, "Blues")[c(4)],
                       '5' = brewer.pal(9, "Greens")[c(6)],
                       '6' = brewer.pal(9, "YlOrRd")[c(6)],
                       '7' = brewer.pal(9, "Blues")[c(3)],
                       '8' = brewer.pal(9, "RdPu")[c(5)],
                       '9' = brewer.pal(9, "Greens")[c(5)],
                       '10' = brewer.pal(9, "Greens")[c(4)],
                       '11' = brewer.pal(9, "Blues")[c(2)],
                       '12' = brewer.pal(9, "RdPu")[c(4)],
                       '13' = brewer.pal(9, "RdPu")[c(3)]),
    GeneClass = c('CUX1' = 'black',
                  'Marker_genes' = 'firebrick1',
                  'MEP_progenitor' =  brewer.pal(8, "Dark2")[2],
                  'Megakaryocyte' = brewer.pal(8, "Dark2")[3],
                  'Platelet' = brewer.pal(8, "Dark2")[4],
                  'Platelet_volume'= brewer.pal(8, "Dark2")[8],
                  'Erythrocyte' = brewer.pal(8, "Dark2")[6],
                  'Mitosis' = brewer.pal(8, "BrBG")[1],
                  'ANAPC_genes' = brewer.pal(8, "BrBG")[3],
                  'CUX1_target_genes' = brewer.pal(8, "Dark2")[1],
                  'PU1_target_genes' = brewer.pal(8, "Accent")[3],
                  'RUNX1_target_genes' = brewer.pal(8, "Accent")[2],
                  'GATA_target_genes' = brewer.pal(8, "Accent")[1],
                  'AML_myeloid' ='gray'))
  
  
  
  cat("ann_colors_0\n")
  cat(str(ann_colors))
  cat("\n")
  
  cat("annotation_col_REMEMBER\n")
  cat(str(annotation_col))
  cat("\n")
  cat(sprintf(as.character(row.names(annotation_col))))
  cat("\n")
  
  #### Actual heatmap ----
  
  heatmap<-pheatmap(lcpm_subset, display_numbers = FALSE, number_format = "%.1e",
                    angle_col = "0",clustering_method="ward.D2",
                    fontsize_row = 1, 
                    fontsize_col = 6,
                    fontsize= 6,
                    breaks=seq(-2,2,length.out=101),
                    color=colorRampPalette(c("blue","white","red"))(100),
                    scale="row",
                    cluster_cols=TRUE,
                    cluster_rows=TRUE,
                    border_color='black',
                    treeheight_row=70, treeheight_col=70, cutree_cols=7,
                    annotation_col = annotation_col, annotation_row = annotation_row,
                    annotation_colors = ann_colors,
                    show_colnames=FALSE)
  
  
  svgname<-paste(paste("Heatmap_ATAC_TSS_and_Linked_Peaks","WT_cells", sep='_'),".svg",sep='')
  
  setwd(path_graphs)
  svglite(svgname, width = 13, height = 13)
  print(heatmap)
  dev.off()
  
  
}

heatmap_ALL = function(option_list)
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
  
  path_graphs = paste(out,'ATAC_heatmaps','/', sep='')
  
  if (file.exists(path_graphs)){
    
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
 
  #### READ and transform CUX1 ----
  
  CUX1 = unlist(strsplit(opt$CUX1, split=','))
  
  cat("CUX1_\n")
  cat(sprintf(as.character(CUX1)))
  cat("\n")
  
  #### READ and transform ery_genes_array ----
  
  ery_genes_array = unlist(strsplit(opt$ery_genes_array, split=','))
  
  cat("ery_genes_array_\n")
  cat(sprintf(as.character(ery_genes_array)))
  cat("\n")
  
  #### READ and transform ery_genes_array ----
  
  ery_genes_array = unlist(strsplit(opt$ery_genes_array, split=','))
  
  cat("ery_genes_array_\n")
  cat(sprintf(as.character(ery_genes_array)))
  cat("\n")
  
  #### READ Master_peak_file_with_SNP_numbered_peaks ----
  
  Master_peak_file_with_SNP_numbered_peaks<-readRDS(file=opt$Master_peak_file_with_SNP_numbered_peaks)
  
  
  cat("Master_peak_file_with_SNP_numbered_peaks\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks))
  cat("\n")
  
  Master_peak_file_with_SNP_numbered_peaks_long<-unique(as.data.frame(cSplit(Master_peak_file_with_SNP_numbered_peaks,sep = ';', direction = "long",
                                                                             splitCols = "Symbol_string"),stringsAsFactors=F))
  
  
  colnames(Master_peak_file_with_SNP_numbered_peaks_long)[which(colnames(Master_peak_file_with_SNP_numbered_peaks_long) == 'Symbol_string')]<-'Symbol'
  
  cat("Master_peak_file_with_SNP_numbered_peaks_long_0\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks_long))
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
  
  new_genes_to_add_string<-unique(c(CUX1,ery_genes_array,megak_genes_array,marker_genes_array,myeloid_genes_array,ANAPC_genes))
  
  cat("new_genes_to_add_string_0\n")
  cat(str(new_genes_to_add_string))
  cat("\n")
  
  
  
  new_genes_to_add <- data.frame(matrix(vector(), length(new_genes_to_add_string), dim(Selected_genes_classified_subset)[2],
                                        dimnames=list(c(),
                                                      colnames(Selected_genes_classified_subset))),stringsAsFactors=F)
  
  
  
  new_genes_to_add$Symbol<-new_genes_to_add_string
  
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%ery_genes_array)]<-'Erythrocyte'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%megak_genes_array)]<-'Megakaryocyte'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%myeloid_genes_array)]<-'AML_myeloid'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%marker_genes_array)]<-'Marker_genes'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%ANAPC_genes)]<-'ANAPC_genes'
  new_genes_to_add$GENE_CLASS[which(new_genes_to_add$Symbol%in%CUX1)]<-'CUX1'
  
  
  
  
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
                                          levels=unique(c(levels_selected_genes,'CUX1','Marker_genes','CUX1_target_genes',
                                                          'MEP_progenitor','Megakaryocyte','Platelet','Platelet_volume',
                                                          'Erythrocyte','Mitosis','ANAPC_genes',
                                                          'PU1_target_genes','RUNX1_target_genes','GATA_target_genes',
                                                          'AML_myeloid','REST_of_genes')),
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
  
  #### Merge ALL_genes with Master_peak_file_with_SNP_numbered_peaks_long ----
  
  Master_peak_file_with_SNP_numbered_peaks_long<-merge(Master_peak_file_with_SNP_numbered_peaks_long,
                                                       ALL_genes,
                                                       by='Symbol')
  
  cat("Master_peak_file_with_SNP_numbered_peaks_long_1\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks_long))
  cat("\n")
  cat(str(unique(Master_peak_file_with_SNP_numbered_peaks_long$Peak_ID)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Master_peak_file_with_SNP_numbered_peaks_long$GENE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Master_peak_file_with_SNP_numbered_peaks_long$GENE_CLASS))))
  cat("\n")
  
  #### order the genes by the WT_cluster of genes ---
  
  
  Master_peak_file_with_SNP_numbered_peaks_long.dt<-data.table(Master_peak_file_with_SNP_numbered_peaks_long, key=c('Peak_CLASS','chr','start','end',
                                                                                                                    'Peak_ID','value_string','feature','activity','rs','VAR','Peak_Number',
                                                                                                                    'Total_peaks'))
  
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed<-as.data.frame(Master_peak_file_with_SNP_numbered_peaks_long.dt[,
                                                                                                                     .(Symbol_string = paste(unique(sort(Symbol)), collapse=';'),
                                                                                                                       GENE_CLASS_string = paste(unique(sort(GENE_CLASS)), collapse=';')
                                                                                                                     ),
                                                                                                                     by=key(Master_peak_file_with_SNP_numbered_peaks_long.dt)])
  
  
  
  
  
  
  
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('CUX1_target_genes;AML_myeloid'))]<-'CUX1_target_genes'
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('MEP_progenitor;Platelet',
                                                                                                                                                       'MEP_progenitor;AML_myeloid',
                                                                                                                                                       'MEP_progenitor;GATA_target_genes',
                                                                                                                                                       'MEP_progenitor;Megakaryocyte'))]<-'MEP_progenitor'
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('Marker_genes;GATA_target_genes'))]<-'Marker_genes'
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('Platelet;AML_myeloid'))]<-'Platelet'
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('GATA_target_genes;AML_myeloid'))]<-'GATA_target_genes'
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('Mitosis;GATA_target_genes'))]<-'Mitosis'
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string%in%c('Megakaryocyte;AML_myeloid'))]<-'Megakaryocyte'
  
  
  
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string<-factor(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string,
                                                                               levels=c('CUX1',
                                                                                        'CUX1_target_genes',
                                                                                        'Erythrocyte',
                                                                                        'GATA_target_genes',
                                                                                        'Marker_genes',
                                                                                        'Megakaryocyte',
                                                                                        'MEP_progenitor',
                                                                                        'Mitosis',
                                                                                        'Platelet',
                                                                                        'RUNX1_target_genes',
                                                                                        'AML_myeloid',
                                                                                        'ANAPC_genes'),
                                                                               ordered=T)
  
  
  
  cat("Master_peak_file_with_SNP_numbered_peaks_collapsed_0\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks_collapsed))
  cat("\n")
  cat(str(unique(Master_peak_file_with_SNP_numbered_peaks_collapsed$Peak_ID)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string)))))
  cat("\n")
  cat(sprintf(as.character(summary(Master_peak_file_with_SNP_numbered_peaks_collapsed$GENE_CLASS_string))))
  cat("\n")
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed<-Master_peak_file_with_SNP_numbered_peaks_collapsed[order(Master_peak_file_with_SNP_numbered_peaks_collapsed$Peak_ID),]
  
  cat("Master_peak_file_with_SNP_numbered_peaks_collapsed_1\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks_collapsed))
  cat("\n")
  
  Master_peak_file_with_SNP_numbered_peaks_collapsed$Peak_ID<-as.character(Master_peak_file_with_SNP_numbered_peaks_collapsed$Peak_ID)
  
  
  Master_peak_file_FEATURE<-Master_peak_file_with_SNP_numbered_peaks_collapsed[which(Master_peak_file_with_SNP_numbered_peaks_collapsed$Peak_CLASS == c('TSS_Peak','Linked_Peak')),]
  
  cat("Master_peak_file_FEATURE_0\n")
  cat(str(Master_peak_file_FEATURE))
  cat("\n")
  cat(str(unique(Master_peak_file_FEATURE$Peak_ID)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Master_peak_file_FEATURE$GENE_CLASS_string)))))
  cat("\n")
  cat(sprintf(as.character(summary(Master_peak_file_FEATURE$GENE_CLASS_string))))
  cat("\n")
  
  #### READ and transform metadata ----
  
  metadata<-readRDS(file=opt$metadata)
  
  cat("metadata_0\n")
  cat(str(metadata))
  cat("\n")
  
  metadata$clone_line<-as.character(metadata$clone_line)
  
  metadata$clone_line<-gsub("^chrGFP_","",metadata$clone_line)
  
  metadata$clone_line<-factor(metadata$clone_line,
                              levels=c("WTA","WTB","WTC","HET","KI_13","KI_27","KI_29","Del_16bp","Del_233","Del_235","Del_287"),
                              ordered=T)
  
  cat("metadata_1\n")
  cat(str(metadata))
  cat("\n")
  
  
  
  metadata<-metadata[order(metadata$seurat_cluster),]
  
  cat("metadata_ordered\n")
  cat(str(metadata))
  cat("\n")
  
  column_order<-metadata$columns_for_matrix
  
  cat("column_order_0\n")
  cat(sprintf(as.character(column_order)))
  cat("\n")
  
  metadata_subset<-droplevels(unique(metadata[,c(which(colnames(metadata) == 'clone_line'),
                                                 which(colnames(metadata) == 'Genotype'))]))
  
  cat("metadata_subset_0\n")
  cat(str(metadata_subset))
  cat("\n")
  
  #### READ and transform ATAC_matrix ----
  
  df_for_matrix<-as.data.frame(fread(file=opt$ATAC_matrix, sep="\t", header=T), stringsAsFactors=F)
  
  cat("df_for_matrix_0\n")
  cat(str(df_for_matrix))
  cat("\n")
  
  df_for_matrix<-merge(df_for_matrix,
                       Master_peak_file_FEATURE,
                       by='Peak_ID')
  
  cat("df_for_matrix_1\n")
  cat(str(df_for_matrix))
  cat("\n")
  
  df_for_matrix_subset<-df_for_matrix[!is.na(df_for_matrix$Peak_Number),]
  
  cat("df_for_matrix_subset_0\n")
  cat(str(df_for_matrix_subset))
  cat("\n")
  
  indx.int<-c(which(colnames(df_for_matrix_subset) == 'Peak_Number'),which(colnames(df_for_matrix_subset)%in%column_order))
  
  
  cat("indx.int\n")
  cat(str(indx.int))
  cat("\n")
  
  df_for_matrix_subset_subset<-unique(df_for_matrix_subset[,indx.int])
  
  cat("df_for_matrix_subset_subset_0\n")
  cat(str(df_for_matrix_subset_subset))
  cat("\n")
  
  indx.reorder<-c(which(colnames(df_for_matrix_subset_subset) == 'Peak_Number'),which(colnames(df_for_matrix_subset_subset)%in%column_order))
  
  
  cat("indx.reorder\n")
  cat(str(indx.reorder))
  cat("\n")
  
  
  # df_for_matrix_subset_reordered[indx.reorder]
  # 
  
  
  df_for_matrix_subset_subset.dt<-setDT(df_for_matrix_subset_subset)
  
  setcolorder(df_for_matrix_subset_subset.dt, c('Peak_Number',column_order))
  
  df_for_matrix_subset_subset_reordered<-as.data.frame(df_for_matrix_subset_subset.dt, stringsAsFactors=F)
  
  cat("df_for_matrix_subset_subset_reordered_0\n")
  cat(str(df_for_matrix_subset_subset_reordered))
  cat("\n")
  
  
  ATAC_matrix<-as.matrix(df_for_matrix_subset_subset_reordered[,-c(which(colnames(df_for_matrix_subset_subset_reordered) == 'Peak_Number'))])
  row.names(ATAC_matrix)<-df_for_matrix_subset_subset_reordered$Peak_Number
  
  cat("ATAC_matrix_0\n")
  cat(str(ATAC_matrix))
  cat("\n")
  
  
  #### Prepare the counts ----
  
  
  lcpm <- cpm(ATAC_matrix, log=TRUE)
  
  cat("lcpm_0\n")
  cat(str(lcpm))
  cat("\n")
  
  
  
  #### subset for genes of interest & columns ----
  
  
  lcpm_subset<-lcpm[which(row.names(lcpm)%in%Master_peak_file_FEATURE$Peak_Number),]
  
  
  cat("lcpm_subset_0\n")
  cat(str(lcpm_subset))
  cat("\n")
  
  cluster<-metadata$seurat_cluster
  
  cat("cluster_0\n")
  cat(str(cluster))
  cat("\n")
  
  
  annot <- data.frame(cluster=paste0("cluster ", cluster))
  
  cat("annot_0\n")
  cat(str(annot))
  cat("\n")
  
  rownames(annot) <- colnames(lcpm_subset)
  
  cat("annot_1\n")
  cat(str(annot))
  cat("\n")
  
  ann_colors <- list(cluster)
  
  cat("ann_colors_0\n")
  cat(str(ann_colors))
  cat("\n")
  
  #### Annotation of rows of heatmap ----
  
  Master_peak_file_FEATURE_subset<-Master_peak_file_FEATURE[!is.na(Master_peak_file_FEATURE$Peak_Number),]
  
  cat("Master_peak_file_FEATURE_subset_0\n")
  cat(str(Master_peak_file_FEATURE_subset))
  cat("\n")
  
  annotation_row = data.frame(GeneClass = Master_peak_file_FEATURE_subset$GENE_CLASS_string)
  
  # ,
  # Symbol = Master_peak_file_FEATURE$Symbol_string
  
  rownames(annotation_row) = Master_peak_file_FEATURE_subset$Peak_Number
  
  
  cat("annotation_row_0\n")
  cat(str(annotation_row))
  cat("\n")
  cat(sprintf(as.character(names(summary(annotation_row$GeneClass)))))
  cat("\n")
  cat(sprintf(as.character(summary(annotation_row$GeneClass))))
  cat("\n")
  
  
  
  
  
  #### Annotation of columns ----
  
  annotation_col<- data.frame(matrix(vector(), length(colnames(lcpm_subset)), 3,
                                     dimnames=list(c(),
                                                   c("time_point","seurat_cluster","clone_line"))),stringsAsFactors=F)
  
  row.names(annotation_col)<-colnames(lcpm_subset)
  
  
  cat("annotation_col_0\n")
  cat(str(annotation_col))
  cat("\n")
  cat(sprintf(as.character(row.names(annotation_col))))
  cat("\n")
  
  
  annotation_col$time_point<-gsub("__.+$","",row.names(annotation_col))
  
  cat("time_point\n")
  cat(sprintf(as.character(names(summary(as.factor(annotation_col$time_point))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(annotation_col$time_point)))))
  cat("\n")
  
  annotation_col$seurat_cluster<-gsub("^[^__]+__","",row.names(annotation_col))
  
  annotation_col$seurat_cluster<-gsub("__.+$","",annotation_col$seurat_cluster)
  
  cat("seurat_cluster\n")
  cat(sprintf(as.character(names(summary(as.factor(annotation_col$seurat_cluster))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(annotation_col$seurat_cluster)))))
  cat("\n")
  
  
  annotation_col$clone_line<-gsub(paste(unique(annotation_col$time_point), collapse='|'),"",row.names(annotation_col))
  annotation_col$clone_line<-gsub(paste(unique(annotation_col$seurat_cluster), collapse='|'),"",annotation_col$clone_line)
  
  annotation_col$clone_line<-gsub("^[^__]+__[^__]+__","",row.names(annotation_col))
  annotation_col$clone_line<-gsub("__.+$","",annotation_col$clone_line)
  
  annotation_col$clone_line<-gsub("^chrGFP_","",annotation_col$clone_line)
  
  
  
  
  cat(sprintf(as.character(names(summary(as.factor(annotation_col$clone_line))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(annotation_col$clone_line)))))
  cat("\n")
  
  annotation_col$Genotype<-NA
  
  annotation_col$Genotype[which(annotation_col$clone_line%in%c('WTA','WTB','WTC'))]<-'G.G'
  annotation_col$Genotype[which(annotation_col$clone_line%in%c('HET'))]<-'A.G'
  annotation_col$Genotype[which(annotation_col$clone_line%in%c('KI_13','KI_27','KI_29'))]<-'A.A'
  annotation_col$Genotype[which(annotation_col$clone_line%in%c('Del_16bp'))]<-'Del16'
  annotation_col$Genotype[which(annotation_col$clone_line%in%c('Del_233','Del_235','Del_287'))]<-'Del80'
  
  cat(sprintf(as.character(names(summary(as.factor(annotation_col$Genotype))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(annotation_col$Genotype)))))
  cat("\n")
  
  
  
  
  annotation_col$time_point<-factor(annotation_col$time_point,
                                    levels=levels(metadata$time_point),
                                    ordered=T)
  
  annotation_col$Genotype<-factor(annotation_col$Genotype,
                                  levels=levels(metadata$Genotype),
                                  ordered=T)
  
  annotation_col$seurat_cluster<-factor(annotation_col$seurat_cluster,
                                       levels=levels(metadata$seurat_cluster),
                                       ordered=T)
  
  annotation_col$clone_line<-factor(annotation_col$clone_line,
                                    levels=levels(metadata$clone_line),
                                    ordered=T)
  
  cat("annotation_col_2\n")
  cat(str(annotation_col))
  cat("\n")
  cat(sprintf(as.character(names(summary(annotation_col$clone_line)))))
  cat("\n")
  cat(sprintf(as.character(summary(annotation_col$clone_line))))
  cat("\n")
  cat(sprintf(as.character(names(summary(annotation_col$seurat_cluster)))))
  cat("\n")
  cat(sprintf(as.character(summary(annotation_col$seurat_cluster))))
  cat("\n")
  cat(sprintf(as.character(names(summary(annotation_col$Genotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(annotation_col$Genotype))))
  cat("\n")
  cat(sprintf(as.character(names(summary(annotation_col$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(annotation_col$time_point))))
  cat("\n")
  
  #### Annotation colors ----
  
  cat("annotation_row_0\n")
  cat(str(annotation_row))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(annotation_row$color))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(annotation_row$color)))))
  cat("\n")
  
  
  ann_colors = list(
    clone_line = c(WTA = brewer.pal(8, "Oranges")[c(3)],
                   WTB = brewer.pal(8, "Oranges")[c(5)],
                   WTC = brewer.pal(8, "Oranges")[c(7)],
                   HET = brewer.pal(8, "Greens")[c(3)],
                   KI_13 = brewer.pal(8, "Greens")[c(6)],
                   KI_27 = brewer.pal(8, "Greens")[c(7)],
                   KI_29 = brewer.pal(8, "Greens")[c(8)],
                   Del_16bp = brewer.pal(8, "Purples")[c(3)],
                   Del_233 = brewer.pal(8, "Purples")[c(6)],
                   Del_235 = brewer.pal(8, "Purples")[c(7)],
                   Del_287 = brewer.pal(8, "Purples")[c(8)]),
    time_point = c(h24 = brewer.pal(9, "Greys")[c(3)],
                   h48 = brewer.pal(9, "Greys")[c(5)],
                   h72 = brewer.pal(9, "Greys")[c(7)],
                   h96 = brewer.pal(9, "Greys")[c(9)]),
    Genotype = c('G.G' = 'black',
                 'A.G' = brewer.pal(8, "Paired")[c(3)],
                 'A.A' =  brewer.pal(8, "Paired")[c(4)],
                 Del16 = brewer.pal(8, "Paired")[c(5)],
                 Del80 = brewer.pal(8, "Paired")[c(6)]),
    seurat_cluster = c('1' = brewer.pal(9, "YlOrRd")[c(7)],
                       '2' = brewer.pal(9, "Blues")[c(5)],
                       '3' = brewer.pal(9, "RdPu")[c(6)],
                       '4' = brewer.pal(9, "Blues")[c(4)],
                       '5' = brewer.pal(9, "Greens")[c(6)],
                       '6' = brewer.pal(9, "YlOrRd")[c(6)],
                       '7' = brewer.pal(9, "Blues")[c(3)],
                       '8' = brewer.pal(9, "RdPu")[c(5)],
                       '9' = brewer.pal(9, "Greens")[c(5)],
                       '10' = brewer.pal(9, "Greens")[c(4)],
                       '11' = brewer.pal(9, "Blues")[c(2)],
                       '12' = brewer.pal(9, "RdPu")[c(4)],
                       '13' = brewer.pal(9, "RdPu")[c(3)]),
    GeneClass = c('CUX1' = 'black',
                  'Marker_genes' = 'firebrick1',
                  'MEP_progenitor' =  brewer.pal(8, "Dark2")[2],
                  'Megakaryocyte' = brewer.pal(8, "Dark2")[3],
                  'Platelet' = brewer.pal(8, "Dark2")[4],
                  'Platelet_volume'= brewer.pal(8, "Dark2")[8],
                  'Erythrocyte' = brewer.pal(8, "Dark2")[6],
                  'Mitosis' = brewer.pal(8, "BrBG")[1],
                  'ANAPC_genes' = brewer.pal(8, "BrBG")[3],
                  'CUX1_target_genes' = brewer.pal(8, "Dark2")[1],
                  'PU1_target_genes' = brewer.pal(8, "Accent")[3],
                  'RUNX1_target_genes' = brewer.pal(8, "Accent")[2],
                  'GATA_target_genes' = brewer.pal(8, "Accent")[1],
                  'AML_myeloid' ='gray'))
  
  
  
  cat("ann_colors_0\n")
  cat(str(ann_colors))
  cat("\n")
  
  cat("annotation_col_REMEMBER\n")
  cat(str(annotation_col))
  cat("\n")
  cat(sprintf(as.character(row.names(annotation_col))))
  cat("\n")
  
  #### Actual heatmap ----
  
  heatmap<-pheatmap(lcpm_subset, display_numbers = FALSE, number_format = "%.1e",
                    angle_col = "0",clustering_method="ward.D2",
                    fontsize_row = 1, 
                    fontsize_col = 6,
                    fontsize= 6,
                    breaks=seq(-2,2,length.out=101),
                    color=colorRampPalette(c("blue","white","red"))(100),
                    scale="row",
                    cluster_cols=TRUE,
                    cluster_rows=TRUE,
                    border_color='black',
                    treeheight_row=70, treeheight_col=70, cutree_cols=7,
                    annotation_col = annotation_col, annotation_row = annotation_row,
                    annotation_colors = ann_colors,
                    show_colnames=FALSE)
  
  
  svgname<-paste(paste("Heatmap_ATAC_TSS_and_Linked_Peaks","ALL_genotypes", sep='_'),".svg",sep='')
  
  setwd(path_graphs)
  svglite(svgname, width = 13, height = 13)
  print(heatmap)
  dev.off()
  
  
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
    make_option(c("--CUX1"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ANAPC_genes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ery_genes_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Master_peak_file_with_SNP_numbered_peaks"), type="character", default=NULL, 
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
    make_option(c("--ATAC_matrix"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--metadata"), type="character", default=NULL, 
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
  
  heatmap_WT(opt)
  heatmap_ALL(opt)

  
}


###########################################################################

system.time( main() )