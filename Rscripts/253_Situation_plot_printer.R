
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
  
  
  #### READ and transform seurat_cluster ----
  
  seurat_cluster = opt$seurat_cluster
  
  cat("seurat_cluster_\n")
  cat(sprintf(as.character(seurat_cluster)))
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
  
  #### Master Peak file
  
 
  #### READ Seurat_object ----
  
  Seurat_object<-readRDS(file=opt$Seurat_object)
  

  # cat("Seurat_object_0\n")
  # cat(str(Seurat_object))
  # cat("\n")
  
  Seurat_object$seurat_clusters_G<-interaction(Seurat_object$seurat_clusters,Seurat_object$Genotype)
  
  Seurat_object$seurat_clusters_G_by_time_point<-interaction(Seurat_object$time_point,Seurat_object$seurat_clusters_G)
  
  metadata<-Seurat_object[[]]
  
  cat("metadata_0\n")
  cat(str(metadata))
  cat("\n")
  
  
  ############### Subset object by variable of interest-----------------
  
  DefaultAssay(Seurat_object) <- 'ATAC'
  peaks <- GetAssayData(Seurat_object,slot="counts")
  
  cell_list=data.frame(colnames(peaks))
  colnames(cell_list)<-'Cell_ID'
  head(cell_list)
  dim(cell_list)
  
  metadata$Cell_ID<-cell_list$Cell_ID
  # str(metadata)
  cells_in_diff_group_1<-metadata[which(metadata$seurat_clusters == '1'),]
  
  cat("cells_in_diff_group_1\n")
  cat(str(cells_in_diff_group_1))
  cat("\n")
  
  cells_in_diff_group_2<-metadata[which(metadata$seurat_clusters == '2'),]
  
  cat("cells_in_diff_group_2\n")
  cat(str(cells_in_diff_group_2))
  cat("\n")
  
  cells_in_diff_group_3<-metadata[which(metadata$seurat_clusters == '3'),]
  
  cat("cells_in_diff_group_3\n")
  cat(str(cells_in_diff_group_3))
  cat("\n")
  
  cells_in_diff_group_4<-metadata[which(metadata$seurat_clusters == '4'),]
  
  cat("cells_in_diff_group_4\n")
  cat(str(cells_in_diff_group_4))
  cat("\n")
  
  cells_in_diff_group_5<-metadata[which(metadata$seurat_clusters == '5'),]
  
  cat("cells_in_diff_group_5\n")
  cat(str(cells_in_diff_group_5))
  cat("\n")
  
  cells_in_diff_group_6<-metadata[which(metadata$seurat_clusters == '6'),]
  
  cat("cells_in_diff_group_6\n")
  cat(str(cells_in_diff_group_6))
  cat("\n")
  
  cells_in_diff_group_7<-metadata[which(metadata$seurat_clusters == '7'),]
  
  cat("cells_in_diff_group_7\n")
  cat(str(cells_in_diff_group_7))
  cat("\n")
  
  cells_in_diff_group_8<-metadata[which(metadata$seurat_clusters == '8'),]
  
  cat("cells_in_diff_group_8\n")
  cat(str(cells_in_diff_group_8))
  cat("\n")
  
  cells_in_diff_group_9<-metadata[which(metadata$seurat_clusters == '9'),]
  
  cat("cells_in_diff_group_9\n")
  cat(str(cells_in_diff_group_9))
  cat("\n")
  
  cells_in_diff_group_10<-metadata[which(metadata$seurat_clusters == '10'),]
  
  cat("cells_in_diff_group_10\n")
  cat(str(cells_in_diff_group_10))
  cat("\n")
  
  cells_in_diff_group_11<-metadata[which(metadata$seurat_clusters == '11'),]
  
  cat("cells_in_diff_group_11\n")
  cat(str(cells_in_diff_group_11))
  cat("\n")
  
  cells_in_diff_group_12<-metadata[which(metadata$seurat_clusters == '12'),]
  
  cat("cells_in_diff_group_12\n")
  cat(str(cells_in_diff_group_12))
  cat("\n")
  
  cells_in_diff_group_13<-metadata[which(metadata$seurat_clusters == '13'),]
  
  cat("cells_in_diff_group_13\n")
  cat(str(cells_in_diff_group_13))
  cat("\n")
  
  
  ### Subdivide the multiome object in different subsets containing the cells from each genotype ----
  
  adata_diff_group_1<-subset(
    Seurat_object,
    cells = cells_in_diff_group_1$Cell_ID)
  
  adata_diff_group_2<-subset(
    Seurat_object,
    cells = cells_in_diff_group_2$Cell_ID)
  
  adata_diff_group_3<-subset(
    Seurat_object,
    cells = cells_in_diff_group_3$Cell_ID)
  
  adata_diff_group_4<-subset(
    Seurat_object,
    cells = cells_in_diff_group_4$Cell_ID)
  
  adata_diff_group_5<-subset(
    Seurat_object,
    cells = cells_in_diff_group_5$Cell_ID)
  
  adata_diff_group_6<-subset(
    Seurat_object,
    cells = cells_in_diff_group_6$Cell_ID)
  
  adata_diff_group_7<-subset(
    Seurat_object,
    cells = cells_in_diff_group_7$Cell_ID)
  
  adata_diff_group_8<-subset(
    Seurat_object,
    cells = cells_in_diff_group_8$Cell_ID)
  
  adata_diff_group_9<-subset(
    Seurat_object,
    cells = cells_in_diff_group_9$Cell_ID)
  
  adata_diff_group_10<-subset(
    Seurat_object,
    cells = cells_in_diff_group_10$Cell_ID)
  
  adata_diff_group_11<-subset(
    Seurat_object,
    cells = cells_in_diff_group_11$Cell_ID)
  
  adata_diff_group_12<-subset(
    Seurat_object,
    cells = cells_in_diff_group_12$Cell_ID)
  
  adata_diff_group_13<-subset(
    Seurat_object,
    cells = cells_in_diff_group_13$Cell_ID)
  
  #### vector_fill for genotypes ----
  
  vector_fill<-brewer.pal(length(levels(metadata$Genotype)),"Dark2")
  
  cat("vector_fill_0\n")
  cat(str(vector_fill))
  cat("\n")
  
  #### READ Genomic_features ----

  Genomic_features<-readRDS(file=opt$Genomic_features)
  
  cat("Genomic_features_0\n")
  cat(str(Genomic_features))
  cat("\n")
  cat(str(unique(Genomic_features$stable_id)))
  cat("\n")
  
  #### Read the Input regions ----
  
  Input_regions<-readRDS(file=opt$Input_regions)
  
  cat("Input_regions_0\n")
  cat(str(Input_regions))
  cat("\n")
  
  #### READ Master_peak_file_with_SNP ----

  Master_peak_file_with_SNP_numbered_peaks<-readRDS(file=opt$Master_peak_file_with_SNP)
  
  cat("Master_peak_file_with_SNP_numbered_peaks_0\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks))
  cat("\n")

  
  Master_peak_file_with_SNP_numbered_peaks$chr<-gsub("^chr","",Master_peak_file_with_SNP_numbered_peaks$chr)
  
  Master_peak_file_with_SNP_numbered_peaks_SAVE<-Master_peak_file_with_SNP_numbered_peaks
  
  cat("Master_peak_file_with_SNP_numbered_peaks_SAVE_0\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks_SAVE))
  cat("\n")
  
  colnames(Master_peak_file_with_SNP_numbered_peaks)[which(colnames(Master_peak_file_with_SNP_numbered_peaks) == 'Peak_ID')]<-'stable_id'
  
  cat("Master_peak_file_with_SNP_numbered_peaks_1\n")
  cat(str(Master_peak_file_with_SNP_numbered_peaks))
  cat("\n")
  
  #### Read metadata GeneEXP ----
  
  
  metadata_GeneEXP<-readRDS(file=opt$metadata_GeneEXP)
  
  # cat("metadata_GeneEXP_0\n")
  # cat(str(metadata_GeneEXP))
  # cat("\n")
  
  
  metadata_GeneEXP$clone_line<-as.character(metadata_GeneEXP$clone_line)
  
  metadata_GeneEXP$clone_line<-gsub("^chrGFP_","",metadata_GeneEXP$clone_line)
  
  metadata_GeneEXP$clone_line<-factor(metadata_GeneEXP$clone_line,
                                      levels=c("WTA","WTB","WTC","HET","KI_13","KI_27","KI_29","Del_16bp","Del_233","Del_235","Del_287"),
                                      ordered=T)
  
  cat("metadata_GeneEXP_1\n")
  cat(str(metadata_GeneEXP))
  cat("\n")
  
  
  
  metadata_GeneEXP<-metadata_GeneEXP[order(metadata_GeneEXP$seurat_cluster),]
  
  cat("metadata_GeneEXP_ordered\n")
  cat(str(metadata_GeneEXP))
  cat("\n")
  
  column_order<-metadata_GeneEXP$columns_for_matrix
  
  # cat("column_order_0\n")
  # cat(sprintf(as.character(column_order)))
  # cat("\n")
  
  metadata_GeneEXP_subset<-droplevels(unique(metadata_GeneEXP[,c(which(colnames(metadata_GeneEXP) == 'clone_line'),
                                                                 which(colnames(metadata_GeneEXP) == 'Genotype'))]))
  
  cat("metadata_GeneEXP_subset_0\n")
  cat(str(metadata_GeneEXP_subset))
  cat("\n")
  
  #### Read the GeneEXP file ----
  
  df_for_matrix<-as.data.frame(fread(file=opt$GeneEXP, sep="\t", header=T), stringsAsFactors=F)
  
  # cat("df_for_matrix_0\n")
  # cat(str(df_for_matrix))
  # cat("\n")
  
  indx.int<-c(which(colnames(df_for_matrix) == 'Symbol'),which(colnames(df_for_matrix)%in%column_order))
  
  
  # cat("indx.int\n")
  # cat(str(indx.int))
  # cat("\n")
  
  df_for_matrix_subset<-unique(df_for_matrix[,indx.int])
  
  # cat("df_for_matrix_subset_0\n")
  # cat(str(df_for_matrix_subset))
  # cat("\n")
  
  indx.reorder<-c(which(colnames(df_for_matrix_subset) == 'Symbol'),which(colnames(df_for_matrix_subset)%in%column_order))
  
  
  # cat("indx.reorder\n")
  # cat(str(indx.reorder))
  # cat("\n")
  
  
  # df_for_matrix_subset_reordered[indx.reorder]
  # 
  
  
  df_for_matrix_subset.dt<-setDT(df_for_matrix_subset)
  
  setcolorder(df_for_matrix_subset.dt, c('Symbol',column_order))
  
  df_for_matrix_subset_reordered<-as.data.frame(df_for_matrix_subset.dt, stringsAsFactors=F)
  
  # cat("df_for_matrix_subset_reordered_0\n")
  # cat(str(df_for_matrix_subset_reordered))
  # cat("\n")
  
  
  GeneEXP_matrix<-as.matrix(df_for_matrix_subset_reordered[,-c(which(colnames(df_for_matrix_subset_reordered) == 'Symbol'))])
  row.names(GeneEXP_matrix)<-df_for_matrix_subset_reordered$Symbol
  
  # cat("GeneEXP_matrix_0\n")
  # cat(str(GeneEXP_matrix))
  # cat("\n")
  
  
  
  
  
  lcpm_GeneEXP <- cpm(GeneEXP_matrix, log=TRUE)
  
  cat("lcpm_GeneEXP_0\n")
  cat(str(lcpm_GeneEXP))
  cat("\n")
 
  
  #### Read metadata ATAC ----
  
  
  metadata_ATAC<-readRDS(file=opt$metadata_ATAC)
  
  # cat("metadata_ATAC_0\n")
  # cat(str(metadata_ATAC))
  # cat("\n")
  
  
  metadata_ATAC$clone_line<-as.character(metadata_ATAC$clone_line)
  
  metadata_ATAC$clone_line<-gsub("^chrGFP_","",metadata_ATAC$clone_line)
  
  metadata_ATAC$clone_line<-factor(metadata_ATAC$clone_line,
                                   levels=c("WTA","WTB","WTC","HET","KI_13","KI_27","KI_29","Del_16bp","Del_233","Del_235","Del_287"),
                                   ordered=T)
  
  # cat("metadata_ATAC_1\n")
  # cat(str(metadata_ATAC))
  # cat("\n")
  
  
  
  metadata_ATAC<-metadata_ATAC[order(metadata_ATAC$seurat_cluster),]
  
  cat("metadata_ATAC_ordered\n")
  cat(str(metadata_ATAC))
  cat("\n")
  
  column_order<-metadata_ATAC$columns_for_matrix
  
  # cat("column_order_0\n")
  # cat(sprintf(as.character(column_order)))
  # cat("\n")
  
  metadata_ATAC_subset<-droplevels(unique(metadata_ATAC[,c(which(colnames(metadata_ATAC) == 'clone_line'),
                                                           which(colnames(metadata_ATAC) == 'Genotype'))]))
  
  cat("metadata_ATAC_subset_0\n")
  cat(str(metadata_ATAC_subset))
  cat("\n")
  
  #### Read the ATAC file -----
  
  df_for_matrix<-as.data.frame(fread(file=opt$ATAC, sep="\t", header=T), stringsAsFactors=F)
  
  cat("df_for_matrix_0\n")
  cat(str(df_for_matrix))
  cat("\n")
  
  df_for_matrix<-merge(df_for_matrix,
                       Master_peak_file_with_SNP_numbered_peaks_SAVE,
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
  
  lcpm_ATAC <- cpm(ATAC_matrix, log=TRUE)
  
  cat("lcpm_ATAC_0\n")
  cat(str(lcpm_ATAC))
  cat("\n")
  #### Get the distribution of the links values ----
  
  indx.Links<-which(Genomic_features$feature_type%in%c('Link'))
  
  Links_df<-Genomic_features[indx.Links,]
  
  max_abs_value<-abs(max(Links_df$zscore[!is.na(Links_df$zscore)]))
  min_abs_value<-abs(min(Links_df$zscore[!is.na(Links_df$zscore)]))
  
  if(max_abs_value > min_abs_value)
  {
    step<-round(abs(max_abs_value--1*max_abs_value)/4,2)
    
    breaks_zscore<-unique(sort(round(c(0,max_abs_value,seq(-1*max_abs_value,max_abs_value, by=step)),1)))
    labels_zscore<-as.character(breaks_zscore)
    
  }else{
    
    step<-round(abs(min_abs_value--1*min_abs_value)/4,2)
    
    breaks_zscore<-unique(sort(round(c(0,min_abs_value,seq(-1*min_abs_value,min_abs_value, by=step)),1)))
    labels_zscore<-as.character(breaks_zscore)
    
  }# max_abs_value > min_abs_value
  
  
  
  
  cat("Links_df_0\n")
  cat(str(Links_df))
  cat("\n")
  cat(str(unique(Links_df$stable_id)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Links_df$Minus_logpval)))))
  cat("\n")
  cat(sprintf(as.character(summary(Links_df$Minus_logpval))))
  cat("\n")
  cat(sprintf(as.character(names(summary(Links_df$zscore)))))
  cat("\n")
  cat(sprintf(as.character(summary(Links_df$zscore))))
  cat("\n")
  
  cat("breaks_zscore\n")
  cat(sprintf(as.character(breaks_zscore)))
  cat("\n")
  cat("labels_zscore\n")
  cat(sprintf(as.character(labels_zscore)))
  cat("\n")
  
  
  #### THE LOOP----
  
  
  DEBUG<-0
  
  cluster_levels<-levels(metadata$seurat_clusters)
  
  cat("cluster_levels_0\n")
  cat(str(cluster_levels))
  cat("\n")
  
  indx_RUSH<-which(Input_regions$Symbol%in%c('EZH2','MOCOS','PDE10A','CUX1','GUSB','SMAD2'))

  # indx_RUSH<-which(Input_regions$Symbol%in%c('SMAD2'))
  
  cat("indx_RUSH_0\n")
  cat(str(indx_RUSH))
  cat("\n")
  
  # [indx_RUSH,]
  
  Input_regions_RUSH<-rbind(Input_regions[indx_RUSH,],Input_regions[-indx_RUSH,])
  # Input_regions_RUSH<-rbind(Input_regions[indx_RUSH,])
  
  
  cat("Input_regions_RUSH_0\n")
  cat(str(Input_regions_RUSH))
  cat("\n")
  
  
  
  START<-1
  
  for(i in START:dim(Input_regions_RUSH)[1])
  {
    Input_regions_RUSH_sel<-Input_regions_RUSH[i,]
    
    if(DEBUG == 1)
    {
      cat("Input_regions_RUSH_sel_0\n")
      cat(str(Input_regions_RUSH_sel))
      cat("\n")
      cat(str(unique(Input_regions_RUSH_sel$Symbol)))
      cat("\n")
    }
    
    Symbol_sel<-unique(Input_regions_RUSH_sel$Symbol)
    Display_genes_sel<-unique(unlist(strsplit(Input_regions_RUSH_sel$Display_genes, split=";")))
    Region_sel<-Input_regions_RUSH_sel$region_name
    chr_sel<-Input_regions_RUSH_sel$chr
    start_sel<-Input_regions_RUSH_sel$start
    end_sel<-Input_regions_RUSH_sel$end
    
    Display_peaks_sel<-unique(unlist(strsplit(Input_regions_RUSH_sel$Display_peaks, split="\\|")))
    
    indx.na<-which(Display_peaks_sel == 'NA')
    
    if(DEBUG == 1)
    {
      cat("indx.na_0\n")
      cat(str(indx.na))
      cat("\n")
    }
    
    if(length(indx.na) >0){
      
      Display_peaks_sel<-Display_peaks_sel[-indx.na]
      
    }else{
      
      Display_peaks_sel<-Display_peaks_sel
    }# length(indx.na) >0
    
    
    ###### Genomic features plot ------ 
    
    Genomic_features_sel<-unique(Genomic_features[which(Genomic_features$region_string == Region_sel &
                                                          Genomic_features$Symbol == Symbol_sel),])
    
    indx.gene.pos.strand<-which(Genomic_features_sel$feature_type == 'gene' & Genomic_features_sel$strand== '+')
    indx.gene.neg.strand<-which(Genomic_features_sel$feature_type == 'gene' & Genomic_features_sel$strand == '-')
    
    
    total_n_genes<-length(which(Genomic_features_sel$feature_type == 'gene'))
    
    random_vec <- round(runif(n=total_n_genes, min=-0.5, max=0.5),3)        
    
    vector_colors_1<-rev(brewer.pal(9, "Blues")[c(3,5,7,9)])
    vector_colors_2<-rev(brewer.pal(9, "Greens")[c(3,5,7,9)])
    vector_colors_3<-rev(brewer.pal(9, "Purples")[c(3,5,7,9)])
    vector_colors_4<-rev(brewer.pal(9, "Reds")[c(3,5,7,9)])
    vector_colors_5<-rev(brewer.pal(9, "Greys")[c(3,5,7,9)])
    
    vector_colors_DA_features<-c(vector_colors_1,vector_colors_2,vector_colors_3,vector_colors_4,vector_colors_5)
    
    indx.features<-which(Genomic_features_sel$feature_type%in%c('Promoter','Enhancer','TFBS','OpenChromatin','CTCF'))
    indx.Peaks<-which(Genomic_features_sel$feature_type%in%c('Peak'))
    indx.Links<-which(Genomic_features_sel$feature_type%in%c('Link'))
    
    if(DEBUG == 1)
    {
      cat("Genomic_features_sel_0\n")
      cat(str(Genomic_features_sel))
      cat("\n")
      cat(str(unique(Genomic_features_sel$stable_id)))
      cat("\n")
      cat(sprintf(as.character(names(summary(Genomic_features_sel$ymin)))))
      cat("\n")
      cat(sprintf(as.character(summary(Genomic_features_sel$ymin))))
      cat("\n")
      cat(sprintf(as.character(names(summary(Genomic_features_sel$ymax)))))
      cat("\n")
      cat(sprintf(as.character(summary(Genomic_features_sel$ymax))))
      cat("\n")  
      cat(sprintf(as.character(names(summary(Genomic_features_sel$strand)))))
      cat("\n")
      cat(sprintf(as.character(summary(Genomic_features_sel$strand))))
      cat("\n") 
      
      cat("Positive_strand_genes\n")
      cat(str(indx.gene.pos.strand))
      cat("\n")
      
      cat("Negative_strand_genes\n")
      cat(str(indx.gene.neg.strand))
      cat("\n")
      
      cat("total_n_genes\n")
      cat(str(total_n_genes))
      cat("\n")
      
      cat("random_vec\n")
      cat(str(random_vec))
      cat("\n")
      
      cat("vector_colors_DA_features\n")
      cat(str(vector_colors_DA_features))
      cat("\n")
      
      cat("indx.features\n")
      cat(str(indx.features))
      cat("\n")
      
      cat("indx.Peaks\n")
      cat(str(indx.Peaks))
      cat("\n")
      
      cat("indx.Links\n")
      cat(str(indx.Links))
      cat("\n")
      
    }
    
    Peaks_sel<-Genomic_features_sel[indx.Peaks,]
    
    if(DEBUG == 1)
    {
      #  cat("Peaks_sel_0\n")
      # cat(str(Peaks_sel))
      # cat("\n")
      #  cat(sprintf(as.character(unique(Peaks_sel$value_string))))
      # cat("\n")
    }
    
    
    Genomic_features_sel$stranded_start<-NA
    Genomic_features_sel$stranded_end<-NA
    
    Genomic_features_sel$ymin[indx.gene.pos.strand]<-Genomic_features_sel$ymin[indx.gene.pos.strand]+random_vec
    Genomic_features_sel$ymax[indx.gene.pos.strand]<-Genomic_features_sel$ymin[indx.gene.pos.strand]+0.2
    
    Genomic_features_sel$stranded_start[indx.gene.pos.strand]<-Genomic_features_sel$start[indx.gene.pos.strand]
    Genomic_features_sel$stranded_end[indx.gene.pos.strand]<-Genomic_features_sel$end[indx.gene.pos.strand]
    
    Genomic_features_sel$ymin[indx.gene.neg.strand]<-Genomic_features_sel$ymin[indx.gene.neg.strand]+random_vec
    Genomic_features_sel$ymax[indx.gene.neg.strand]<-Genomic_features_sel$ymin[indx.gene.neg.strand]+0.2
    
    Genomic_features_sel$stranded_start[indx.gene.neg.strand]<-Genomic_features_sel$end[indx.gene.neg.strand]
    Genomic_features_sel$stranded_end[indx.gene.neg.strand]<-Genomic_features_sel$start[indx.gene.neg.strand]
    
    
    breaks.y<-unique(sort(c(max(Genomic_features_sel$ymax), seq(0, max(Genomic_features_sel$ymax), by=1))))
    labels.y<-as.character(breaks.y)
    
    if(DEBUG == 1)
    {
      cat(sprintf(as.character(labels.y)))
      cat("\n")
    }
    
    ### Master Peak_file_sel -----------------------------------------
    
    Master_peak_file_with_SNP_numbered_peaks_sel<-Master_peak_file_with_SNP_numbered_peaks[which(Master_peak_file_with_SNP_numbered_peaks$Peak_Number%in%Display_peaks_sel),]
    
    Master_peak_file_with_SNP_numbered_peaks_sel<-merge(Master_peak_file_with_SNP_numbered_peaks_sel,
                                                        Peaks_sel,
                                                        by=c('chr','start','end','stable_id'))
    
    
    
    peak_starts<-Master_peak_file_with_SNP_numbered_peaks_sel$start
    peak_ends<-Master_peak_file_with_SNP_numbered_peaks_sel$end
    
    ##### GR object for Peaks_selected ------------
    
    gr_Peaks_selected <- GRanges(
      seqnames = as.character(Master_peak_file_with_SNP_numbered_peaks_sel$chr), 
      color="#707070",
      ranges=IRanges(
        start=peak_starts,
        end=peak_ends,
        name=Master_peak_file_with_SNP_numbered_peaks_sel$Peak_ID))
    
    if(DEBUG == 1)
    {
      # cat("gr_Peaks_selected_0\n")
      # cat(str(gr_Peaks_selected))
      # cat("\n")
    }
    
    if(DEBUG == 1)
    {
      cat("Master_peak_file_with_SNP_numbered_peaks_sel_0\n")
      cat(str(Master_peak_file_with_SNP_numbered_peaks_sel))
      cat("\n")
      cat(str(unique(Master_peak_file_with_SNP_numbered_peaks_sel$Peak_Number)))
      cat("\n")
    }
    
    Elements_selected<-unique(unlist(strsplit(Master_peak_file_with_SNP_numbered_peaks_sel$value_string, split=";")))
    
    indx.ENSR<-grep("^ENSR",Elements_selected)
    
    Elements_selected_reg<-Elements_selected[indx.ENSR]
    
    Elements_selected_reg<-gsub("\\|.+$","",Elements_selected_reg)
    
    if(DEBUG == 1)
    {
      cat("Elements_selected_0\n")
      cat(str(Elements_selected))
      cat("\n")
      
      cat("Elements_selected_reg_0\n")
      cat(sprintf(as.character(Elements_selected_reg)))
      cat("\n")
      
    }
    
    
    #### Actual graph_Genomic_features ----
    
    graph_Genomic_features<-ggplot()+
      scale_x_continuous(name=NULL,breaks=c(start_sel,end_sel),
                         labels=as.character(c(start_sel,end_sel)),
                         limits=c(start_sel-1,end_sel+1))+
      scale_y_continuous(name=NULL, 
                         breaks=c(0.5,1,2,3,4),
                         labels=c(paste('Links','Gene Exp', sep="\n"),paste('ATAC','Peaks', sep="\n"),paste('Regulatory','Elements', sep="\n"),paste('Genes','strand -', sep="\n"),paste('Genes','strand +', sep="\n")),
                         limits=c(breaks.y[1],breaks.y[length(breaks.y)]+0.3))+
      geom_segment(data=Genomic_features_sel[which(Genomic_features_sel$feature_type == 'gene'),],
                   aes(y=ymin,
                       yend=ymin,
                       x=stranded_start,
                       xend=stranded_end),color='black',size=1,
                   lineend = "round", # See available arrow types in example above
                   linejoin = "bevel",
                   size = 2, 
                   arrow = arrow(length=unit(1,"mm")))+
      geom_rect(data=Genomic_features_sel[which(Genomic_features_sel$feature_type%in%c('Promoter','Enhancer','TFBS','OpenChromatin','CTCF')),],
                aes(ymin=ymin,
                    ymax=ymax,
                    xmin=start,
                    xmax=end,
                    fill=interaction),color='black',size=0.05)+    
      scale_fill_manual(name=paste("K-562 Genomic","Feature activity",sep="\n"),values=vector_colors_DA_features, drop=F)+   
      geom_rect(data=Peaks_sel,
                aes(ymin=ymin,
                    ymax=ymax,
                    xmin=start,
                    xmax=end),color='black',size=0.05)+  
      geom_curve(data=Genomic_features_sel[indx.Links,],
                 aes(x=start,
                     xend=end,
                     y=ymax,
                     yend=ymax,
                     color=zscore),
                 curvature = 0.2, size=0.75)+
      scale_color_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        breaks=breaks_zscore,labels=labels_zscore,
        limits=c(breaks_zscore[1],breaks_zscore[length(breaks_zscore)]),name=paste('Links',"zscore",sep="\n"),na.value = "gray")+
      geom_text_repel(data=Genomic_features_sel[which(Genomic_features_sel$gene_name%in%Display_genes_sel |
                                                        Genomic_features_sel$stable_id%in%Display_genes_sel),],
                      aes(x=stranded_start,
                          y=ymin,
                          label=paste(stable_id,gene_name,sep='|')),size=2,ylim=c(3,4), family="sans",color='black', fontface='plain',box.padding = unit(0.5, "lines"),                      
                      segment.size  = 0.2,max.overlaps = Inf,segment.color = "black")+     
      geom_text_repel(data=Genomic_features_sel[which(Genomic_features_sel$feature_type%in%c('Promoter','Enhancer','TFBS','OpenChromatin','CTCF') &
                                                        Genomic_features_sel$stable_id%in%Elements_selected_reg),],
                      aes(x=start,
                          y=ymax,                 
                          label=paste(stable_id,feature_type,activity,sep='|')),size=2,ylim=c(2,3), family="sans",color='black', fontface='plain',box.padding = unit(0.5, "lines"),                    
                      segment.size  = 0.2,max.overlaps = Inf,
                      segment.color = "black")+    
      geom_text_repel(data=Master_peak_file_with_SNP_numbered_peaks_sel,
                      aes(x=start,
                          y=ymin,
                          label=Peak_Number),size=2,ylim=c(1,2), family="sans",color='black', fontface='plain',box.padding = unit(0.5, "lines"),                    
                      segment.size  = 0.2,max.overlaps = Inf,
                      segment.color = "black")+        
      theme_classic()+
      theme(plot.title=element_text(size=6, color="black", family="sans"),
            axis.title.y=element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
            axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
            axis.line.x = element_line(size = 0.4),
            axis.ticks.x = element_line(size = 0.4),
            axis.ticks.y = element_line(size = 0.4),
            axis.line.y = element_line(size = 0.4))+
      theme(legend.title = element_text(size=6, color="black", family="sans"),
            legend.text = element_text(size=6, color="black", family="sans"),
            legend.key.size = unit(0.25, 'cm'), #change legend key size
            legend.key.height = unit(0.25, 'cm'), #change legend key height
            legend.key.width = unit(0.25, 'cm'), #change legend key width
            legend.position="right")+        
      ggeasy::easy_center_title()
    
    # [!is.na(Peaks_sel$value_string),]
    # geom_text_repel(data=Genomic_features_sel[which(Genomic_features_sel$gene_name%in%Display_genes_sel),],
    
    if(Symbol_sel == 'CUX1'){  
      graph_Genomic_features <-graph_Genomic_features +
        geom_vline(xintercept=c(101856650,peak_starts,peak_ends), color="#363636",linetype="dashed",linewidth=0.1)
    }else{
      
      graph_Genomic_features <-graph_Genomic_features +
        geom_vline(xintercept=c(peak_starts,peak_ends), color="#363636",linetype="dashed",linewidth=0.1)
    }
    
    for(k in 1:length(cluster_levels))
    {
      seurat_clusters_sel<-cluster_levels[k]
      
      cat("-------------------------------------->\t")
      cat(sprintf(as.character(paste(i,Region_sel,Display_genes_sel,Display_peaks_sel, collapse=" "))))
      cat("\t")
      cat(sprintf(as.character(paste(k,seurat_clusters_sel, collapse=" "))))
      cat("\n")
      
      
      ##### define which Seurat_object object I am going to use ----
      
      
      selected_adata<-NULL
      
      if(seurat_clusters_sel == '1')
      {
        selected_adata<-adata_diff_group_1
      }else{
        if(seurat_clusters_sel == '2'){
          
          selected_adata<-adata_diff_group_2
          
        }else{
          
          if(seurat_clusters_sel == '3'){
            
            selected_adata<-adata_diff_group_3
            
          }else{
            
            if(seurat_clusters_sel == '4'){
              
              selected_adata<-adata_diff_group_4
              
            }else{
              
              if(seurat_clusters_sel == '5'){
                
                selected_adata<-adata_diff_group_5
                
              }else{
                if(seurat_clusters_sel == '6'){
                  
                  selected_adata<-adata_diff_group_6
                  
                }else{
                  if(seurat_clusters_sel == '7'){
                    
                    selected_adata<-adata_diff_group_7
                    
                  }else{
                    if(seurat_clusters_sel == '8'){
                      
                      selected_adata<-adata_diff_group_8
                      
                    }else{
                      if(seurat_clusters_sel == '9'){
                        
                        selected_adata<-adata_diff_group_9
                        
                      }else{
                        if(seurat_clusters_sel == '10'){
                          
                          selected_adata<-adata_diff_group_10
                          
                        }else{
                          if(seurat_clusters_sel == '11'){
                            
                            selected_adata<-adata_diff_group_11
                            
                          }else{
                            
                            if(seurat_clusters_sel == '12'){
                              
                              selected_adata<-adata_diff_group_12
                              
                            }else{
                              
                              if(seurat_clusters_sel == '13'){
                                
                                selected_adata<-adata_diff_group_13
                                
                              }else{
                                
                                cat("ERROR_no_adata_subject\n")
                                break  
                                
                              }# seurat_clusters_sel == '13'
                            }# seurat_clusters_sel == '12'
                          }#seurat_clusters_sel == '11'
                        }#seurat_clusters_sel == '10'
                      }#seurat_clusters_sel == '9'
                    }# seurat_clusters_sel == '8'
                  }#seurat_clusters_sel == '7'
                }#seurat_clusters_sel == '6'
              }#seurat_clusters_sel == '5'
            }#seurat_clusters_sel == '4'
          }#seurat_clusters_sel == '3'
        }#seurat_clusters_sel == '2'
      }#seurat_clusters_sel == '1'
      
      selected_adata[['seurat_clusters']]<-droplevels(selected_adata[['seurat_clusters']])
      selected_adata[['seurat_clusters_G_by_time_point']]<-droplevels(selected_adata[['seurat_clusters_G_by_time_point']])
      
      vector_fill_expanded<-c(rep(vector_fill[1],4),rep(vector_fill[2],4),rep(vector_fill[3],4),rep(vector_fill[4],4),rep(vector_fill[5],4))
      
      if(DEBUG == 1)
      {
        cat("selected_adata_metadata\n")
        cat(str(selected_adata[['seurat_clusters']]))
        cat("\n")
        cat(str(selected_adata[['Genotype']]))
        cat("\n")
        cat(sprintf(as.character(names(summary(selected_adata[['Genotype']])))))
        cat("\n")
        cat(sprintf(as.character(summary(selected_adata[['Genotype']]))))
        cat("\n")
        cat(sprintf(as.character(vector_fill_expanded)))
        cat("\n")
      }
      
      ##### extract the region of coverage-----
      
      coverage_region<-paste(chr_sel,start_sel,end_sel,sep='-')
      
      if(DEBUG == 1)
      {
        cat("coverage_region_0\n")
        cat(sprintf(as.character(coverage_region)))
        cat("\n")
      }
      
      ###### Idents(selected_adata) & DefaultAssay(Seurat_object) <- "SCT"-------
      
      Idents(selected_adata) = 'seurat_clusters_G_by_time_point'
      DefaultAssay(selected_adata) <- "ATAC"
      
      ###### Cov plot ------    
      
      cov_plot <-CoveragePlot(
        object = selected_adata,
        region = coverage_region,        
        annotation = FALSE,
        peaks = FALSE,
        links=FALSE,
        region.highlight = gr_Peaks_selected,
      )+
        scale_fill_manual(values=vector_fill_expanded, drop=F)+
        theme_classic()+
        theme(axis.title.y=element_text(size=8, color="black", family="sans"),
              axis.title.x=element_blank(),
              axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
              axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
              axis.line.x = element_line(size = 0.4),
              axis.ticks.x = element_line(size = 0.4),
              axis.ticks.y = element_line(size = 0.4),
              axis.line.y = element_line(size = 0.4))+
        theme(legend.title =element_blank(),
              legend.text = element_text(size=6, color="black", family="sans"),
              legend.key.size = unit(0.45, 'cm'), #change legend key size
              legend.key.height = unit(0.45, 'cm'), #change legend key height
              legend.key.width = unit(0.45, 'cm'), #change legend key width
              legend.position="right")+        
        ggeasy::easy_center_title()
      
      if(Symbol_sel == 'CUX1'){  
        # cov_plot <-cov_plot +
        #     geom_vline(xintercept=101856650, color="#363636",linetype="dashed",linewidth=0.1)
      }else{
        
        # Do nothing
      }
      
      #### GeneEXP file ------
      
      metadata_GeneEXP_sel<-droplevels(metadata_GeneEXP[which(metadata_GeneEXP$seurat_cluster == seurat_clusters_sel),])
      
      if(DEBUG == 1)
      {
        cat("metadata_GeneEXP_sel_0\n")
        cat(str(metadata_GeneEXP_sel))
        cat("\n")
      }
      
      metadata_GeneEXP_sel<-metadata_GeneEXP_sel[order(metadata_GeneEXP_sel$Genotype, metadata_GeneEXP_sel$time_point),]
      
      if(DEBUG == 1)
      {
        cat("metadata_GeneEXP_sel_ORDERED\n")
        cat(str(metadata_GeneEXP_sel))
        cat("\n")
      }
      
      
      column_order<-metadata_GeneEXP_sel$columns_for_matrix
      
      if(DEBUG == 1)
      {
        cat("column_order_0\n")
        cat(sprintf(as.character(column_order)))
        cat("\n")
      }
      
      
      
      
      lcpm_GeneEXP_sel<-lcpm_GeneEXP[which(row.names(lcpm_GeneEXP)%in%Display_genes_sel),
                                     which(colnames(lcpm_GeneEXP)%in%column_order), drop=F]
      
      
      lcpm_GeneEXP_sel_df <- as.data.frame(lcpm_GeneEXP_sel %>%
                                             as_tibble())
      
      lcpm_GeneEXP_sel_df$Symbol<-as.character(row.names(lcpm_GeneEXP_sel))
      
      if(DEBUG == 1)
      {
        # cat("lcpm_GeneEXP_sel\n")
        # cat(str(lcpm_GeneEXP_sel))
        # cat("\n")
        
        #  #  cat("names(lcpm_GeneEXP_sel)\n")
        #  # cat(sprintf(as.character(names(lcpm_GeneEXP_sel))))
        #  # cat("\n")
        
        
        #  cat("row.names(lcpm_GeneEXP_sel)\n")
        #  cat(sprintf(as.character(row.names(lcpm_GeneEXP_sel))))
        #  cat("\n")
        #    cat("colnames(lcpm_GeneEXP_sel)\n")
        #  cat(sprintf(as.character(colnames(lcpm_GeneEXP_sel))))
        #  cat("\n")
        
        #  cat("lcpm_GeneEXP_sel_df\n")
        # cat(str(lcpm_GeneEXP_sel_df))
        # cat("\n")
        
      }
      
      lcpm_GeneEXP_sel_df.m<-melt(lcpm_GeneEXP_sel_df, id.vars='Symbol', variable.name='columns_for_matrix', value.name='cpm')
      
      if(DEBUG == 1)
      {
        cat("lcpm_GeneEXP_sel_df.m_0\n")
        cat(str(lcpm_GeneEXP_sel_df.m))
        cat("\n")
      }
      
      lcpm_GeneEXP_sel_df.m<-merge(lcpm_GeneEXP_sel_df.m,
                                   metadata_GeneEXP_sel,
                                   by='columns_for_matrix')
      
      if(DEBUG == 1)
      {
        cat("lcpm_GeneEXP_sel_df.m_1\n")
        cat(str(lcpm_GeneEXP_sel_df.m))
        cat("\n")
      }        
      
      lcpm_GeneEXP_sel_df.m.dt<-data.table(lcpm_GeneEXP_sel_df.m,
                                           key=c('Symbol','Genotype','time'))
      
      
      
      Mean_GeneEXP_df<-as.data.frame(lcpm_GeneEXP_sel_df.m.dt[,.(cpm=mean(cpm),
                                                                 sd=sd(cpm)), by=key(lcpm_GeneEXP_sel_df.m.dt)], stringsAsFactors=F)
      
      Mean_GeneEXP_df$cpm_max<-Mean_GeneEXP_df$cpm+Mean_GeneEXP_df$sd
      Mean_GeneEXP_df$cpm_min<-Mean_GeneEXP_df$cpm-Mean_GeneEXP_df$sd
      
      if(DEBUG == 1)
      {
        cat("--------------------------------------------------------------------------->Mean_GeneEXP_df_0\n")
        cat(str(Mean_GeneEXP_df))
        cat("\n")
      }
      
      Mean_GeneEXP_df.dt<-data.table(Mean_GeneEXP_df, key=c('Symbol','time'))
      
      summary_cpm<-unique(c(Mean_GeneEXP_df$cpm_max[!is.na(Mean_GeneEXP_df$cpm_max)],Mean_GeneEXP_df$cpm_min[!is.na(Mean_GeneEXP_df$cpm_min)]))
      
      Max_GeneEXP_cpm<-max(summary_cpm)
      Min_GeneEXP_cpm<-min(summary_cpm)
      
      step<-round((Max_GeneEXP_cpm-Min_GeneEXP_cpm)/4,2)
      
      if(step == 0)
      {
        step<-0.1
      }
      
      if(DEBUG == 1)
      {
        cat("Max_GeneEXP_cpm:\t")
        cat(sprintf(as.character(Max_GeneEXP_cpm)))
        cat("\n")
        cat("Min_GeneEXP_cpm:\t")
        cat(sprintf(as.character(Min_GeneEXP_cpm)))
        cat("\n")
        cat("step:\t")
        cat(sprintf(as.character(step)))
        cat("\n")
      }
      
      breaks_cpm<-unique(sort(c(Max_GeneEXP_cpm,seq(Min_GeneEXP_cpm,Max_GeneEXP_cpm, by=step))))
      labels_cpm<-as.character(round(breaks_cpm,1))
      
      if(DEBUG == 1)
      {
        cat("breaks_cpm:\t")
        cat(sprintf(as.character(breaks_cpm)))
        cat("\n")
        cat("labels_cpm:\t")
        cat(sprintf(as.character(labels_cpm)))
        cat("\n")
      }
      
      breaks_time<-sort(unique(Mean_GeneEXP_df$time))
      labels_time<-as.character(round(breaks_time,0))
      
      if(DEBUG == 1)
      {
        cat("breaks_time:\t")
        cat(sprintf(as.character(breaks_time)))
        cat("\n")
        cat("labels_time:\t")
        cat(sprintf(as.character(labels_time)))
        cat("\n")
      }
      
      order_of_symbol_df<-Genomic_features_sel[which(Genomic_features_sel$gene_name%in%Display_genes_sel |
                                                       Genomic_features_sel$stable_id%in%Display_genes_sel),]
      
      order_of_symbol_df<-order_of_symbol_df[order(order_of_symbol_df$start),]
      
      symbol_terms<-unique(order_of_symbol_df$Symbol)
      
      if(DEBUG == 1)
      {
        cat("order_of_symbol_df_0\n")
        cat(str(order_of_symbol_df))
        cat("\n")
        cat(str(symbol_terms))
        cat("\n")
        
      }
      
      Mean_GeneEXP_df$Symbol<-factor(Mean_GeneEXP_df$Symbol,
                                     levels=symbol_terms,
                                     ordered=T)
      
      if(DEBUG == 1)
      {
        cat("Mean_GeneEXP_df_IMMEDIATELY_PRE\n")
        cat(str(Mean_GeneEXP_df))
        cat("\n")
        
      }
      
      
      
      
      graph_GeneEXP_cpm<-ggplot(data=Mean_GeneEXP_df,
                                aes(x=time,
                                    y=cpm))+   
        new_scale("color")+
        geom_line(aes(group=interaction(Symbol, Genotype), 
                      color=Genotype), 
                  size=1, alpha=0.8)+
        scale_color_manual(values=vector_fill, drop=F)+
        new_scale("fill")+
        new_scale("color")+
        geom_point(aes(fill=Genotype, color=Genotype),
                   size=1.5, shape=21, stroke=0.5, alpha=0.8)+
        scale_color_manual(values=vector_fill, drop=F)+
        scale_fill_manual(values=vector_fill, drop=F)+
        scale_x_continuous(name="Time",breaks=breaks_time,labels=labels_time, limits=c(breaks_time[1],breaks_time[length(breaks_time)]), expand = c(0.1, 0.1))+     
        scale_y_continuous(name="Gene Expression (cpm)",breaks=breaks_cpm,labels=labels_cpm, limits=c(breaks_cpm[1],breaks_cpm[length(breaks_cpm)]), expand = c(0.1, 0.1))
      
      
      graph_GeneEXP_cpm <-graph_GeneEXP_cpm+
        theme_cowplot(font_size = 2,
                      font_family = "sans")+
        facet_grid(Symbol ~ ., scales='free_x', space='free_x', switch="y", drop=F)+
        theme( strip.background = element_blank(),
               strip.placement = "outside",
               strip.text = element_text(size=5,color="black", family="sans"),
               panel.spacing = unit(0.2, "lines"),
               panel.background=element_rect(fill="white"),
               panel.border=element_rect(colour="white",size=0,5),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())+
        theme_classic()+
        theme(axis.title=element_blank(),
              axis.title.y=element_text(size=8,color="black", family="sans"),
              axis.title.x=element_text(size=8,color="black", family="sans"),
              axis.text.y=element_text(size=6,color="black", family="sans"),
              axis.text.x=element_text(size=6,color="black", family="sans"),
              axis.line.x = element_line(size = 0.4),
              axis.ticks.x = element_line(size = 0.4),
              axis.ticks.y = element_line(size = 0.4),
              axis.line.y = element_line(size = 0.4))+
        theme(legend.title = element_text(size=6),
              legend.text = element_text(size=6),
              legend.key.size = unit(0.35, 'cm'), #change legend key size
              legend.key.height = unit(0.35, 'cm'), #change legend key height
              legend.key.width = unit(0.35, 'cm'), #change legend key width
              legend.position="bottom")+
        guides(fill=guide_legend(nrow=2,byrow=TRUE))+
        ggeasy::easy_center_title()
      
      
      
      
      #### ATAC file ------
      
      FLAG_Display_peaks_sel<-sum(Display_peaks_sel == 'NA')
      
      if(DEBUG == 1)
      {
        cat("FLAG_Display_peaks_sel_0\n")
        cat(str(FLAG_Display_peaks_sel))
        cat("\n")
      }
      
      if(FLAG_Display_peaks_sel == 0){
        
        metadata_ATAC_sel<-droplevels(metadata_ATAC[which(metadata_ATAC$seurat_cluster == seurat_clusters_sel),])
        
        if(DEBUG == 1)
        {
          cat("metadata_ATAC_sel_0\n")
          cat(str(metadata_ATAC_sel))
          cat("\n")
        }
        
        if(dim(metadata_ATAC_sel)[1] >0){
          
          metadata_ATAC_sel<-metadata_ATAC_sel[order(metadata_ATAC_sel$Genotype, metadata_ATAC_sel$time_point),]
          
          if(DEBUG == 1)
          {
            cat("metadata_ATAC_sel_ORDERED\n")
            cat(str(metadata_ATAC_sel))
            cat("\n")
          }
          
          
          column_order<-metadata_ATAC_sel$columns_for_matrix
          
          if(DEBUG == 1)
          {
            cat("column_order_0\n")
            cat(sprintf(as.character(column_order)))
            cat("\n")
          }
          
          
          lcpm_ATAC_sel<-lcpm_ATAC[which(row.names(lcpm_ATAC)%in%Display_peaks_sel),
                                   which(colnames(lcpm_ATAC)%in%column_order), drop=F]
          
          if(dim(lcpm_ATAC_sel)[1] > 0){
            
            lcpm_ATAC_sel_df <- as.data.frame(lcpm_ATAC_sel %>%
                                                as_tibble())
            
            lcpm_ATAC_sel_df$Peak_Number<-as.character(row.names(lcpm_ATAC_sel))
            
            if(DEBUG == 1)
            {
              # cat("lcpm_ATAC_sel\n")
              # cat(str(lcpm_ATAC_sel))
              # cat("\n")
              
              #  #  cat("names(lcpm_ATAC_sel)\n")
              #  # cat(sprintf(as.character(names(lcpm_ATAC_sel))))
              #  # cat("\n")
              
              
              #  cat("row.names(lcpm_ATAC_sel)\n")
              #  cat(sprintf(as.character(row.names(lcpm_ATAC_sel))))
              #  cat("\n")
              #    cat("colnames(lcpm_ATAC_sel)\n")
              #  cat(sprintf(as.character(colnames(lcpm_ATAC_sel))))
              #  cat("\n")
              
              #  cat("lcpm_ATAC_sel_df\n")
              # cat(str(lcpm_ATAC_sel_df))
              # cat("\n")
              
            }
            
            lcpm_ATAC_sel_df.m<-melt(lcpm_ATAC_sel_df, id.vars='Peak_Number', variable.name='columns_for_matrix', value.name='cpm')
            
            if(DEBUG == 1)
            {
              cat("lcpm_ATAC_sel_df.m_0\n")
              cat(str(lcpm_ATAC_sel_df.m))
              cat("\n")
            }
            
            lcpm_ATAC_sel_df.m<-merge(lcpm_ATAC_sel_df.m,
                                      metadata_ATAC_sel,
                                      by='columns_for_matrix')
            
            if(DEBUG == 1)
            {
              cat("lcpm_ATAC_sel_df.m_1\n")
              cat(str(lcpm_ATAC_sel_df.m))
              cat("\n")
            }
            
            
            
            
            lcpm_ATAC_sel_df.m.dt<-data.table(lcpm_ATAC_sel_df.m,
                                              key=c('Peak_Number','Genotype','time'))
            
            
            
            Mean_ATAC_df<-as.data.frame(lcpm_ATAC_sel_df.m.dt[,.(cpm=mean(cpm),
                                                                 sd=sd(cpm)), by=key(lcpm_ATAC_sel_df.m.dt)], stringsAsFactors=F)
            
            Mean_ATAC_df$cpm_max<-Mean_ATAC_df$cpm+Mean_ATAC_df$sd
            Mean_ATAC_df$cpm_min<-Mean_ATAC_df$cpm-Mean_ATAC_df$sd
            
            if(DEBUG == 1)
            {
              cat("--------------------------------------------------------------------------->Mean_ATAC_df_0\n")
              cat(str(Mean_ATAC_df))
              cat("\n")
            }
            Master_peak_file_with_SNP_numbered_peaks_sel_ordered<-Master_peak_file_with_SNP_numbered_peaks_sel[order(Master_peak_file_with_SNP_numbered_peaks_sel$start),]
            
            Peak_Number_order<-unique(Master_peak_file_with_SNP_numbered_peaks_sel_ordered$Peak_Number)
            
            
            if(DEBUG == 1)
            {
              cat("Master_peak_file_with_SNP_numbered_peaks_sel_ordered_0\n")
              cat(str(Master_peak_file_with_SNP_numbered_peaks_sel_ordered))
              cat("\n")
              
              cat("Peak_Number_order_0\n")
              cat(str(Peak_Number_order))
              cat("\n")
            }
            
            
            
            colnames(Mean_ATAC_df)[which(colnames(Mean_ATAC_df) == 'Peak_Number')]<-'Feature_ID'
            Mean_ATAC_df$Feature_Type<-'scATAC'
            colnames(Mean_GeneEXP_df)[which(colnames(Mean_GeneEXP_df) == 'Symbol')]<-'Feature_ID'
            Mean_GeneEXP_df$Feature_Type<-'scRNA-seq'
            
            
            
            
            Mean_features_df<-rbind(Mean_GeneEXP_df,Mean_ATAC_df)
            
            
            if(DEBUG == 1)
            {
              cat("---------------KEY RBIND------------------------------------------------------------->Mean_features_df_0\n")             
              cat(str(Mean_features_df))
              cat("\n")           
              
            }
            
            
            
            
            summary_cpm<-unique(c(Mean_features_df$cpm_max[!is.na(Mean_features_df$cpm_max)],Mean_features_df$cpm_min[!is.na(Mean_features_df$cpm_min)]))
            
            Max_GeneEXP_cpm<-max(summary_cpm)
            Min_GeneEXP_cpm<-min(summary_cpm)
            
            step<-round((Max_GeneEXP_cpm-Min_GeneEXP_cpm)/4,2)
            
            if(step == 0)
            {
              step<-0.1
            }
            
            if(DEBUG == 1)
            {
              cat("Max_GeneEXP_cpm:\t")
              cat(sprintf(as.character(Max_GeneEXP_cpm)))
              cat("\n")
              cat("Min_GeneEXP_cpm:\t")
              cat(sprintf(as.character(Min_GeneEXP_cpm)))
              cat("\n")
              cat("step:\t")
              cat(sprintf(as.character(step)))
              cat("\n")
            }
            
            breaks_cpm<-unique(sort(c(Max_GeneEXP_cpm,seq(Min_GeneEXP_cpm,Max_GeneEXP_cpm, by=step))))
            labels_cpm<-as.character(round(breaks_cpm,1))
            
            if(DEBUG == 1)
            {
              cat("breaks_cpm:\t")
              cat(sprintf(as.character(breaks_cpm)))
              cat("\n")
              cat("labels_cpm:\t")
              cat(sprintf(as.character(labels_cpm)))
              cat("\n")
            }
            
            breaks_time<-sort(unique(Mean_features_df$time))
            labels_time<-as.character(round(breaks_time,0))
            
            if(DEBUG == 1)
            {
              cat("breaks_time:\t")
              cat(sprintf(as.character(breaks_time)))
              cat("\n")
              cat("labels_time:\t")
              cat(sprintf(as.character(labels_time)))
              cat("\n")
            } 
            
            
            Mean_features_df$Feature_ID<-factor(Mean_features_df$Feature_ID,
                                                levels=c(symbol_terms,Peak_Number_order),
                                                ordered=T)
            
            Mean_features_df$Feature_Type<-factor(Mean_features_df$Feature_Type,
                                                  levels=c('scRNA-seq','scATAC'),
                                                  ordered=T)
            
            
            if(DEBUG == 1)
            {
              cat(">Mean_features_df_IMMEDIATELY_PRE\n")             
              cat(str(Mean_features_df))
              cat("\n")           
              
            }
            
            pd <- position_dodge(0.2)
            
            
            graph_ALL_cpm<-ggplot(data=Mean_features_df,
                                  aes(x=time,
                                      y=cpm))+   
              new_scale("color")+
              geom_line(aes(group=interaction(Feature_ID, 
                                              Genotype), 
                            color=Genotype), size=0.5, alpha=0.8)+
              geom_errorbar(aes(ymin=cpm-sd, ymax=cpm+sd,
                                group=interaction(Feature_ID, 
                                                  Genotype), 
                                color=Genotype), width=3)+
              scale_color_manual(values=vector_fill, drop=F)+
              new_scale("fill")+
              new_scale("color")+
              geom_point(aes(fill=Genotype, 
                             color=Genotype),
                         size=1.5, shape=21, stroke=0.5,alpha=0.8)+
              scale_color_manual(values=vector_fill, drop=F)+
              scale_fill_manual(values=vector_fill, drop=F)+     
              scale_x_continuous(name="Time",breaks=breaks_time,labels=labels_time, limits=c(breaks_time[1],breaks_time[length(breaks_time)]), expand = c(0.1, 0.1))+     
              scale_y_continuous(name="CPM",breaks=breaks_cpm,labels=labels_cpm, limits=c(breaks_cpm[1],breaks_cpm[length(breaks_cpm)]), expand = c(0.1, 0.1))
            
            
            
            graph_ALL_cpm <-graph_ALL_cpm+
              theme_cowplot(font_size = 2,
                            font_family = "sans")+
              facet_grid(. ~ Feature_Type+Feature_ID, scales='free_x', space='free_x', switch="y", labeller=labeller(paste0(Mean_features_df$Feature_Type, "\n",Mean_features_df$Feature_ID)))+
              theme( strip.background = element_blank(),
                     strip.placement = "outside",
                     strip.text = element_text(size=5,color="black", family="sans"),
                     panel.spacing = unit(0.2, "lines"),
                     panel.background=element_rect(fill="white"),
                     panel.border=element_rect(colour="white",size=0,5),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
              theme_classic()+
              theme(axis.title=element_blank(),
                    axis.title.y=element_text(size=8,color="black", family="sans"),
                    axis.title.x=element_text(size=8,color="black", family="sans"),
                    axis.text.y=element_text(size=6,color="black", family="sans"),
                    axis.text.x=element_text(size=6,color="black", family="sans"),
                    axis.line.x = element_line(size = 0.4),
                    axis.ticks.x = element_line(size = 0.4),
                    axis.ticks.y = element_line(size = 0.4),
                    axis.line.y = element_line(size = 0.4))+
              theme(legend.title = element_text(size=6),
                    legend.text = element_text(size=6),
                    legend.key.size = unit(0.35, 'cm'), #change legend key size
                    legend.key.height = unit(0.35, 'cm'), #change legend key height
                    legend.key.width = unit(0.35, 'cm'), #change legend key width
                    legend.position="bottom")+
              guides(fill=guide_legend(nrow=2,byrow=TRUE))+
              ggeasy::easy_center_title()
            
            
            
            graph_GeneEXP_and_ATAC<-graph_ALL_cpm
            
          }else{
            graph_GeneEXP_and_ATAC<-graph_GeneEXP_cpm
            
          }#dim(lcpm_ATAC_sel)[1] > 0
          
        }else{
          
          graph_GeneEXP_and_ATAC<-graph_GeneEXP_cpm
          
        }#dim(metadata_ATAC_sel)[1] >0
        
      }else{
        graph_GeneEXP_and_ATAC<-graph_GeneEXP_cpm
      }# FLAG_Display_peaks_sel == 0
      
      
      
      
      
      ##### Put together graphs -----
      
      graph_cov<-plot_grid(cov_plot,NA,
                           ncol=2,
                           nrow=1,
                           rel_widths=c(1,0.01))
      
      graph_DEF<-plot_grid(graph_cov,graph_Genomic_features,graph_GeneEXP_and_ATAC,
                           ncol=1,
                           nrow=3,
                           rel_heights=c(1.5,0.5,0.5))
      
      
      setwd(path_graphs)
      
      svgname<-paste(paste("Situation_plot",seurat_clusters_sel,Region_sel,sep='_'),".svg",sep='')
      svglite(svgname, width = 8, height = 13)
      print(graph_DEF)
      dev.off()                
    }# k in 1:length(cluster_levels)          
  }#i in 1:dim(Input_regions_RUSH)[1]
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
    make_option(c("--Master_peak_file_with_SNP_numbered_peaks"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Genomic_features"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Input_regions"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--metadata_GeneEXP"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GeneEXP"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--metadata_ATAC"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ATAC"), type="character", default=NULL, 
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
  