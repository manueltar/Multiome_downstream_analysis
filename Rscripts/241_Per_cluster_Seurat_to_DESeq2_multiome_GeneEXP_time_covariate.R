
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
  
  #### Extract RNA counts ----
  
  matrix_RNA<-GetAssayData(object = Seurat_object, assay = "RNA", layer = "counts")
  
  cat("matrix_RNA\n")
  cat(str(matrix_RNA))
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
  
  
  cat(sprintf(as.character(names(summary(metadata$seurat_clusters)))))
  cat("\n")
  cat(sprintf(as.character(summary(metadata$seurat_clusters))))
  cat("\n")
  
  
  cat(sprintf(as.character(names(summary(as.factor(metadata$Assigned_GFPbc))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(metadata$Assigned_GFPbc)))))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(metadata$sample_id))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(metadata$sample_id)))))
  cat("\n")
  
  
  
  #### create new Seurat object ----
  
  RNA_object <- CreateSeuratObject(counts = matrix_RNA, assay = "RNA",
                                   meta.data=metadata)
  
  # cat("RNA_object_0\n")
  # cat(str(RNA_object))
  # cat("\n")
  
  #### Pseudobulk 1 Aggregate sample_id by seurat_clusters ----
  
  cluster_names <- levels(metadata[,which(colnames(metadata) == 'seurat_clusters')])
  
  cat("cluster_names_0\n")
  cat(str(cluster_names))
  cat("\n")
  
  sample_names <- levels(metadata[,which(colnames(metadata) == 'sample_id')])
  
  cat("sample_names_0\n")
  cat(str(sample_names))
  cat("\n")
  
  groups <- metadata[,c(which(colnames(metadata) == 'sample_id'),which(colnames(metadata) == 'seurat_clusters'))]
  
  cat("groups_0\n")
  cat(str(groups))
  cat("\n")
  
  aggr_counts <- Seurat2PB(RNA_object, sample="sample_id", cluster="seurat_clusters")
  
  cat("aggr_counts_0\n")
  cat(str(aggr_counts))
  cat("\n")
  
  
  #### Aggregate per cluster----
  
  
  # As a reminder, we stored our cell types in a vector called cluster_names
  cluster_names
  
  
  
  # Loop over all cell types to extract corresponding counts, and store information in a list
  
  ## Initiate empty list
  counts_ls <- list()
  
  DEBUG<-1
  
  for (i in 1:length(cluster_names)) {
    
    cluster_names[i]
    
    ## Extract indexes of columns in the global matrix that match a given cluster
    column_idx <- which(tstrsplit(colnames(aggr_counts), "_cluster")[[2]] == cluster_names[i])
    
    sub_aggr<- aggr_counts[, column_idx]
    
    if(DEBUG == 1)
    {
      cat("sub_aggr_0\n")
      cat(str(sub_aggr))
      cat("\n")
      
    }
    
    ## Store corresponding sub-matrix as one element of a list
    counts_ls[[i]] <-sub_aggr
    names(counts_ls)[i] <- cluster_names[i]
    
  }
  
  # Explore the different components of the list
  
  cat("counts_ls_0\n")
  cat(str(counts_ls))
  cat("\n")
  
  #### Create group level variables -----
  
  # Extract sample-level variables
  metadata_NEW <- metadata %>% 
    as.data.frame() %>% 
    dplyr::select(Genotype, Assigned_GFPbc, time_point, sample_id)
  
  
  # Exclude duplicated rows
  metadata_NEW <- metadata_NEW[!duplicated(metadata_NEW), ]
  
 
  # Rename rows
  rownames(metadata_NEW) <- metadata_NEW$sample_id
  
  cat("metadata_NEW_0\n")
  cat(str(metadata_NEW))
  cat("\n")
  
  #### aggregate cell counts per diff group and sample id -----
  
  
  t <- table(metadata$sample_id,
             metadata$seurat_clusters)
  
  cat("t_0\n")
  cat(str(t))
  cat("\n")
  
  ##### Creating metadata list ----------------
  
  ## Initiate empty list
  metadata_ls <- list()
  
  DEBUG<-1
  
  for (i in 1:length(counts_ls)) {
    
    ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
    
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$cluster_id <- tstrsplit(df$cluster_sample_id, "_cluster")[[2]]
    df$sample_id  <- tstrsplit(df$cluster_sample_id, "_cluster")[[1]]
    
    
    
    
    
    ## Retrieve cell count information for this cluster from global cell count table
    idx <- which(colnames(t) == unique(df$cluster_id))
    cell_counts <- t[, idx]
    
    #   cat("cell_counts_0\n")
    # cat(str(cell_counts))
    # cat("\n")
    
    
    
    
    
    ## Remove samples with zero cell contributing to the cluster
    cell_counts <- cell_counts[cell_counts > 0]
    
    #     cat("cell_counts_1\n")
    # cat(str(cell_counts))
    # cat("\n")
    
    
    ## Match order of cell_counts and sample_ids
    sample_order <- match(df$sample_id, names(cell_counts))
    cell_counts <- cell_counts[sample_order]
    
    
    
    ## Append cell_counts to data frame
    df$cell_count <- cell_counts
    
    
    ## Join data frame (capturing metadata specific to cluster) to generic metadata
    df <- plyr::join(df, metadata_NEW, 
                     by = intersect(names(df), names(metadata_NEW)))
    
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df) <- df$cluster_sample_id
    
    ## Store complete metadata for cluster i in list
    metadata_ls[[i]] <- df
    names(metadata_ls)[i] <- unique(df$cluster_id)
    
  }
  
  # Explore the different components of the list
  
  cat("metadata_ls_0\n")
  cat(str(metadata_ls))
  cat("\n")
  
  
  #### DESEq2 analysis -----
  
  
  # # Select cell type of interest
  # cluster_names
  
  # Double-check that both lists have same names
  check_1<-all(names(counts_ls) == names(metadata_ls))
  
  cat("check_1\n")
  cat(sprintf(as.character(check_1)))
  cat("\n")
  
  
  idx <- which(names(counts_ls) == seurat_cluster)
  cluster_counts <- counts_ls[[idx]]
  
  cat("cluster_counts_0\n")
  cat(str(cluster_counts))
  cat("\n")
  
  
  cluster_metadata <- metadata_ls[[idx]]
  cluster_metadata$Genotype<-as.character(cluster_metadata$Genotype)
  cluster_metadata$Genotype<-factor(cluster_metadata$Genotype)
  
  cluster_metadata$Genotype<-relevel(cluster_metadata$Genotype, ref='G/G') ### Nothing works with ordered factors
  
  cluster_metadata$time_point<-as.character(cluster_metadata$time_point)
  
  cluster_metadata$time<-NA
  
  cluster_metadata$time[which(cluster_metadata$time_point == 'h24')]<-24
  cluster_metadata$time[which(cluster_metadata$time_point == 'h48')]<-48
  cluster_metadata$time[which(cluster_metadata$time_point == 'h72')]<-72
  cluster_metadata$time[which(cluster_metadata$time_point == 'h96')]<-96
  
  
  cat("cluster_metadata_0\n")
  cat(str(cluster_metadata))
  cat("\n")
  cat(str(unique(cluster_metadata$cluster_id)))
  cat("\n")
  
 
  
  check_2<-all(colnames(cluster_counts) == rownames(cluster_metadata))
  
  cat("check_2\n")
  cat(sprintf(as.character(check_1)))
  cat("\n")

#### LOOP FOR COMPARISON OF GENOTYPES ------
  
  
  Results_genotype_comparisons<-data.frame()
  
  array_Genotypes<-levels(cluster_metadata$Genotype)
  
  cat("cluster_metadata_0\n")
  cat(str(cluster_metadata))
  cat("\n")
  
  for(i in 2:length(array_Genotypes))
  {
    comparison_genotype<-array_Genotypes[1]
    array_Genotypes_sel<-array_Genotypes[i]
    
    cat("------------------------->\t")
    cat(sprintf(as.character(comparison_genotype)))
    cat("\t")
    cat(sprintf(as.character(array_Genotypes_sel)))
    cat("\n")
    
    
    cluster_metadata_sel<-cluster_metadata[which(cluster_metadata$Genotype%in%c(comparison_genotype,array_Genotypes_sel)),]
    
    
    cluster_metadata_sel$Genotype<-factor(cluster_metadata_sel$Genotype)
    
    cluster_metadata_sel$Genotype<-relevel(cluster_metadata_sel$Genotype, ref='G/G') ### Nothing works with ordered factors
    
    
    
    cluster_metadata_sel$time<-NA
    
    cluster_metadata_sel$time[which(cluster_metadata_sel$time_point == 'h24')]<-24
    cluster_metadata_sel$time[which(cluster_metadata_sel$time_point == 'h48')]<-48
    cluster_metadata_sel$time[which(cluster_metadata_sel$time_point == 'h72')]<-72
    cluster_metadata_sel$time[which(cluster_metadata_sel$time_point == 'h96')]<-96
    
    # str(cluster_metadata_sel)
    # str(unique(cluster_metadata_sel$cluster_id))
    # str(unique(cluster_metadata_sel$Genotype))
    #  str(unique(cluster_metadata_sel$sample_id))
    
    indx.selected_cols<-which(colnames(cluster_counts)%in%rownames(cluster_metadata_sel))
    
    cluster_counts_sel<-cluster_counts[,indx.selected_cols]
    
    # str(cluster_counts_sel)
    
    # Check matching of matrix columns and metadata rows
    all(colnames(cluster_counts_sel) == rownames(cluster_metadata_sel))
    
    ############# Create DESeq2 object -----------------------
    
    dds <- DESeqDataSetFromMatrix(cluster_counts_sel, 
                                  colData = cluster_metadata_sel, 
                                  design = ~ Genotype + time)
    
    ############# Test object -----------------------
    
    dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ time)
    
    names<-resultsNames(dds_lrt)
    
    cat("Names\n")
    cat(sprintf(as.character(names)))
    cat("\n")
    
    coeff_sel<-names[-which(names%in%c('time','Intercept'))]
    
    cat("coeff_sel\n")
    cat(sprintf(as.character(coeff_sel)))
    cat("\n")
    
    res_LRT <- results(dds_lrt,
                       name = coeff_sel,
                       alpha = 0.05)
    
    # cat("res_LRT_0\n")
    # str(res_LRT)
    # cat("\n")
    
    
    
    res_LRT <- lfcShrink(dds_lrt, 
                         coef = coeff_sel,
                         res=res_LRT,
                         type = "apeglm")
    
    #   cat("res_LRT_1\n")
    # str(res_LRT)
    # cat("\n")
    
    
    res_LRT_tb <- res_LRT %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>%
      as_tibble() %>%
      arrange(padj)
    
    colnames(res_LRT_tb)[which(colnames(res_LRT_tb) == 'gene')]<-'Symbol'
    
    # Check results output
    
    res_LRT_tb$comparison<-coeff_sel    #paste('Genotype_',gsub("\\/","\\.",array_Genotypes_sel),'_vs_',gsub("\\/","\\.",comparison_genotype),sep='')
    res_LRT_tb$Minus_logpval<-round(-1*log10(res_LRT_tb$padj),2)
    
    res_LRT_tb<-as.data.frame(res_LRT_tb)
    
    # cat("res_LRT_tb_0\n")
    # str(res_LRT_tb)
    # cat("\n")
    
    
    Results_genotype_comparisons<-rbind(res_LRT_tb,Results_genotype_comparisons)
    
    
    
   
    #### PCA and heatmap ----
    
    # Transform counts for data visualization
    rld <- rlog(dds_lrt, blind=TRUE)
    
    setwd(path_graphs)
    
    svgname<-paste('PCA_Genotype_',coeff_sel,'.svg', sep='')
    
    svglite(svgname, width = 6, height = 6)
    DESeq2::plotPCA(rld, ntop = 500, intgroup = "Genotype")
    dev.off()
    
    setwd(path_graphs)
    
    svgname<-paste('PCA_time_point_',coeff_sel,'.svg', sep='')
    
    svglite(svgname, width = 6, height = 6)
    DESeq2::plotPCA(rld, ntop = 500, intgroup = "time_point")
    dev.off()
    
    setwd(path_graphs)
    
    svgname<-paste('PCA_sample_id_',coeff_sel,'.svg', sep='')
    
    svglite(svgname, width = 10, height = 10)
    DESeq2::plotPCA(rld, ntop = 500, intgroup = "sample_id")
    dev.off()
    
    setwd(path_graphs)
    
    svgname<-paste('PCA_cell_count_',coeff_sel,'.svg', sep='')

    svglite(svgname, width = 6, height = 6)
    DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count")
    dev.off()
    
    # Extract the rlog matrix from the object and compute pairwise correlation values
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)
    
    heatmap<-pheatmap(rld_cor, annotation = cluster_metadata[, c("Genotype"), drop=F])
    
    
    # Plot heatmap
    setwd(path_graphs)
    
    svgname<-paste('Heatmap_correlation_',coeff_sel,'.svg',sep='')
    
    svglite(svgname, width = 10, height = 10)
    print(heatmap)
    dev.off()
    
  
    #### Plot dispersion ----
    
    setwd(path_graphs)
    
    svgname<-'Dispersion.svg'
    
    svglite(svgname, width = 10, height = 10)
    plotDispEsts(dds_lrt)
    dev.off()
    
    #### Retrieve normalised counts to plot later for all genes -----
    
    array_genes<-row.names(dds_lrt@assays@data@listData$counts)
    
    cat("array_genes_0\n")
    cat(str(array_genes))
    cat("\n")
    
    filename<-paste('Normalised_counts_',coeff_sel,'.tsv',sep='')
    setwd(out)
    
    cat(paste(c('count','Genotype','time','sample_id','time_point','cluster_id','Assigned_GFPbc','Symbol'), collapse="\t"),
        file=filename,sep="\n")
    
    ALL_counts<-data.frame()
    
    for(k in 1:length(array_genes))
    {
      array_genes_sel<-array_genes[k]
      
      cat("------------------------------------->\t")
      cat(sprintf(as.character(k)))
      cat("\t")
      cat(sprintf(as.character(array_genes_sel)))
      cat("\n")
      
      d <- plotCounts(dds_lrt, gene=array_genes_sel,
                      intgroup=c("Genotype","time","sample_id","time_point","cluster_id","Assigned_GFPbc"), 
                      returnData=TRUE)
      
      d$Symbol<-array_genes_sel
      
      write.table(d, 
                  file=filename, append = TRUE, sep="\t", 
                  quote=F,col.names = F, row.names = F, eol="\n")
      
      
    }#k in 1:length(array_genes)
    
    
  }#i in 2:length(array_Genotypes)
  
  #### SAVE results ----
  
  if(dim(Results_genotype_comparisons)[1] > 0)
  {
    setwd(out)
    # Write all results to file
    write.table(Results_genotype_comparisons,file="DE_genes.tsv",sep="\t",          
                quote = FALSE, 
                row.names = FALSE)
    
  }#dim(Results_genotype_comparisons)[1] > 0
  

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
    make_option(c("--seurat_cluster"), type="character", default=NULL, 
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
  