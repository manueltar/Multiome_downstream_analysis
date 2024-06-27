
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


opt = NULL

options(warn = 1)

collects_GeneEXP = function(option_list)
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
  
  #### READ and transform Genotype_array ----
  
  Genotype_array = unlist(strsplit(opt$Genotype_array, split=","))
  
  cat("comparison_array_\n")
  cat(sprintf(as.character(Genotype_array)))
  cat("\n")
  
  #### READ and transform clone_line_array ----
  
  clone_line_array = unlist(strsplit(opt$clone_line_array, split=","))
  
  cat("comparison_array_\n")
  cat(sprintf(as.character(clone_line_array)))
  cat("\n")
  
  #### READ and transform time_point_array ----
  
  time_point_array = unlist(strsplit(opt$time_point_array, split=","))
  
  cat("comparison_array_\n")
  cat(sprintf(as.character(time_point_array)))
  cat("\n")

 
  #### READ and transform cluster_array ----
  
  cluster_array = unlist(strsplit(opt$cluster_array, split=","))
  
  cat("cluster_array_\n")
  cat(sprintf(as.character(cluster_array)))
  cat("\n")
  
 
  
  #### LOOP comparisons ----
  
  DEBUG <- 1
  
  list_results<-list()
  

  for(i in 1:length(cluster_array))
  {
    cluster_array_sel<-cluster_array[i]
    
    cat("--------comparison-------------->\t")
    cat(sprintf(as.character(cluster_array_sel)))
    cat("\t")
    

    path_comparisons_seurat_cluster<-paste(out,cluster_array_sel,'/',sep='')
    
    cat(sprintf(as.character(path_comparisons_seurat_cluster)))
    cat("\n")
    
   
    #### Open Normalised_counts file ----
    
    if (file.exists(path_comparisons_seurat_cluster)){
      
      setwd(path_comparisons_seurat_cluster)
      
      comparisons_df<-data.frame()
      
      for(iteration_Genotype_array in 2:length(Genotype_array))
      {
        wt_reference<-Genotype_array[1]
        Genotype_array_sel<-Genotype_array[iteration_Genotype_array]
        
        
        genotype_comparison<-paste('Genotype_',Genotype_array_sel,'_vs_',wt_reference, sep='')
        
        cat("Genotype_comparison>\t")
        cat(sprintf(as.character(genotype_comparison)))
        cat("\n")
        
        filename<-paste("Normalised_counts_",genotype_comparison,'.tsv',sep='')
        
        if (file.exists(filename)){
          
          Normalised_counts<-as.data.frame(fread(file=filename, sep="\t", header=T), stringsAsFactors=F)
          
          if(DEBUG == 1)
          {
            cat("Normalised_counts_0\n")
            cat(str(Normalised_counts))
            cat("\n")
          }
          
          if(dim(Normalised_counts)[1] >0)
          {
            colnames(Normalised_counts)[which(colnames(Normalised_counts) == 'cluster_id')]<-'seurat_cluster'
            colnames(Normalised_counts)[which(colnames(Normalised_counts) == 'Assigned_GFPbc')]<-'clone_line'
            
            Normalised_counts$Genotype<-gsub("\\/","\\.",Normalised_counts$Genotype)
           
            
            comparisons_df<-rbind(Normalised_counts,comparisons_df)
            
          }#dim(Normalised_counts)[1] >0
          
        }#file.exists(filename
        
      }#iteration_Genotype_array in 1:length(Genotype_array)
      
      if(dim(comparisons_df)[1] >0)
      {
        
        # Modelling cluster group and genotype vs wt with time as a continuous covariate implies that G.G counts are estimated in every model
        # The values for each of the four comparisons are not the same for G/G
        # I am going to keep the maximum value
        
        comparisons_df.dt<-data.table(comparisons_df, key = c('Peak_ID','sample_id','seurat_cluster'))
        
        comparisons_df_MAX<-as.data.frame(comparisons_df.dt[,.SD[which.max(count)], by=key(comparisons_df.dt)])
        
        if(DEBUG == 1)
        {
          cat("comparisons_df_MAX_0\n")
          cat(str(comparisons_df_MAX))
          cat("\n")
          cat(sprintf(as.character(names(summary(as.factor((unique(comparisons_df_MAX$Genotype))))))))
          cat("\n")
          cat(sprintf(as.character(summary(as.factor((unique(comparisons_df_MAX$Genotype)))))))
          cat("\n")
        }
        
        list_results[[i]]<-comparisons_df_MAX
        
      }#dim(comparisons_df)[1] >0
    }#file.exists(path_comparisons_seurat_cluster)
  }# i in 1:length(cluster_array)
    
   
 
  
  
  if(length(list_results) >0)
  {
    FINAL_df = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
   
    
    FINAL_df$seurat_cluster<-factor(FINAL_df$seurat_cluster,
                                levels=cluster_array,
                                ordered=T)
    
    FINAL_df$time_point<-factor(FINAL_df$time_point,
                              levels=time_point_array,
                              ordered=T)
    
    FINAL_df$Genotype<-factor(FINAL_df$Genotype,
                                   levels=Genotype_array,
                                   ordered=T)
    
    FINAL_df$clone_line<-factor(FINAL_df$clone_line,
                              levels=clone_line_array,
                              ordered=T)
    
    FINAL_df<-FINAL_df[order(FINAL_df$time_point,FINAL_df$seurat_cluster,FINAL_df$clone_line,FINAL_df$Peak_ID),]
    
    cat("FINAL_df_0\n")
    cat(str(FINAL_df))
    cat("\n")
    cat(str(unique(FINAL_df$Peak_ID)))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df$seurat_cluster)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df$seurat_cluster))))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df$time_point)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df$time_point))))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df$Genotype)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df$Genotype))))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df$clone_line)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df$clone_line))))
    cat("\n")
    
    FINAL_df$columns_for_matrix<-paste(FINAL_df$time_point,FINAL_df$seurat_cluster,FINAL_df$clone_line, sep='__')
    
    cat("FINAL_df_1\n")
    cat(str(FINAL_df))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(FINAL_df$columns_for_matrix))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(FINAL_df$columns_for_matrix)))))
    
  
    
    #### pivot to wider ----
    
    
    
    df_for_matrix<-as.data.frame(pivot_wider(FINAL_df,
                                             id_cols=Peak_ID,
                                             names_from=columns_for_matrix,
                                             values_from=count), stringsAsFactors=F)
    
    cat("df_for_matrix_0\n")
    cat(str(df_for_matrix))
    cat("\n")
    
   
    
    setwd(out)
    
    write.table(df_for_matrix, file='Normalised_count_matrix.tsv', sep="\t", row.names=F, quote=F)
    
    
    metadata<-unique(FINAL_df[,-c(which(colnames(FINAL_df) == 'count'),
                                  which(colnames(FINAL_df) == 'Peak_ID'))])
    
    cat("metadata_0\n")
    cat(str(metadata))
    cat("\n")
    
    setwd(out)
    
    saveRDS(metadata, file='Metadata.rds')
    

  }#length(list_results) >0
    
    
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
    make_option(c("--cluster_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Genotype_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--clone_line_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--time_point_array"), type="character", default=NULL, 
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
  
  collects_GeneEXP(opt)

  
}


###########################################################################

system.time( main() )