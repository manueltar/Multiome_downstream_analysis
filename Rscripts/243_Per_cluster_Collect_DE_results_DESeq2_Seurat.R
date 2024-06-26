
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

collects_DE = function(option_list)
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
  
  #### READ and transform comparison_array ----
  
  comparison_array = unlist(strsplit(opt$comparison_array, split=","))
  
  cat("comparison_array_\n")
  cat(sprintf(as.character(comparison_array)))
  cat("\n")

 
  #### READ and transform cluster_array ----
  
  cluster_array = unlist(strsplit(opt$cluster_array, split=","))
  
  cat("seurat_cluster_array_\n")
  cat(sprintf(as.character(cluster_array)))
  cat("\n")
  
  #### LOOP comparisons ----
  
  DEBUG <- 1
  
  list_results<-list()
  

  for(i in 1:length(cluster_array))
  {
    seurat_cluster_array_sel<-cluster_array[i]
    
    cat("--------comparison-------------->\t")
    cat(sprintf(as.character(seurat_cluster_array_sel)))
    cat("\t")
    

    path_comparisons_seurat_cluster<-paste(out,seurat_cluster_array_sel,'/',sep='')
    
    cat(sprintf(as.character(path_comparisons_seurat_cluster)))
    cat("\n")
    
    if (file.exists(path_comparisons_seurat_cluster)){
      
      
    }else{
      
      
      
    }#path_comparisons_seurat_cluster
    
    #### Open DE_results file ----
    
    if (file.exists(path_comparisons_seurat_cluster)){
      setwd(path_comparisons_seurat_cluster)
      
      filename<-paste("DE_genes",'.tsv',sep='')
      
      if (file.exists(filename)){
        
        DE_results<-as.data.frame(fread(file=filename, sep="\t", header=T), stringsAsFactors=F)
        
        if(DEBUG == 1)
        {
          cat("DE_results_0\n")
          cat(str(DE_results))
          cat("\n")
        }
        
        if(dim(DE_results)[1] >0)
        {
          

          DE_results$seurat_cluster<-seurat_cluster_array_sel
          
          if(DEBUG == 1)
          {
            cat("DE_results_0\n")
            cat(str(DE_results))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor((unique(DE_results$comparison))))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor((unique(DE_results$comparison)))))))
            cat("\n")
          }
          
          list_results[[i]]<-DE_results
          
        }#dim(DE_results)[1] >0
        
      }#file.exists(filename
    }
  }# i in 1:length(cluster_array)
    
   
 
  
  
  if(length(list_results) >0)
  {
    FINAL_df = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
   
    
    FINAL_df$seurat_cluster<-factor(FINAL_df$seurat_cluster,
                                levels=cluster_array,
                                ordered=T)
    
    FINAL_df$comparison<-factor(FINAL_df$comparison,
                                   levels=comparison_array,
                                   ordered=T)
    
    cat("FINAL_df_0\n")
    cat(str(FINAL_df))
    cat("\n")
    cat(str(unique(FINAL_df$id)))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df$comparison)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df$comparison))))
    cat("\n")
    
     
    #### Perform a DE analysis summary-----
    
    FINAL_df.dt<-data.table(FINAL_df, key=c("seurat_cluster","comparison"))
    
    
    
    Summary_table<-as.data.frame(FINAL_df.dt[,.(n_genes=.N), by=key(FINAL_df.dt)])
    
    
    cat("Summary_table_0\n")
    cat(str(Summary_table))
    cat("\n")
   
    Summary_table_DE <- as.data.frame(FINAL_df.dt[Minus_logpval >=1.3,.(n_DE_genes=.N), by = key(FINAL_df.dt)])
    
    
    
    cat("Summary_table_DE_0\n")
    cat(str(Summary_table_DE))
    cat("\n")
    
    Summary_table_DE_UP <- as.data.frame(FINAL_df.dt[Minus_logpval >=1.3 & log2FoldChange >=0.25,.(n_DE_genes_UP=.N), by = key(FINAL_df.dt)])
    
    
    
    cat("Summary_table_DE_0\n")
    cat(str(Summary_table_DE_UP))
    cat("\n")
    
    
    Summary_table_DE_DOWN <- as.data.frame(FINAL_df.dt[Minus_logpval >=1.3 & log2FoldChange <= -0.25,.(n_DE_genes_DOWN=.N), by = key(FINAL_df.dt)])
    
    
    
    cat("Summary_table_DE_0\n")
    cat(str(Summary_table_DE_DOWN))
    cat("\n")
    
    Summary_table<-merge(Summary_table,
                         Summary_table_DE,
                         by=c("seurat_cluster","comparison"), all.x=TRUE)
    
    cat("Summary_table_1\n")
    cat(str(Summary_table))
    cat("\n")
    
    Summary_table<-merge(Summary_table,
                         Summary_table_DE_UP,
                         by=c("seurat_cluster","comparison"), all.x=TRUE)
    
    cat("Summary_table_2\n")
    cat(str(Summary_table))
    cat("\n")
    
    Summary_table<-merge(Summary_table,
                         Summary_table_DE_DOWN,
                         by=c("seurat_cluster","comparison"), all.x=TRUE)
    
    cat("Summary_table_3\n")
    cat(str(Summary_table))
    cat("\n")
    
    
    Summary_table[is.na(Summary_table)]<-0
    
    cat("Summary_table_4\n")
    cat(str(Summary_table))
    cat("\n")
    
    #### Open DE file ----
    
    setwd(out)
    
    write.table(FINAL_df,file="DE_results.tsv", sep="\t", quote=F, row.names = F)
    saveRDS(FINAL_df,file="DE_results.rds")
    
    write.table(Summary_table,file="Summary_table_DE.tsv", sep="\t", quote=F, row.names = F)
    
    
  }else{
    
    cat("No significant results\n")
    
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
    make_option(c("--comparison_array"), type="character", default=NULL, 
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
  
  collects_DE(opt)

  
}


###########################################################################

system.time( main() )