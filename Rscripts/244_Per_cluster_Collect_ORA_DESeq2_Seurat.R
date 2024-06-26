
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

collects_ORA = function(option_list)
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
  
  cat("cluster_array_\n")
  cat(sprintf(as.character(cluster_array)))
  cat("\n")
  
  #### LOOP comparisons ----
  
  DEBUG <- 0
  
  list_results<-list()
  

  for(i in 1:length(cluster_array))
  {
    cluster_array_sel<-cluster_array[i]
    
    cat("--------comparison-------------->\t")
    cat(sprintf(as.character(cluster_array_sel)))
    cat("\t")
    

    path_comparisons_cluster<-paste(out,cluster_array_sel,'/','MySigDb','/',sep='')
    
    cat(sprintf(as.character(path_comparisons_cluster)))
    cat("\n")
    
    if (file.exists(path_comparisons_cluster)){
      
      
    }else{
      
      
      
    }#path_comparisons_cluster
    
    #### Open MySigDB_results file ----
    
    if (file.exists(path_comparisons_cluster)){
      setwd(path_comparisons_cluster)
      
      filename<-paste("results",'.tsv',sep='')
      
      if (file.exists(filename)){
        
        MySigDB_results<-as.data.frame(fread(file=filename, sep="\t", header=T), stringsAsFactors=F)
        
        if(DEBUG == 1)
        {
          cat("MySigDB_results_0\n")
          cat(str(MySigDB_results))
          cat("\n")
        }
        
        if(dim(MySigDB_results)[1] >0)
        {
          MySigDB_results_subset<-unique(MySigDB_results[,-which(colnames(MySigDB_results) == 'collection')])
          
          if(DEBUG == 1)
          {
            cat("MySigDB_results_subset_0\n")
            cat(str(MySigDB_results_subset))
            cat("\n")
            cat(str(unique(MySigDB_results_subset$id)))
            cat("\n")
          }
          
          
          
          array_id<-unique(MySigDB_results_subset$id)
          
          if(DEBUG == 1)
          {
            cat("array_id_0\n")
            cat(str(array_id))
            cat("\n")
            
          }
          
          
          result_per_id<-data.frame()
          

          for(k in 1:length(array_id))
          {
            array_id_sel<-array_id[k]
            
            
            cat("--------array_id_sel-------------->\t")
            cat(sprintf(as.character(array_id_sel)))
            cat("\n")
            
            MySigDB_results_subset_sel<-MySigDB_results_subset[which(MySigDB_results_subset$id == array_id_sel),]
            
            if(DEBUG == 1)
            {
              cat("MySigDB_results_subset_0\n")
              cat(str(MySigDB_results_subset_sel))
              cat("\n")
              cat(str(unique(MySigDB_results_subset_sel$id)))
              cat("\n")
            }
            
            
            
            
            result_per_id<-rbind(MySigDB_results_subset_sel,result_per_id)
            # quit(status = 1)
            
          }#k in 1:length(array_id)
          

          result_per_id$seurat_cluster<-cluster_array_sel
          
          if(DEBUG == 1)
          {
            cat("result_per_id_0\n")
            cat(str(result_per_id))
            cat("\n")
            cat(str(unique(result_per_id$id)))
            cat("\n")
          }
          
          list_results[[i]]<-result_per_id
          
        }#dim(MySigDB_results)[1] >0
        
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
    
    
    
    indx.dep<-grep(paste(c("NEURON","HEPATOCYTE"), collapse="|"),FINAL_df$id)

    cat("indx.dep_0\n")
    cat(str(indx.dep))
    cat("\n")
    
    FINAL_df_keep<-droplevels(FINAL_df[-indx.dep,])
      
      #
    
    cat("FINAL_df_keep_0\n")
    cat(str(FINAL_df_keep))
    cat("\n")
    cat(str(unique(FINAL_df_keep$id)))
    cat("\n")
    cat(sprintf(as.character(names(summary(FINAL_df_keep$comparison)))))
    cat("\n")
    cat(sprintf(as.character(summary(FINAL_df_keep$comparison))))
    cat("\n")
    
    
    #### Open DE file ----
    
    setwd(out)
    
    write.table(FINAL_df_keep,file="ORA_ActivePathways_results.tsv", sep="\t", quote=F, row.names = F)
    saveRDS(FINAL_df_keep,file="ORA_ActivePathways_results.rds")
    
    
    
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
  
  collects_ORA(opt)

  
}


###########################################################################

system.time( main() )