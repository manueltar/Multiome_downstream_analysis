
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
  
  DEBUG<-1
 #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform marker_genes ----
  
  marker_genes = unlist(strsplit(opt$marker_genes, split=","))
  
  cat("marker_genes_0\n")
  cat(str(marker_genes))
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
  
  #### Add the metadata column seurat_clusters_by_diff_groups----
  
  order_cluster_by_diff_groups<-c(c('7','4','2'),c('12','11'),c('10','13','8','3'),c('9','5','6','1'))
  
  
  
  Seurat_object$seurat_clusters_by_diff_groups<-factor(Seurat_object$seurat_clusters,
                                                              levels=order_cluster_by_diff_groups,
                                                              ordered=T)
  
  Seurat_object$seurat_clusters<-factor(Seurat_object$seurat_clusters,
                                                       levels=as.character(seq(1,13,by=1)),
                                                       ordered=T)
  cat("seurat_clusters\n")
  cat(sprintf(as.character(names(summary(Seurat_object$seurat_clusters)))))
  cat("\n")
  cat(sprintf(as.character(summary(Seurat_object$seurat_clusters))))
  cat("\n")
  
  cat("seurat_clusters_by_diff_groups\n")
  cat(sprintf(as.character(names(summary(Seurat_object$seurat_clusters_by_diff_groups)))))
  cat("\n")
  cat(sprintf(as.character(summary(Seurat_object$seurat_clusters_by_diff_groups))))
  cat("\n")
  
  ##### Do not exclude 16bp del-----
  
  

  Seurat_object_subset<-Seurat_object

  Seurat_object_subset$Genotype<-droplevels(Seurat_object_subset$Genotype)
  
  # cat("Seurat_object_subset_0\n")
  # cat(str(Seurat_object_subset))
  # cat("\n")
  
  
 
  setwd(path_graphs)
  
  Idents(Seurat_object_subset) = 'seurat_clusters'
  
  
  p2 = DimPlot(object = Seurat_object_subset, reduction ='umap.wnn' ,label = TRUE,
               cols=c(brewer.pal(9, "YlOrRd")[c(7)],
                      brewer.pal(9, "Blues")[c(5)],
                      brewer.pal(9, "RdPu")[c(6)],
                      brewer.pal(9, "Blues")[c(4)],
                      brewer.pal(9, "Greens")[c(6)],
                      brewer.pal(9, "YlOrRd")[c(6)],
                      brewer.pal(9, "Blues")[c(3)],
                      brewer.pal(9, "RdPu")[c(5)],
                      brewer.pal(9, "Greens")[c(5)],
                      brewer.pal(9, "Greens")[c(4)],
                      brewer.pal(9, "Blues")[c(2)],
                      brewer.pal(9, "RdPu")[c(4)],
                      brewer.pal(9, "RdPu")[c(3)]),
               label.size = 3,
               label.color = "black",
               label.box = TRUE,
               repel=TRUE)+
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
          legend.position="hidden")+        
    ggeasy::easy_center_title()
  
  
  setwd(path_graphs)
  
  svgname<-paste(paste("wnn_","clusters", sep='_'),".svg",sep='')
  svglite(svgname, width = 3, height = 3)
  print(p2)
  dev.off()
  
  
  ##### data wrangling for Stacked barplot -----
  
  
  
  
  df<-as.data.frame(cbind(as.character(Seurat_object_subset$Genotype),
                          as.character(Seurat_object_subset$time_point),
                          as.character(Seurat_object_subset$seurat_clusters),
                          as.character(Seurat_object_subset$diff_groups)), stringsAsFactors=T)
  
  colnames(df)<-c('Genotype','time_point','seurat_clusters','diff_groups')
  
  cat("df_0\n")
  cat(str(df))
  cat("\n")
  
  df$Genotype<-factor(df$Genotype,
                      levels=levels(Seurat_object_subset$Genotype), ordered=T)
  df$diff_groups<-factor(df$diff_groups,
                         levels=levels(Seurat_object_subset$diff_groups), ordered=T)
  
  df$seurat_clusters_by_diff_groups<-factor(df$seurat_clusters,
                                            levels=order_cluster_by_diff_groups, ordered=T)
  
  df$seurat_clusters<-factor(df$seurat_clusters,
                                        levels=as.character(seq(1,13,by=1)),
                                        ordered=T)
  
  cat(sprintf(as.character(names(summary(df$Genotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(df$Genotype))))
  cat("\n")
  cat(sprintf(as.character(names(summary(df$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(df$time_point))))
  cat("\n")
  cat(sprintf(as.character(names(summary(df$seurat_clusters)))))
  cat("\n")
  cat(sprintf(as.character(summary(df$seurat_clusters))))
  cat("\n")
  cat(sprintf(as.character(names(summary(df$diff_groups)))))
  cat("\n")
  cat(sprintf(as.character(summary(df$diff_groups))))
  cat("\n")
  cat(sprintf(as.character(names(summary(df$seurat_clusters_by_diff_groups)))))
  cat("\n")
  cat(sprintf(as.character(summary(df$seurat_clusters_by_diff_groups))))
  cat("\n")
 

  
  

  df.dt<-data.table(df, key=c('time_point','Genotype','seurat_clusters'))
  
  
  Freq.table<-as.data.frame(df.dt[,.(Freq=.N), by=key(df.dt)], stringsAsFactors=F)
  
  cat("Freq.table_0\n")
  cat(str(Freq.table))
  cat("\n")
  
  df.dt<-data.table(df, key=c('time_point','Genotype'))
  
  
  Freq.TOTAL<-as.data.frame(df.dt[,.(TOTAL=.N), by=key(df.dt)], stringsAsFactors=F)
  
  cat("Freq.TOTAL_0\n")
  cat(str(Freq.TOTAL))
  cat("\n")
  
  Freq.table<-merge(Freq.table,
                    Freq.TOTAL,
                    by=c('time_point','Genotype'))
  
  cat("Freq.table_1\n")
  cat(str(Freq.table))
  cat("\n")
  
  Freq.table$Perc<-round(100*(Freq.table$Freq/Freq.table$TOTAL),2)
  
  cat("Freq.table_2\n")
  cat(str(Freq.table))
  cat("\n")
  
  df.dt<-data.table(df, key=c('seurat_clusters'))
  
  
  Freq.Cluster<-as.data.frame(df.dt[,.(Freq=.N), by=key(df.dt)], stringsAsFactors=F)
  
  cat("Freq.table_3\n")
  cat(str(Freq.table))
  cat("\n")
  
  df.dt<-data.table(df, key=c('Genotype'))
  
  
  Freq.Genotype<-as.data.frame(df.dt[,.(nGenotype=.N), by=key(df.dt)], stringsAsFactors=F)
  
  cat("Freq.table_4\n")
  cat(str(Freq.table))
  cat("\n")
  
  Freq.table<-merge(Freq.table,
                    Freq.Genotype,
                    by=c('Genotype'))
  
  cat("Freq.table_5\n")
  cat(str(Freq.table))
  cat("\n")
  
 # order_cluster_by_diff_groups<-c(c('7','4','2'),c('12','11'),c('10','13','8','3'),c('6','1','9','5'))
 
  vector_colors_clusters<-c(brewer.pal(9, "Blues")[c(3)],
                            brewer.pal(9, "Blues")[c(4)],brewer.pal(9, "Blues")[c(5)],brewer.pal(9, "RdPu")[c(4)],
                            brewer.pal(9, "Blues")[c(2)],brewer.pal(9, "Greens")[c(4)],brewer.pal(9, "RdPu")[c(3)],
                            brewer.pal(9, "RdPu")[c(5)],brewer.pal(9, "RdPu")[c(6)],brewer.pal(9, "YlOrRd")[c(6)],
                            brewer.pal(9, "YlOrRd")[c(7)],brewer.pal(9, "Greens")[c(5)],brewer.pal(9, "Greens")[c(6)])
  
  vector_colors_clusters<-c(brewer.pal(9, "YlOrRd")[c(7)],
    brewer.pal(9, "Blues")[c(5)],
    brewer.pal(9, "RdPu")[c(6)],
    brewer.pal(9, "Blues")[c(4)],
    brewer.pal(9, "Greens")[c(6)],
    brewer.pal(9, "YlOrRd")[c(6)],
    brewer.pal(9, "Blues")[c(3)],
    brewer.pal(9, "RdPu")[c(5)],
    brewer.pal(9, "Greens")[c(5)],
    brewer.pal(9, "Greens")[c(4)],
    brewer.pal(9, "Blues")[c(2)],
    brewer.pal(9, "RdPu")[c(4)],
    brewer.pal(9, "RdPu")[c(3)])
  
  cat("vector_colors_clusters_0\n")
  cat(str(vector_colors_clusters))
  cat("\n")
  

  #### Stacked Graph ---------------
  
  breaks.Rank<-(seq(0,100,by=25))
  labels.Rank<-as.character(breaks.Rank)
  
  cat("-------------------------------------->\t")
  cat(sprintf(as.character(labels.Rank)))
  cat("\n")
  
  
  stacked_barplot<-Freq.table %>%
    mutate(myaxis = paste0("\n", "n=", TOTAL), drop=F) %>%
    mutate(myaxis=fct_reorder(myaxis,as.numeric(time_point)), drop=F) %>%
    ggplot(aes(x=myaxis, y=Perc, fill=seurat_clusters)) +
    geom_bar(stat="identity",colour='white')+
    scale_y_continuous(name=paste("Percentage of cells",sep=" "),breaks=breaks.Rank,labels=labels.Rank,
                       limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]+1))+
    scale_fill_manual(values=vector_colors_clusters,
                      drop=F,
                      name="Clusters", breaks=Freq.Cluster$seurat_clusters,
                      labels=paste(Freq.Cluster$seurat_clusters,
                                   Freq.Cluster$Freq, sep =' n= '))
  
  
  stacked_barplot<-stacked_barplot+
    theme_cowplot(font_size = 2)+
    facet_grid(. ~ Genotype+nGenotype, scales='free_x', space='free_x', switch="y", labeller=labeller(paste0(Freq.table$Genotype, "\n", "n=", Freq.table$nGenotype)))+   
    scale_x_discrete(name="Seurat clusters", drop=T)+
    theme( strip.background = element_blank(),
           strip.placement = "outside",
           strip.text = element_text(size=6),
           panel.spacing = unit(0.2, "lines"),
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="white",size=0,5),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    theme_classic()+
    theme(axis.title.y=element_text(size=8, color="black", family="sans"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
          axis.text.x=element_text(angle=45,size=6,vjust=1,hjust=1,color="black", family="sans"),
          axis.line.x = element_line(size = 0.4),
          axis.ticks.x = element_line(size = 0.4),
          axis.ticks.y = element_line(size = 0.4),
          axis.line.y = element_line(size = 0.4))+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=6, color="black", family="sans"),
          legend.key.size = unit(0.35, 'cm'), #change legend key size
          legend.key.height = unit(0.35, 'cm'), #change legend key height
          legend.key.width = unit(0.35, 'cm'), #change legend key width
          legend.position="right")+
    ggeasy::easy_center_title()
  
  setwd(path_graphs)
  
  svgname<-paste(paste("Stacked_barplot","by_G", sep='_'),".svg",sep='')
  svglite(svgname, width = 5, height = 3)
  print(stacked_barplot)
  dev.off()
  
  
  ###### Gather data of marker genes ---------
  
  Results_marker_genes<-data.frame()
  
  Idents(Seurat_object_subset) = 'Genotype'
  
  Seurat_object_subset_WT<-subset(x = Seurat_object_subset, idents = c('G/G'))
  
  
  Idents(Seurat_object_subset_WT) = 'seurat_clusters_by_diff_groups'
  
  Dotplot_object <- DotPlot(object = Seurat_object_subset_WT, features = marker_genes, assay="SCT")
  
  # cat("Dotplot_object_0\n")
  # cat(str(Dotplot_object))
  # cat("\n")
  
  df_dotplot<-as.data.frame(Dotplot_object$data)
  
  cat("df_dotplot_0\n")
  cat(str(df_dotplot))
  cat("\n")
  
  
  if(dim(df_dotplot)[1] >0){
    
    colnames(df_dotplot)[which(colnames(df_dotplot) == 'id')]<-'seurat_clusters_by_diff_groups'
    colnames(df_dotplot)[which(colnames(df_dotplot) == 'features.plot')]<-'Symbol'
    
    df_dotplot$Symbol<-as.character(df_dotplot$Symbol)
    df_dotplot$Assay<-'scRNA-seq'
    
    df_dotplot$Symbol<-factor(df_dotplot$Symbol,
                              levels=rev(marker_genes),
                              ordered=T)
      
    
  
    cat("df_dotplot_1\n")
    cat(str(df_dotplot))
    cat("\n")
        
    
    
    Results_marker_genes<-rbind(df_dotplot,Results_marker_genes)
    
    cat("Results_marker_genes_0\n")
    cat(str(Results_marker_genes))
    cat("\n")
    
    
    Results_marker_genes$GENE_CLASS<-NA
    
    Results_marker_genes$GENE_CLASS[which(Results_marker_genes$Symbol%in%c('ITGA2B','GYPA'))]<-'Flow cyt markers'
    Results_marker_genes$GENE_CLASS[which(Results_marker_genes$Symbol%in%c("IL1B","IL33","CCR7"))]<-'Myeloid leukemia markers'
    Results_marker_genes$GENE_CLASS[which(Results_marker_genes$Symbol%in%c("GP6","FYB1"))]<-'Megakaryocyte markers'
    Results_marker_genes$GENE_CLASS[which(Results_marker_genes$Symbol%in%c("HBQ1","HBZ","HBA2"))]<-'Hemoglobin genes'
    Results_marker_genes$GENE_CLASS[which(Results_marker_genes$Symbol%in%c("ANAPC7","STIL","KIF15"))]<-'Meakaryocyte polyploidization'
    
    Results_marker_genes$GENE_CLASS<-factor(Results_marker_genes$GENE_CLASS,
                                            levels=c('Meakaryocyte polyploidization','Hemoglobin genes','Megakaryocyte markers','Myeloid leukemia markers','Flow cyt markers'),
                                            ordered=T)
    
    cat("Results_marker_genes_1\n")
    cat(str(Results_marker_genes))
    cat("\n")
    
    
    
    #### Dotplot ----
    
    
    breaks_pct.exp<-unique(sort(unique(c(0,seq(0,100, by=25)))))
    labels_pct.exp<-as.character(round(breaks_pct.exp,2))
    
    if(DEBUG == 1)
    {      
      cat("labels_pct.exp\n")
      cat(sprintf(as.character(labels_pct.exp)))
      cat("\n")
    }
    
    ### graph parameters_avg.exp ----
    
    indx_avg.exp<-which(colnames(df_dotplot) == 'avg.exp.scaled')
    
    A_avg.exp<-summary(df_dotplot[,indx_avg.exp])
    
    
    if(DEBUG == 1)
    {
      cat("A_avg.exp\n")
      cat(sprintf(as.character(names(A_avg.exp))))
      cat("\n")
      cat(sprintf(as.character(A_avg.exp)))
      cat("\n")
    }
    
    
    
    max_abs_value<-abs(A_avg.exp[6])
    min_abs_value<-abs(A_avg.exp[1])
    
    if(max_abs_value > min_abs_value)
    {
      step<-round(abs(max_abs_value--1*max_abs_value)/4,2)
      
      breaks_avg.exp<-unique(sort(round(c(0,max_abs_value,seq(-1*max_abs_value,max_abs_value, by=step)),1)))
      labels_avg.exp<-as.character(breaks_avg.exp)
      
    }else{
      
      step<-round(abs(min_abs_value--1*min_abs_value)/4,2)
      
      breaks_avg.exp<-unique(sort(round(c(0,min_abs_value,seq(-1*min_abs_value,min_abs_value, by=step)),1)))
      labels_avg.exp<-as.character(breaks_avg.exp)
      
    }# max_abs_value > min_abs_value
    
    
    
    if(DEBUG == 1)
    {
      cat("labels_avg.exp\n")
      cat(sprintf(as.character(labels_avg.exp)))
      cat("\n")
    }
    
    
    # order_cluster_by_diff_groups<-c(c('7','4','2'),c('12','11'),c('10','13','8','3'),c('6','1','9','5'))
    
    
    Dotplot<-ggplot()+
      geom_point(data=df_dotplot,
                 aes(y=Symbol,
                     x=seurat_clusters_by_diff_groups,                       
                     fill=avg.exp.scaled,
                     size=pct.exp),
                 stroke=0.25, shape=21,color='black')+
      scale_size(range = c(0,6), name='% of cells positive',
                 breaks=breaks_pct.exp, labels=labels_pct.exp, limits=c(breaks_pct.exp[1],breaks_pct.exp[length(breaks_pct.exp)]))+
      scale_fill_gradient2(
        low = "gray", 
        mid = "white", 
        high = "chartreuse3", 
        midpoint = -0.5,
        breaks=breaks_avg.exp,labels=labels_avg.exp,
        limits=c(breaks_avg.exp[1]-0.01,breaks_avg.exp[length(breaks_avg.exp)]+0.01),
        name=paste('Av. Exp.','scaled',sep="\n"),na.value = NA)+       
      scale_x_discrete(name=NULL,drop=F)+
      scale_y_discrete(name=NULL,drop=F)+
      geom_vline(xintercept=c(3.5,5.5,9.5), color="black", linetype='dashed',size=0.2)+
      theme_classic()+
      theme(axis.title=element_blank(),
            axis.title.y=element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_text(size=6, color="black", family="sans", face ="italic"),
            axis.text.x=element_text(angle=0,size=6,color="black", family="sans", face ="plain"),
            axis.line.x = element_line(size = 0.3),
            axis.ticks.x = element_line(size = 0.3),
            axis.ticks.y = element_line(size = 0.3),
            axis.line.y = element_line(size = 0.3))+
      theme(legend.title = element_text(size=6),
            legend.text = element_text(size=6),
            legend.key.size = unit(0.25, 'cm'), #change legend key size
            legend.key.height = unit(0.25, 'cm'), #change legend key height
            legend.key.width = unit(0.25, 'cm'), #change legend key width
            legend.position="bottom")+
      ggeasy::easy_center_title()
    
    
    # Dotplot <-Dotplot+
    #   geom_text_repel(data=Results_marker_genes,
    #                   aes(y=Symbol,
    #                       x=seurat_clusters_by_diff_groups,
    #                       label=Symbol,
    #                       color=GENE_CLASS),
    #                   family="sans",
    #                   fontface='italic',
    #                   segment.size  = 0.2,
    #                   segment.color = "black",
    #                   force=25,
    #                   size=2,
    #                   box.padding = 1,
    #                   max.overlaps = Inf)
    
    # guides(fill=guide_legend(nrow=2,byrow=TRUE))+
      
    setwd(path_graphs)
    
    svgname<-paste(paste("Dotplot","marker_genes", sep='_'),".svg",sep='')
    svglite(svgname, width = 3, height = 3)
    print(Dotplot)
    dev.off()
    
  }#dim(df_dotplot)[1] >0)         
  
  
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
    make_option(c("--marker_genes"), type="character", default=NULL, 
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
  