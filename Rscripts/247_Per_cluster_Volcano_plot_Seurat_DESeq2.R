
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
suppressMessages(library("stringi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

suppressMessages(library("viridis", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggrepel", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

volcano_function = function(option_list)
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
  
 
  
  #### READ and transform platelet_volume_genes ----
  
  platelet_volume_genes = unique(unlist(strsplit(opt$platelet_volume_genes, split=',')))
  
  cat("platelet_volume_genes_\n")
  cat(sprintf(as.character(platelet_volume_genes)))
  cat("\n")
  #### READ and transform platelet_genes ----
  
  platelet_genes = unique(unlist(strsplit(opt$platelet_genes, split=',')))
  
  cat("platelet_genes_\n")
  cat(sprintf(as.character(platelet_genes)))
  cat("\n")
  #### READ and transform EZH2_signature ----
  
  EZH2_signature = unique(unlist(strsplit(opt$EZH2_signature, split=',')))
  
  cat("EZH2_signature_\n")
  cat(sprintf(as.character(EZH2_signature)))
  cat("\n")
  #### READ and transform ZEB1_signature ----
  
  ZEB1_signature = unique(unlist(strsplit(opt$ZEB1_signature, split=',')))
  
  cat("ZEB1_signature_\n")
  cat(sprintf(as.character(ZEB1_signature)))
  cat("\n")
  #### READ and transform AKT_signature ----
  
  AKT_signature = unique(unlist(strsplit(opt$AKT_signature, split=',')))
  
  cat("AKT_signature_\n")
  cat(sprintf(as.character(AKT_signature)))
  cat("\n")
  
  
  Gather<- data.frame(matrix(vector(), length(c(platelet_volume_genes,platelet_genes,EZH2_signature,ZEB1_signature,AKT_signature)), 3,
                             dimnames=list(c(),
                                           c("Symbol","Path","color"))),stringsAsFactors=F)
  
  Gather$Symbol<-c(platelet_volume_genes,platelet_genes,EZH2_signature,ZEB1_signature,AKT_signature)
  
  Gather$Path<-c(rep('Platelet volume',length(platelet_volume_genes)),
                 rep('Platelet genes',length(platelet_genes)),
                 rep('EZH2 platelet signature',length(EZH2_signature)),
                 rep('ZEB1 signature',length(ZEB1_signature)),
                 rep('AKT signature',length(AKT_signature)))
  
  Gather$color[which(Gather$Path == 'Platelet volume')]<-brewer.pal(8, "Dark2")[8]
  Gather$color[which(Gather$Path == 'Platelet genes')]<-brewer.pal(8, "Dark2")[4]
  Gather$color[which(Gather$Path == 'EZH2 platelet signature')]<-brewer.pal(8, "Accent")[3]
  Gather$color[which(Gather$Path == 'ZEB1 signature')]<-brewer.pal(8, "Accent")[2]
  Gather$color[which(Gather$Path == 'AKT signature')]<-brewer.pal(8, "Dark2")[1]
  
  cat("Gather_0\n")
  cat(str(Gather))
  cat("\n") 
  
  
 
  
  #### Read the DESeq2 results file ----
  
  DE_results<-readRDS(file=opt$DE_results)
  
  cat("DE_results_0\n")
  cat(str(DE_results))
  cat("\n") 
 
  DE_results$Significance<-NA

  DE_results$Significance[which(DE_results$Minus_logpval < 1.3)]<-'NO'
  DE_results$Significance[which(DE_results$Minus_logpval >= 1.3)]<-'YES'

  DE_results$Significance<-factor(DE_results$Significance,
                                             levels=c('NO','YES'),
                                             ordered=T)

  cat("DE_results_1\n")
  cat(str(DE_results))
  cat("\n")
  
  
  DE_results<-merge(DE_results,
        Gather,
        by='Symbol',
        all.x=T)
  
  cat("DE_results_2\n")
  cat(str(DE_results))
  cat("\n")
  
 

  DEBUG <- 1

 

 

  path_graphs<-paste(out,'volcano_plots','/',sep='')

  if (file.exists(path_graphs)){


  }else{

    dir.create(file.path(path_graphs))

  }#path_graphs



  indx_log2FoldChange<-which(colnames(DE_results) == 'log2FoldChange')
  
  A_log2FoldChange<-summary(DE_results[,indx_log2FoldChange])
  
  cat("A_log2FoldChange\n")
  cat(sprintf(as.character(names(A_log2FoldChange))))
  cat("\n")
  cat(sprintf(as.character(A_log2FoldChange)))
  cat("\n")
  
  max_abs_value<-abs(A_log2FoldChange[6])
  min_abs_value<-abs(A_log2FoldChange[1])
  
  if(max_abs_value > min_abs_value)
  {
    step<-round(abs(max_abs_value--1*max_abs_value)/4,2)
    
    breaks_log2FoldChange<-unique(sort(round(c(0,max_abs_value,seq(-1*max_abs_value,max_abs_value, by=step)),1)))
    labels_log2FoldChange<-as.character(breaks_log2FoldChange)
    
  }else{
    
    step<-round(abs(min_abs_value--1*min_abs_value)/4,2)
    
    breaks_log2FoldChange<-unique(sort(round(c(0,min_abs_value,seq(-1*min_abs_value,min_abs_value, by=step)),1)))
    labels_log2FoldChange<-as.character(breaks_log2FoldChange)
    
  }# max_abs_value > min_abs_value
  
  
  cat("step_log2FoldChange\n")
  cat(sprintf(as.character(step)))
  cat("\n")
  cat("labels_log2FoldChange\n")
  cat(sprintf(as.character(labels_log2FoldChange)))
  cat("\n")
  
  
  
  indx_Minus_logpval<-which(colnames(DE_results) == 'Minus_logpval')
  
  A_Minus_logpval<-summary(DE_results[,indx_Minus_logpval])
  
  
  
  cat("A_Minus_logpval\n")
  cat(sprintf(as.character(names(A_Minus_logpval))))
  cat("\n")
  cat(sprintf(as.character(A_Minus_logpval)))
  cat("\n")
  
  
  max_value<-A_Minus_logpval[6]
  min_value<-A_Minus_logpval[1]
  
  
  step<-round(abs(max_value-min_value)/3,0)
  
  if(step == 0)
  {
    
    step<-1
  }
  breaks_Minus_logpval<-unique(sort(unique(c(0,max_value,seq(min_value,max_value, by=step)))))
  labels_Minus_logpval<-as.character(round(breaks_Minus_logpval,1))
  
  
  cat("step_Minus_logpval\n")
  cat(sprintf(as.character(step)))
  cat("\n")
  cat("labels_Minus_logpval\n")
  cat(sprintf(as.character(labels_Minus_logpval)))
  cat("\n")
  
  
  
  if(DEBUG == 1)
  {
    cat("Volcano_START:\n")
    
  }
  
  
  
  
  volcano_plot<-ggplot(data=DE_results,
                       aes(x=log2FoldChange,
                           y=Minus_logpval)) +
    geom_vline(xintercept=c(-0.25,0.25), color="gray", linetype='dashed',size=0.2)+
    geom_point(data=DE_results[which(DE_results$Significance == 'NO'),],
               color="black",fill="gray", stroke=0.2, shape=21, size=0.75)+
    geom_point(data=DE_results[which(DE_results$Significance == 'YES' &
                                           DE_results$log2FoldChange >=0.25),],
               color="black",fill="red", stroke=0.2, shape=21, size=0.75)+
    geom_point(data=DE_results[which(DE_results$Significance == 'YES' &
                                           DE_results$log2FoldChange <= -0.25),],
               color="black",fill="blue", stroke=0.2, shape=21, size=0.75)+
    scale_y_continuous(name='-log10pval',
                       breaks=breaks_Minus_logpval,
                       labels=labels_Minus_logpval,
                       limits=c(breaks_Minus_logpval[1],breaks_Minus_logpval[length(breaks_Minus_logpval)]))+
    scale_x_continuous(name='log2FoldChange',
                       breaks=breaks_log2FoldChange,
                       labels=labels_log2FoldChange,
                       limits=c(breaks_log2FoldChange[1],breaks_log2FoldChange[length(breaks_log2FoldChange)]))
  
  volcano_plot <-volcano_plot+
    theme_cowplot(font_size = 2,
                  font_family = "sans")+
    facet_grid(comparison ~ seurat_cluster, scales='free_x', space='free_x', switch="y", drop=F)+
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
          legend.position="hidden")+
    ggeasy::easy_center_title()
  
  
  
  volcano_plot <-volcano_plot+
    geom_label_repel(data=DE_results[which(DE_results$Symbol%in%Gather$Symbol &
                                             DE_results$Significance == 'YES' &
                                             DE_results$log2FoldChange >=0.25),],
                     aes(x=log2FoldChange,
                         y=Minus_logpval,
                         label=Symbol,
                         color=Path),
                     family="sans",
                     fontface='italic',
                     segment.size  = 0.3,
                     segment.color = "red",
                     force=25,
                     size=2,
                     xlim=c(0.25,breaks_log2FoldChange[length(breaks_log2FoldChange)]),
                     box.padding = 1,
                     max.overlaps = Inf)+
    geom_label_repel(data=DE_results[which(DE_results$Symbol%in%Gather$Symbol &
                                             DE_results$Significance == 'YES' &
                                             DE_results$log2FoldChange <= -0.25),],
                     aes(x=log2FoldChange,
                         y=Minus_logpval,
                         label=Symbol,
                         color=Path),
                     family="sans",
                     fontface='italic',
                     force=25,
                     segment.size  = 0.3,
                     segment.color = "blue",
                     xlim=c(breaks_log2FoldChange[1],-0.25),
                     size=2,
                     box.padding = 1,
                     max.overlaps = Inf)
  
  
  cat("END part volcano_plot\n")
  
  
  setwd(path_graphs)
  
  svgname<-paste(paste("Volcano","ALL_comparisons", sep='_'),".svg",sep='')
  
  setwd(path_graphs)
  svglite(svgname, width = 13, height = 13)
  print(volcano_plot)
  dev.off()
  
 
}

volcano_subset_function = function(option_list)
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
  
  
  
  #### READ and transform platelet_volume_genes ----
  
  platelet_volume_genes = unique(unlist(strsplit(opt$platelet_volume_genes, split=',')))
  
  cat("platelet_volume_genes_\n")
  cat(sprintf(as.character(platelet_volume_genes)))
  cat("\n")
  #### READ and transform platelet_genes ----
  
  platelet_genes = unique(unlist(strsplit(opt$platelet_genes, split=',')))
  
  cat("platelet_genes_\n")
  cat(sprintf(as.character(platelet_genes)))
  cat("\n")
  #### READ and transform EZH2_signature ----
  
  EZH2_signature = unique(unlist(strsplit(opt$EZH2_signature, split=',')))
  
  cat("EZH2_signature_\n")
  cat(sprintf(as.character(EZH2_signature)))
  cat("\n")
  #### READ and transform ZEB1_signature ----
  
  ZEB1_signature = unique(unlist(strsplit(opt$ZEB1_signature, split=',')))
  
  cat("ZEB1_signature_\n")
  cat(sprintf(as.character(ZEB1_signature)))
  cat("\n")
  #### READ and transform AKT_signature ----
  
  AKT_signature = unique(unlist(strsplit(opt$AKT_signature, split=',')))
  
  cat("AKT_signature_\n")
  cat(sprintf(as.character(AKT_signature)))
  cat("\n")
  
  
  Gather<- data.frame(matrix(vector(), length(c(platelet_volume_genes,platelet_genes,EZH2_signature,ZEB1_signature,AKT_signature)), 3,
                             dimnames=list(c(),
                                           c("Symbol","Path","color"))),stringsAsFactors=F)
  
  Gather$Symbol<-c(platelet_volume_genes,platelet_genes,EZH2_signature,ZEB1_signature,AKT_signature)
  
  Gather$Path<-c(rep('Platelet volume',length(platelet_volume_genes)),
                 rep('Platelet genes',length(platelet_genes)),
                 rep('EZH2 platelet signature',length(EZH2_signature)),
                 rep('ZEB1 signature',length(ZEB1_signature)),
                 rep('AKT signature',length(AKT_signature)))
  
  Gather$color[which(Gather$Path == 'Platelet volume')]<-brewer.pal(8, "Dark2")[8]
  Gather$color[which(Gather$Path == 'Platelet genes')]<-brewer.pal(8, "Dark2")[4]
  Gather$color[which(Gather$Path == 'EZH2 platelet signature')]<-brewer.pal(8, "Accent")[3]
  Gather$color[which(Gather$Path == 'ZEB1 signature')]<-brewer.pal(8, "Accent")[2]
  Gather$color[which(Gather$Path == 'AKT signature')]<-brewer.pal(8, "Dark2")[1]
  
  cat("Gather_0\n")
  cat(str(Gather))
  cat("\n") 
  
  #### READ and transform selected_clusters ----
  
  selected_clusters = unique(unlist(strsplit(opt$selected_clusters, split='_')))
  
  cat("selected_clusters_\n")
  cat(sprintf(as.character(selected_clusters)))
  cat("\n")
  
  #### READ and transform selected_comparisons ----
  
  selected_comparisons = unique(unlist(strsplit(opt$selected_comparisons, split=',')))
  
  cat("selected_comparisons_0\n")
  cat(sprintf(as.character(selected_comparisons)))
  cat("\n")
  
  selected_comparisons<-gsub('G_G','G.G', selected_comparisons)
  selected_comparisons<-gsub('A_G','A.G', selected_comparisons)
  selected_comparisons<-gsub('A_A','A.A', selected_comparisons)
  
  cat("selected_comparisons_1\n")
  cat(sprintf(as.character(selected_comparisons)))
  cat("\n")
  
  #### Read the DESeq2 results file ----
  
  DE_results<-readRDS(file=opt$DE_results)
  
  cat("DE_results_0\n")
  cat(str(DE_results))
  cat("\n") 
  
  DE_results$Significance<-NA
  
  DE_results$Significance[which(DE_results$Minus_logpval < 1.3)]<-'NO'
  DE_results$Significance[which(DE_results$Minus_logpval >= 1.3)]<-'YES'
  
  DE_results$Significance<-factor(DE_results$Significance,
                                  levels=c('NO','YES'),
                                  ordered=T)
  
  cat("DE_results_1\n")
  cat(str(DE_results))
  cat("\n")
  
  
  
  DEBUG <- 1
  
  
  DE_results_subset<-droplevels(DE_results[which(DE_results$seurat_cluster%in%selected_clusters &
                                                   DE_results$comparison%in%selected_comparisons),])
  
  cat("DE_results_subset_0\n")
  cat(str(DE_results_subset))
  cat("\n") 
  
  DE_results_subset$seurat_cluster<-factor(DE_results_subset$seurat_cluster,
                                           levels=selected_clusters,
                                           ordered = T)
  
  DE_results_subset$comparison<-factor(DE_results_subset$comparison,
                                           levels=selected_comparisons,
                                           ordered = T)
  
  cat("DE_results_subset_1\n")
  cat(str(DE_results_subset))
  cat("\n")
  
  DE_results_subset<-merge(DE_results_subset,
                           Gather,
                           all.x=T)
  
  cat("DE_results_subset_2\n")
  cat(str(DE_results_subset))
  cat("\n")
  
  
  path_graphs<-paste(out,'volcano_plots','/',sep='')
  
  if (file.exists(path_graphs)){
    
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  
  
  indx_log2FoldChange<-which(colnames(DE_results_subset) == 'log2FoldChange')
  
  A_log2FoldChange<-summary(DE_results_subset[,indx_log2FoldChange])
  
  cat("A_log2FoldChange\n")
  cat(sprintf(as.character(names(A_log2FoldChange))))
  cat("\n")
  cat(sprintf(as.character(A_log2FoldChange)))
  cat("\n")
  
  max_abs_value<-abs(A_log2FoldChange[6])
  min_abs_value<-abs(A_log2FoldChange[1])
  
  if(max_abs_value > min_abs_value)
  {
    step<-round(abs(max_abs_value--1*max_abs_value)/4,2)
    
    breaks_log2FoldChange<-unique(sort(round(c(0,max_abs_value,seq(-1*max_abs_value,max_abs_value, by=step)),1)))
    labels_log2FoldChange<-as.character(breaks_log2FoldChange)
    
  }else{
    
    step<-round(abs(min_abs_value--1*min_abs_value)/4,2)
    
    breaks_log2FoldChange<-unique(sort(round(c(0,min_abs_value,seq(-1*min_abs_value,min_abs_value, by=step)),1)))
    labels_log2FoldChange<-as.character(breaks_log2FoldChange)
    
  }# max_abs_value > min_abs_value
  
  
  cat("step_log2FoldChange\n")
  cat(sprintf(as.character(step)))
  cat("\n")
  cat("labels_log2FoldChange\n")
  cat(sprintf(as.character(labels_log2FoldChange)))
  cat("\n")
  
  
  
  indx_Minus_logpval<-which(colnames(DE_results_subset) == 'Minus_logpval')
  
  A_Minus_logpval<-summary(DE_results_subset[,indx_Minus_logpval])
  
  
  
  cat("A_Minus_logpval\n")
  cat(sprintf(as.character(names(A_Minus_logpval))))
  cat("\n")
  cat(sprintf(as.character(A_Minus_logpval)))
  cat("\n")
  
  
  max_value<-A_Minus_logpval[6]
  min_value<-A_Minus_logpval[1]
  
  
  step<-round(abs(max_value-min_value)/3,0)
  
  if(step == 0)
  {
    
    step<-1
  }
  breaks_Minus_logpval<-unique(sort(unique(c(0,max_value,seq(min_value,max_value, by=step)))))
  labels_Minus_logpval<-as.character(round(breaks_Minus_logpval,1))
  
  
  cat("step_Minus_logpval\n")
  cat(sprintf(as.character(step)))
  cat("\n")
  cat("labels_Minus_logpval\n")
  cat(sprintf(as.character(labels_Minus_logpval)))
  cat("\n")
  
  if(DEBUG == 1)
  {
    cat("Subsample_start:\n")
    
  }
  
  
  DE_results_subset.dt<-data.table(DE_results_subset, key=c("comparison","seurat_cluster"))
  NO_SIG_DE_results_subset<-as.data.frame(DE_results_subset.dt[,.SD[Minus_logpval < 1.3], by=key(DE_results_subset.dt)])
  
  if(DEBUG == 1)
  {
    cat("NO_SIG_DE_results_subset_0\n")
    cat(str(NO_SIG_DE_results_subset))
    cat("\n")
  }
  
  NO_SIG_DE_results_subset.dt<-data.table(NO_SIG_DE_results_subset, key=c("comparison","seurat_cluster"))
  test<-as.data.frame(NO_SIG_DE_results_subset.dt[,.(Symbol=sample(Symbol, round(0.1*length(Symbol)))), by=key(NO_SIG_DE_results_subset.dt)])
  
  test$Flag_subset<-1
  
  if(DEBUG == 1)
  {
    cat("test_0\n")
    cat(str(test))
    cat("\n")
  }
  
  NO_SIG_DE_results_subset<-merge(NO_SIG_DE_results_subset, test, by=c("comparison","seurat_cluster","Symbol"), all.x=T)
  
  if(DEBUG == 1)
  {
    cat("NO_SIG_DE_results_subset_1\n")
    cat(str(NO_SIG_DE_results_subset))
    cat("\n")
  }
  
  
  if(DEBUG == 1)
  {
    cat("Volcano_START:\n")
    
  }
  
  # geom_point(data=DE_results_subset[which(DE_results_subset$Significance == 'NO'),],
  #            color="black",fill="gray", stroke=0.2, shape=21, size=0.75)+
  
  
  volcano_plot<-ggplot(data=DE_results_subset,
                       aes(x=log2FoldChange,
                           y=Minus_logpval)) +
    geom_vline(xintercept=c(-0.25,0.25), color="gray", linetype='dashed',size=0.2)+
    geom_point(data=NO_SIG_DE_results_subset[which(NO_SIG_DE_results_subset$Flag_subset == 1),],
               color="black",fill="gray", stroke=0.2, shape=21, size=0.75)+
    geom_point(data=DE_results_subset[which(DE_results_subset$Significance == 'YES' &
                                       DE_results_subset$log2FoldChange >=0.25),],
               color="black",fill="red", stroke=0.2, shape=21, size=0.75)+
    geom_point(data=DE_results_subset[which(DE_results_subset$Significance == 'YES' &
                                       DE_results_subset$log2FoldChange <= -0.25),],
               color="black",fill="blue", stroke=0.2, shape=21, size=0.75)+
    scale_y_continuous(name='-log10pval',
                       breaks=breaks_Minus_logpval,
                       labels=labels_Minus_logpval,
                       limits=c(breaks_Minus_logpval[1],breaks_Minus_logpval[length(breaks_Minus_logpval)]))+
    scale_x_continuous(name='log2FoldChange',
                       breaks=breaks_log2FoldChange,
                       labels=labels_log2FoldChange,
                       limits=c(breaks_log2FoldChange[1],breaks_log2FoldChange[length(breaks_log2FoldChange)]))
  
  volcano_plot <-volcano_plot+
    theme_cowplot(font_size = 2,
                  font_family = "sans")+
    facet_grid(comparison ~ seurat_cluster, scales='free_x', space='free_x', switch="y", drop=F)+
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
          legend.position="right")+
    ggeasy::easy_center_title()
  
  
  
  volcano_plot <-volcano_plot+
    geom_text_repel(data=DE_results_subset[which(DE_results_subset$Symbol%in%Gather$Symbol &
                                             DE_results_subset$Significance == 'YES' &
                                             DE_results_subset$log2FoldChange >=0.25),],
                     aes(x=log2FoldChange,
                         y=Minus_logpval,
                         label=Symbol,
                         color=Path),
                     family="sans",
                     fontface='italic',
                     segment.size  = 0.2,
                     segment.color = "black",
                     force=25,
                     size=2,
                     xlim=c(0.25,breaks_log2FoldChange[length(breaks_log2FoldChange)]),
                     box.padding = 1,
                     max.overlaps = Inf)+
    geom_text_repel(data=DE_results_subset[which(DE_results_subset$Symbol%in%Gather$Symbol &
                                             DE_results_subset$Significance == 'YES' &
                                             DE_results_subset$log2FoldChange <= -0.25),],
                     aes(x=log2FoldChange,
                         y=Minus_logpval,
                         label=Symbol,
                         color=Path),
                     family="sans",
                     fontface='italic',
                     force=25,
                     segment.size  = 0.2,
                     segment.color = "black",
                     xlim=c(breaks_log2FoldChange[1],-0.25),
                     size=2,
                     box.padding = 1,
                     max.overlaps = Inf)
  
  volcano_plot <-volcano_plot+
    geom_text_repel(data=DE_results_subset[which(DE_results_subset$Symbol== 'CUX1' &
                                                   DE_results_subset$Significance == 'YES' &
                                                   DE_results_subset$log2FoldChange <= -0.25),],
                    aes(x=log2FoldChange,
                        y=Minus_logpval,
                        label=Symbol),
                    color='black',
                    family="sans",
                    fontface='italic',
                    segment.size  = 0.2,
                    segment.color = "black",
                    force=25,
                    size=2,
                    xlim=c(breaks_log2FoldChange[1],-0.25),
                    box.padding = 1,
                    max.overlaps = Inf)
  
  volcano_plot <-volcano_plot+
    geom_text_repel(data=DE_results_subset[which(DE_results_subset$Symbol== 'CUX1' &
                                                   DE_results_subset$Significance == 'YES'  &
                                                   DE_results_subset$log2FoldChange >=0.25),],
                    aes(x=log2FoldChange,
                        y=Minus_logpval,
                        label=Symbol),
                    color='black',
                    family="sans",
                    fontface='italic',
                    segment.size  = 0.2,
                    segment.color = "black",
                    force=25,
                    size=2,
                    xlim=c(0.25,breaks_log2FoldChange[length(breaks_log2FoldChange)]),
                    box.padding = 1,
                    max.overlaps = Inf)
  
  cat("END part volcano_plot\n")
  
  
  setwd(path_graphs)
  
  svgname<-paste(paste("Volcano","selected_comparisons","subset", sep='_'),".svg",sep='')
  
  setwd(path_graphs)
  svglite(svgname, width = 4, height = 4)
  print(volcano_plot)
  dev.off()
  
  
}

printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}



#### main script ----

# make_option(c("--comparison_array"), type="character", default=NULL, 
#             metavar="type", 
#             help="Path to tab-separated input file listing regions to analyze. Required."),

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--platelet_volume_genes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--platelet_genes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--EZH2_signature"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ZEB1_signature"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--AKT_signature"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selected_clusters"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selected_comparisons"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--DE_results"), type="character", default=NULL, 
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
  
  volcano_function(opt)
  volcano_subset_function(opt)

  
}


###########################################################################

system.time( main() )