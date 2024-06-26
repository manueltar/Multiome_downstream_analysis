
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggupset", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
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
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggrepel", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)


data_wrangling_upsetr = function(option_list)
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
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ ORA_ActivePathways_results ----
  
  ORA_ActivePathways_results<-readRDS(file=opt$ORA_ActivePathways_results)
  

  cat("ORA_ActivePathways_results_0\n")
  cat(str(ORA_ActivePathways_results))
  cat("\n")
  
  
  #### READ and transform selected_pathways ----
  
  selected_pathways = unlist(strsplit(opt$selected_pathways, split=","))
  
  cat("selected_pathways_\n")
  cat(sprintf(as.character(selected_pathways)))
  cat("\n")
  
  
  ##### Subset ORA_ActivePathways_results ----
  
  ORA_ActivePathways_results_sel<-ORA_ActivePathways_results[which(ORA_ActivePathways_results$id%in%selected_pathways),]
  
  cat("ORA_ActivePathways_results_sel_0\n")
  cat(str(ORA_ActivePathways_results_sel))
  cat("\n")
  
  
 
  
  #### classification of pathways -----
  
  ORA_ActivePathways_results_sel$CLASS<-NA
  
  MEP_progenitor<-c('GSE15330_LYMPHOID_MULTIPOTENT_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_UP',
                    'GSE15330_LYMPHOID_MULTIPOTENT_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_IKAROS_KO_UP',
                    'GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_UP',
                    'GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_DN',
                    'ZHENG_CORD_BLOOD_C3_MEGAKARYOCYTE_ERYTHROID_PROGENITOR',
                    'ZHENG_CORD_BLOOD_C1_PUTATIVE_MEGAKARYOCYTE_PROGENITOR',
                    'GSE15330_MEGAKARYOCYTE_ERYTHROID_VS_GRANULOCYTE_MONOCYTE_PROGENITOR_IKAROS_KO_DN',
                    'GSE15330_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_VS_PRO_BCELL_UP',
                    'GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_IKAROS_KO_DN')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%MEP_progenitor)]<-"MEP_progenitor"
  
  Megakaryocyte<-c('DESCARTES_FETAL_SPLEEN_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_LUNG_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_KIDNEY_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_HEART_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_ADRENAL_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_MUSCLE_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_LIVER_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_CEREBRUM_MEGAKARYOCYTES',
                   'DESCARTES_MAIN_FETAL_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_PLACENTA_MEGAKARYOCYTES',
                   'REACTOME_FACTORS_INVOLVED_IN_MEGAKARYOCYTE_DEVELOPMENT_AND_PLATELET_PRODUCTION')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Megakaryocyte)]<-"Megakaryocyte"
  
  Platelet<-c('GNATENKO_PLATELET_SIGNATURE',
              'HP_IMPAIRED_COLLAGEN_INDUCED_PLATELET_AGGREGATION',
              'HP_IMPAIRED_PLATELET_AGGREGATION',
              'HP_ABNORMAL_PLATELET_FUNCTION',
              'GAVISH_3CA_MALIGNANT_METAPROGRAM_34_PLATELET_ACTIVATION',
              'RAGHAVACHARI_PLATELET_SPECIFIC_GENES',
              'GOBP_PLATELET_ACTIVATION',
              'HP_ABNORMAL_PLATELET_VOLUME')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Platelet)]<-"Platelet"
  
  
  Erythrocyte<-c("HP_ABNORMAL_ERYTHROCYTE_PHYSIOLOGY")
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Erythrocyte)]<-"Erythrocyte"
  
  Platelet_volume<-c('HP_ABNORMAL_PLATELET_VOLUME')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Platelet_volume)]<-"Platelet_volume"
  
  Mitosis<-c('GOBP_MICROTUBULE_CYTOSKELETON_ORGANIZATION_INVOLVED_IN_MITOSIS',
             'REICHERT_MITOSIS_LIN9_TARGETS')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Mitosis)]<-"Mitosis"
  
  CUX1_target_genes<-"CUX1_TARGET_GENES"
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%CUX1_target_genes)]<-"CUX1_target_genes"
  
  GATA_target_genes<-c('GSE40273_GATA1_KO_VS_WT_TREG_DN',
                       'GATA2_01',
                       'GATA1_02',
                       'HUANG_GATA2_TARGETS_DN',
                       'GSE40277_EOS_AND_LEF1_TRANSDUCED_VS_GATA1_AND_SATB1_TRANSDUCED_CD4_TCELL_DN',
                       'GSE40274_GATA1_VS_FOXP3_AND_GATA1_TRANSDUCED_ACTIVATED_CD4_TCELL_UP',
                       'GATA1_05',
                       'GATA1_03',
                       'GATA1_01')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%GATA_target_genes)]<-"GATA_target_genes"
  
  RUNX1_target_genes<-c('RACCACAR_AML_Q6','REACTOME_RUNX1_REGULATES_GENES_INVOLVED_IN_MEGAKARYOCYTE_DIFFERENTIATION_AND_PLATELET_FUNCTION',
                        'MORF_RUNX1','AML1_Q6','AML1_01','AML_Q6')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%RUNX1_target_genes)]<-"RUNX1_target_genes"
  
  PU1_target_genes<-'PU1_Q6'
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%PU1_target_genes)]<-"PU1_target_genes"
  
  
  AML<-c('HP_ABNORMAL_MYELOID_LEUKOCYTE_MORPHOLOGY',
         'GSE23502_WT_VS_HDC_KO_MYELOID_DERIVED_SUPPRESSOR_CELL_COLON_TUMOR_UP',
         'GSE10325_CD4_TCELL_VS_MYELOID_UP',
         'RUBENSTEIN_SKELETAL_MUSCLE_MYELOID_CELLS',
         'NAKAYA_MYELOID_DENDRITIC_CELL_FLUMIST_AGE_18_50YO_7DY_UP',
         'GSE10325_LUPUS_CD4_TCELL_VS_LUPUS_MYELOID_UP',
         'GSE10325_LUPUS_CD4_TCELL_VS_LUPUS_MYELOID_DN',
         'VALK_AML_CLUSTER_15',
         'VALK_AML_CLUSTER_3',
         'ROSS_ACUTE_MYELOID_LEUKEMIA_CBF',
         'ALCALAY_AML_BY_NPM1_LOCALIZATION_DN')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%AML)]<-"AML_and_myeloid"
  
  
  ORA_ActivePathways_results_sel$CLASS<-droplevels(factor(ORA_ActivePathways_results_sel$CLASS,
                                                           levels=c("MEP_progenitor","Megakaryocyte","Platelet","Platelet_volume","Erythrocyte","Mitosis",
                                                                    "CUX1_target_genes","PU1_target_genes","RUNX1_target_genes","GATA_target_genes","AML_and_myeloid"),
                                                           ordered=T))
  
  cat("ORA_ActivePathways_results_sel_0\n")
  cat(str(ORA_ActivePathways_results_sel))
  cat("\n")
  cat(sprintf(as.character(names(summary(ORA_ActivePathways_results_sel$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(ORA_ActivePathways_results_sel$CLASS))))
  cat("\n")
  
  check.NA<-droplevels(ORA_ActivePathways_results_sel[is.na(ORA_ActivePathways_results_sel$CLASS),])
  
  cat("check.NA_0\n")
  cat(str(check.NA))
  cat("\n")
  cat(sprintf(as.character(names(summary(check.NA$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(check.NA$CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check.NA$id))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(check.NA$id)))))
  cat("\n")
  
  
  
  
  
  
  
  
  
  ####Split LONG ----
  
  ORA_ActivePathways_results_sel_long<-unique(as.data.frame(cSplit(ORA_ActivePathways_results_sel,sep = '|', direction = "long",
                                                                       splitCols = "overlap_symbol"),stringsAsFactors=F))
  
  colnames(ORA_ActivePathways_results_sel_long)[which(colnames(ORA_ActivePathways_results_sel_long) == 'overlap_symbol')]<-'Symbol'
  
  cat("ORA_ActivePathways_results_sel_long_0\n")
  cat(str(ORA_ActivePathways_results_sel_long))
  cat("\n")
  
  
  Subset_symbol<-unique(ORA_ActivePathways_results_sel_long[,c(which(colnames(ORA_ActivePathways_results_sel_long) == 'Symbol'),
                                                          which(colnames(ORA_ActivePathways_results_sel_long) == 'CLASS'))])
  
  Subset_symbol$CLASS<-as.character(Subset_symbol$CLASS)
  
  Subset_symbol$CLASS<-droplevels(factor(Subset_symbol$CLASS,
                              levels=c("PU1_target_genes","Platelet_volume","Erythrocyte","Mitosis","CUX1_target_genes","RUNX1_target_genes","Megakaryocyte","MEP_progenitor","GATA_target_genes","Platelet","AML_and_myeloid"),
                              ordered=T))
  
  Subset_symbol$Presence<-1
  
  cat("Subset_symbol_0\n")
  cat(str(Subset_symbol))
  cat("\n")
  cat(sprintf(as.character(names(summary(Subset_symbol$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Subset_symbol$CLASS))))
  cat("\n")
  
  ######path graphs --------
  
  path_graphs<-paste(out,'upsetR_ORA','/',sep='')
  
  if (file.exists(path_graphs)){
    
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  ######################### subgraph--------------
  
 
  subgraph.dt<-data.table(Subset_symbol,key=c('CLASS'))
  
  cat("subgraph.dt\n")
  cat(str(subgraph.dt))
  cat("\n")
  
  
  subgraph_CLASS_summarised<-as.data.frame(subgraph.dt[,.(instances=.N), by=key(subgraph.dt)],stringsAsFactors=F)
  
  
  
  cat("subgraph_CLASS_summarised\n")
  cat(str(subgraph_CLASS_summarised))
  cat("\n")
  cat(sprintf(as.character(names(summary(subgraph_CLASS_summarised$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(subgraph_CLASS_summarised$CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(subgraph_CLASS_summarised$instances)))))
  cat("\n")
  cat(sprintf(as.character(summary(subgraph_CLASS_summarised$instances))))
  cat("\n")
  
  
  step<-round(max(subgraph_CLASS_summarised$instances)/4,0)
  
  step<-25
  
  cat("--------->\t")
  cat(sprintf(as.character(step)))
  cat("\n")
  
  if(step == 0){
    
    step=1
  }
  breaks.x<-rev(sort(unique(c(max(subgraph_CLASS_summarised$instances),seq(0,max(subgraph_CLASS_summarised$instances), by=step)))))
  labels.x<-as.character(round(breaks.x,0))
  
  
  cat(sprintf(as.character(breaks.x)))
  cat("\n")
  
  vector.fill<-c(brewer.pal(length(levels(subgraph_CLASS_summarised$CLASS)), "Set3"),brewer.pal(length(levels(subgraph_CLASS_summarised$CLASS)), "Set2"))
  
  cat("vector.fill\n")
  cat(str(vector.fill))
  cat("\n")
  
  
  
  subgraph<-ggplot(data=subgraph_CLASS_summarised,
                   aes(y=CLASS,
                       x=instances,
                       fill=CLASS)) +
    geom_bar(stat="identity",colour='black')+
    theme_bw()+
    theme_classic()+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
          axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
          axis.line.x = element_line(size = 0.4),
          axis.line.y = element_line(size = 0.4),
          axis.ticks.x = element_line(size = 0.4),
          axis.ticks.y = element_line(size = 0.4))+
    scale_x_reverse(name="Number of genes",breaks=breaks.x,labels=labels.x,
                    limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
    scale_y_discrete(position="right",name=NULL, drop=T)+
    scale_fill_manual(values=vector.fill, drop=F)+
    theme(legend.position="hidden")+
    ggeasy::easy_center_title()
  
  cat("subgraph genes DONE\n")
  
 
  
  setwd(path_graphs)
  
  svgname<-'test.svg'
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= subgraph,
           device="svg",
           height=3, width=3)
  }
  
  #### UpsetR ---------------
  
  DEBUG<-1
  
  Subset_symbol_wide<-unique(as.data.frame(pivot_wider(Subset_symbol,
                                                         id_cols=Symbol,
                                                         names_from=CLASS,
                                                         values_from=Presence), stringsAsFactors=F))
  
  
  if(DEBUG == 1)
  {
    cat("Subset_symbol_wide_0\n")
    cat(str(Subset_symbol_wide))
    cat("\n")
    cat(str(unique(Subset_symbol_wide$Symbol)))
    cat("\n")
  }
  
  Subset_symbol_wide[is.na(Subset_symbol_wide)]<-0
  
  if(DEBUG == 1)
  {
    cat("Subset_symbol_wide_1\n")
    cat(str(Subset_symbol_wide))
    cat("\n")
    cat(str(unique(Subset_symbol_wide$Symbol)))
    cat("\n")
  }
  
  indx.dep<-c(which(colnames(Subset_symbol_wide) == 'Symbol'))
  Subset_symbol_wide.t<-t(Subset_symbol_wide[,-indx.dep])
  colnames(Subset_symbol_wide.t)<-Subset_symbol_wide$Symbol
  
  
  if(DEBUG == 1)
  {
    cat("Subset_symbol_wide.t_0\n")
    cat(str(Subset_symbol_wide.t))
    cat("\n")
  }
  
  
  
  Subset_symbol_wide.t.logit<-Subset_symbol_wide.t == 1
  
  
  if(DEBUG == 1)
  {
    cat("Subset_symbol_wide.t.logit_0\n")
    cat(str(Subset_symbol_wide.t.logit))
    cat("\n")
  }
  
 
  
  tidy_tag<- Subset_symbol_wide.t.logit %>%
    as_tibble(rownames = "CLASS") %>%
    gather(Symbol, Member, -CLASS) %>%
    filter(Member) %>%
    select(- Member)
  
  UpsetR_plot<-tidy_tag %>%
    group_by(Symbol) %>%
    summarize(CLASSs = list(CLASS)) %>%
    ggplot(aes(x = CLASSs)) +
    geom_bar() +
    scale_x_upset()+
    geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
    theme_combmatrix(
      combmatrix.label.make_space = TRUE,
      combmatrix.label.width = NULL,
      combmatrix.label.height = NULL,
      combmatrix.label.extra_spacing = 3,
      combmatrix.label.total_extra_spacing = unit(10, "pt"),
      combmatrix.label.text = NULL,
      combmatrix.panel.margin = unit(c(1.5, 1.5), "pt"),
      combmatrix.panel.striped_background = TRUE,
      combmatrix.panel.striped_background.color.one = "white",
      combmatrix.panel.striped_background.color.two = "#F7F7F7",
      combmatrix.panel.point.size = 1,
      combmatrix.panel.line.size = 1,
      combmatrix.panel.point.color.fill = "black",
      combmatrix.panel.point.color.empty = "#E0E0E0")
  
  UpsetR_plot<-UpsetR_plot+
    ggtitle(paste("n unique genes =",dim(Subset_symbol_wide)[1],sep=' '))+
    theme_classic()+
    theme(plot.title=element_text(size=8, color="black", family="sans"),
          axis.title.y=element_text(size=8, color="black", family="sans"),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=8, color="black", family="sans"),
          axis.text.x=element_blank(),
          axis.line.x = element_line(size = 0.4),
          axis.line.y = element_line(size = 0.4),
          axis.ticks.x = element_line(size = 0.4),
          axis.ticks.y = element_line(size = 0.4))+
    theme(legend.key.size = unit(0.25, 'cm'), #change legend key size
          legend.key.height = unit(0.25, 'cm'), #change legend key height
          legend.key.width = unit(0.25, 'cm'), #change legend key width
          legend.title = element_text(size=8, family="sans"), #change legend title font size
          legend.text = element_text(size=6, family="sans"),
          legend.position="hidden")+ #change legend text font size
    ggeasy::easy_center_title()
  
  subgraph<-subgraph+
    theme(axis.text.y=element_blank())
  
  subgraph_FINAL<-plot_grid(NULL,subgraph,
                            nrow = 2,
                            ncol=1,
                            rel_heights = c(0.7, 0.325))
  
  
  graph_FINAL<-plot_grid(subgraph_FINAL,UpsetR_plot,
                         nrow = 1,
                         ncol=2,
                         rel_widths=c(0.15,1.05))
  
  cat("UpsetR genes DONE\n")
  
 
  
  setwd(path_graphs)
  
  
  svgname<-paste("UpsetR_","genes",".svg",sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= graph_FINAL,
           device="svg",
           height=6, width=13)
  }
  
}

data_wrangling_lolliplot = function(option_list)
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
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ ORA_ActivePathways_results ----
  
  ORA_ActivePathways_results<-readRDS(file=opt$ORA_ActivePathways_results)
  
  
  cat("ORA_ActivePathways_results_0\n")
  cat(str(ORA_ActivePathways_results))
  cat("\n")
  
  
  #### READ and transform selected_pathways ----
  
  selected_pathways = unlist(strsplit(opt$selected_pathways, split=","))
  
  cat("selected_pathways_\n")
  cat(sprintf(as.character(selected_pathways)))
  cat("\n")
  
  
  ##### Subset ORA_ActivePathways_results ----
  
  ORA_ActivePathways_results_sel<-ORA_ActivePathways_results[which(ORA_ActivePathways_results$id%in%selected_pathways),]
  
  cat("ORA_ActivePathways_results_sel_0\n")
  cat(str(ORA_ActivePathways_results_sel))
  cat("\n")
  
  #### classification of pathways -----
  
  ORA_ActivePathways_results_sel$CLASS<-NA
  
  MEP_progenitor<-c('GSE15330_LYMPHOID_MULTIPOTENT_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_UP',
                    'GSE15330_LYMPHOID_MULTIPOTENT_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_IKAROS_KO_UP',
                    'GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_UP',
                    'GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_DN',
                    'ZHENG_CORD_BLOOD_C3_MEGAKARYOCYTE_ERYTHROID_PROGENITOR',
                    'ZHENG_CORD_BLOOD_C1_PUTATIVE_MEGAKARYOCYTE_PROGENITOR',
                    'GSE15330_MEGAKARYOCYTE_ERYTHROID_VS_GRANULOCYTE_MONOCYTE_PROGENITOR_IKAROS_KO_DN',
                    'GSE15330_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_VS_PRO_BCELL_UP',
                    'GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_IKAROS_KO_DN')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%MEP_progenitor)]<-"MEP_progenitor"
  
  Megakaryocyte<-c('DESCARTES_FETAL_SPLEEN_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_LUNG_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_KIDNEY_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_HEART_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_ADRENAL_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_MUSCLE_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_LIVER_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_CEREBRUM_MEGAKARYOCYTES',
                   'DESCARTES_MAIN_FETAL_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_PLACENTA_MEGAKARYOCYTES',
                   'REACTOME_FACTORS_INVOLVED_IN_MEGAKARYOCYTE_DEVELOPMENT_AND_PLATELET_PRODUCTION')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Megakaryocyte)]<-"Megakaryocyte"
  
  Platelet<-c('GNATENKO_PLATELET_SIGNATURE',
              'HP_IMPAIRED_COLLAGEN_INDUCED_PLATELET_AGGREGATION',
              'HP_IMPAIRED_PLATELET_AGGREGATION',
              'HP_ABNORMAL_PLATELET_FUNCTION',
              'GAVISH_3CA_MALIGNANT_METAPROGRAM_34_PLATELET_ACTIVATION',
              'RAGHAVACHARI_PLATELET_SPECIFIC_GENES',
              'GOBP_PLATELET_ACTIVATION',
              'HP_ABNORMAL_PLATELET_VOLUME')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Platelet)]<-"Platelet"
  
  
  Erythrocyte<-c("HP_ABNORMAL_ERYTHROCYTE_PHYSIOLOGY")
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Erythrocyte)]<-"Erythrocyte"
  
  Platelet_volume<-c('HP_ABNORMAL_PLATELET_VOLUME')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Platelet_volume)]<-"Platelet_volume"
  
  Mitosis<-c('GOBP_MICROTUBULE_CYTOSKELETON_ORGANIZATION_INVOLVED_IN_MITOSIS',
             'REICHERT_MITOSIS_LIN9_TARGETS')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Mitosis)]<-"Mitosis"
  
  CUX1_target_genes<-"CUX1_TARGET_GENES"
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%CUX1_target_genes)]<-"CUX1_target_genes"
  
  GATA_target_genes<-c('GSE40273_GATA1_KO_VS_WT_TREG_DN',
                       'GATA2_01',
                       'GATA1_02',
                       'HUANG_GATA2_TARGETS_DN',
                       'GSE40277_EOS_AND_LEF1_TRANSDUCED_VS_GATA1_AND_SATB1_TRANSDUCED_CD4_TCELL_DN',
                       'GSE40274_GATA1_VS_FOXP3_AND_GATA1_TRANSDUCED_ACTIVATED_CD4_TCELL_UP',
                       'GATA1_05',
                       'GATA1_03',
                       'GATA1_01')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%GATA_target_genes)]<-"GATA_target_genes"
  
  RUNX1_target_genes<-c('RACCACAR_AML_Q6','REACTOME_RUNX1_REGULATES_GENES_INVOLVED_IN_MEGAKARYOCYTE_DIFFERENTIATION_AND_PLATELET_FUNCTION',
                        'MORF_RUNX1','AML1_Q6','AML1_01','AML_Q6')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%RUNX1_target_genes)]<-"RUNX1_target_genes"
  
  PU1_target_genes<-'PU1_Q6'
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%PU1_target_genes)]<-"PU1_target_genes"
  
  
  AML<-c('HP_ABNORMAL_MYELOID_LEUKOCYTE_MORPHOLOGY',
         'GSE23502_WT_VS_HDC_KO_MYELOID_DERIVED_SUPPRESSOR_CELL_COLON_TUMOR_UP',
         'GSE10325_CD4_TCELL_VS_MYELOID_UP',
         'RUBENSTEIN_SKELETAL_MUSCLE_MYELOID_CELLS',
         'NAKAYA_MYELOID_DENDRITIC_CELL_FLUMIST_AGE_18_50YO_7DY_UP',
         'GSE10325_LUPUS_CD4_TCELL_VS_LUPUS_MYELOID_UP',
         'GSE10325_LUPUS_CD4_TCELL_VS_LUPUS_MYELOID_DN',
         'VALK_AML_CLUSTER_15',
         'VALK_AML_CLUSTER_3',
         'ROSS_ACUTE_MYELOID_LEUKEMIA_CBF',
         'ALCALAY_AML_BY_NPM1_LOCALIZATION_DN')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%AML)]<-"AML_and_myeloid"
  
  
  ORA_ActivePathways_results_sel$CLASS<-droplevels(factor(ORA_ActivePathways_results_sel$CLASS,
                                                          levels=c("MEP_progenitor","Megakaryocyte","Platelet","Platelet_volume","Erythrocyte","Mitosis",
                                                                   "CUX1_target_genes","PU1_target_genes","RUNX1_target_genes","GATA_target_genes","AML_and_myeloid"),
                                                          ordered=T))
  
  cat("ORA_ActivePathways_results_sel_0\n")
  cat(str(ORA_ActivePathways_results_sel))
  cat("\n")
  cat(sprintf(as.character(names(summary(ORA_ActivePathways_results_sel$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(ORA_ActivePathways_results_sel$CLASS))))
  cat("\n")
  
  check.NA<-droplevels(ORA_ActivePathways_results_sel[is.na(ORA_ActivePathways_results_sel$CLASS),])
  
  cat("check.NA_0\n")
  cat(str(check.NA))
  cat("\n")
  cat(sprintf(as.character(names(summary(check.NA$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(check.NA$CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check.NA$id))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(check.NA$id)))))
  cat("\n")
  
  
  
  
  
  
  
  
  
  ##### Select only SIG results -----
  
  DEBUG <- 1
  ORA_ActivePathways_results_sel_SIG<-ORA_ActivePathways_results_sel[which(ORA_ActivePathways_results_sel$Significance == 'YES'),]
  
  
  
  ORA_ActivePathways_results_sel_SIG<-ORA_ActivePathways_results_sel_SIG[order(ORA_ActivePathways_results_sel_SIG$CLASS,
                                                                                 ORA_ActivePathways_results_sel_SIG$id,
                                                                                 ORA_ActivePathways_results_sel_SIG$comparison),]
  
  if(DEBUG == 1)
  {
    cat("ORA_ActivePathways_results_sel_SIG_0\n")
    cat(str(ORA_ActivePathways_results_sel_SIG))
    cat("\n")
  }
  
  ORA_ActivePathways_results_sel_SIG$DUMMY<-interaction(ORA_ActivePathways_results_sel_SIG$id,
                                                         ORA_ActivePathways_results_sel_SIG$seurat_cluster,
                                                         ORA_ActivePathways_results_sel_SIG$comparison)
  
  if(DEBUG == 1)
  {
    cat("ORA_ActivePathways_results_sel_SIG$DUMMY\n")
    cat(str(unique(ORA_ActivePathways_results_sel_SIG$DUMMY)))
    cat("\n")
  }
  
  ordered_DUMMY<-as.character(ORA_ActivePathways_results_sel_SIG$DUMMY)
  
  ORA_ActivePathways_results_sel_SIG$DUMMY<-factor(ORA_ActivePathways_results_sel_SIG$DUMMY,
                                          levels=rev(ordered_DUMMY),
                                          ordered=T)
  
  if(DEBUG == 1)
  {
    cat("ORA_ActivePathways_results_sel_SIG_1\n")
    cat(str(ORA_ActivePathways_results_sel_SIG))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(ORA_ActivePathways_results_sel_SIG$id))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(ORA_ActivePathways_results_sel_SIG$id)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(ORA_ActivePathways_results_sel_SIG$DUMMY))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(ORA_ActivePathways_results_sel_SIG$DUMMY)))))
    cat("\n")
  }
  
  
  
  
  ### graph parameters_Minus_logpval
  
  indx_Minus_logpval<-which(colnames(ORA_ActivePathways_results_sel_SIG) == 'Minus_logpval')
  
  A_Minus_logpval<-summary(ORA_ActivePathways_results_sel_SIG[,indx_Minus_logpval])
  
  
  if(DEBUG == 1)
  {
    cat("A_Minus_logpval\n")
    cat(sprintf(as.character(names(A_Minus_logpval))))
    cat("\n")
    cat(sprintf(as.character(A_Minus_logpval)))
    cat("\n")
  }
  
  max_value<-A_Minus_logpval[6]
  min_value<-0
  
  
  step<-round(abs(max_value-min_value)/3,1)
  
  if(step == 0)
  {
    
    step<-1
  }
  breaks_Minus_logpval<-unique(sort(unique(c(0,max_value,seq(min_value,max_value, by=step)))))
  labels_Minus_logpval<-as.character(round(breaks_Minus_logpval,1))
  
  if(DEBUG == 1)
  {
    cat("step_Minus_logpval\n")
    cat(sprintf(as.character(step)))
    cat("\n")
    cat("labels_Minus_logpval\n")
    cat(sprintf(as.character(labels_Minus_logpval)))
    cat("\n")
  }
  
  
  
  vector.fill<-c(brewer.pal(8, "Dark2"),
                   brewer.pal(12, "Set3"),
                   brewer.pal(8, "Set1"))
  
  if(DEBUG == 1)
  {
    cat("vector.fill\n")
    cat(str(vector.fill))
    cat("\n")
  }
  
  breaks_gene_sets<-as.numeric(ORA_ActivePathways_results_sel_SIG$DUMMY)
  labels_gene_sets<-as.character(gsub("\\..+$","",ORA_ActivePathways_results_sel_SIG$DUMMY))
  
  if(DEBUG == 1)
  {
    cat("breaks_gene_sets\n")
    cat(sprintf(as.character(breaks_gene_sets)))
    cat("\n")
    cat("labels_gene_sets\n")
    cat(sprintf(as.character(labels_gene_sets)))
    cat("\n")
  }
  
  #### barplot
  
  if(DEBUG == 1)
  {
    cat("Lolliplot_START:\n")
    
  }
  
  
  Gene_set_lolliplot<-ggplot(data=ORA_ActivePathways_results_sel_SIG, 
                             aes(y=as.numeric(DUMMY),
                                 x=Minus_logpval)) +
    geom_hline(yintercept=as.numeric(ORA_ActivePathways_results_sel_SIG$DUMMY), color="gray", linetype='dashed',size=0.1)+
    geom_segment(data=ORA_ActivePathways_results_sel_SIG,
                 aes(y=as.numeric(DUMMY),
                     yend=as.numeric(DUMMY),
                     x=min(breaks_Minus_logpval),
                     xend=Minus_logpval),color="black", size=0.4)+
    geom_point(data=ORA_ActivePathways_results_sel_SIG,
               aes(fill=comparison,
                   color=comparison), size=3, stroke=1, shape=21)+
    geom_text(data=ORA_ActivePathways_results_sel_SIG,
              aes(x=Minus_logpval, y=as.numeric(DUMMY), label=n_genes_in_overlap),color="black",size=2, family="sans",fontface="bold")+
    scale_fill_manual(values=vector.fill, drop=F)+
    scale_color_manual(values=vector.fill, drop=F)
  
  if(DEBUG == 1)
  {
    cat("Lolliplot_facet_grid:\n")
    
  }
  
 
 
  
  Gene_set_lolliplot <-Gene_set_lolliplot+
    theme_cowplot(font_size = 2,
                  font_family = "sans")+
    facet_grid(. ~ seurat_cluster+comparison, scales='free_x', space='free_x', switch="y", drop=F)+
      scale_y_continuous(name=NULL, breaks=breaks_gene_sets,
                         labels=labels_gene_sets)+
    scale_x_continuous(name='-log10pval',
                       breaks=breaks_Minus_logpval,
                       labels=labels_Minus_logpval)+
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
          axis.title.y=element_blank(),
          axis.title.x=element_text(size=8,color="black", family="sans"),
          axis.text.y=element_text(size=6,color="black", family="sans", face='bold'),
          axis.text.x=element_text(size=6,color="black", family="sans"),
          axis.line.x = element_line(size = 0.4),
          axis.ticks.x = element_line(size = 0.4),
          axis.ticks.y = element_line(size = 0.4),
          axis.line.y = element_line(size = 0.4))+
    theme(legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.key.size = unit(0.35, 'cm'), #change legend key size
          legend.key.height = unit(0.35, 'cm'), #change legend key height
          legend.key.width = unit(0.35, 'cm'), #change legend key width
          legend.position="bottom")+
    guides(fill=guide_legend(nrow=1,byrow=TRUE))+
    ggeasy::easy_center_title()
  
  if(DEBUG == 1)
  {
    cat("Lolliplot_END:\n")
    
  }
  
  path_graphs<-paste(out,'Lolliplot','/',sep='')
  
  if (file.exists(path_graphs)){
    
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  setwd(path_graphs)
  
  svgname<-paste(paste("Lolliplot",'Figure',sep='_'),".svg",sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= Gene_set_lolliplot,
           device="svg",
           height=13, width=13)
  }
  
  if(DEBUG == 1)
  {
    cat("THE_END:\n")
    
  }
  
  
}

data_wrangling_classification_of_genes = function(option_list)
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
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ ORA_ActivePathways_results ----
  
  ORA_ActivePathways_results<-readRDS(file=opt$ORA_ActivePathways_results)
  
  
  cat("ORA_ActivePathways_results_0\n")
  cat(str(ORA_ActivePathways_results))
  cat("\n")
  
  
  #### READ and transform selected_pathways ----
  
  selected_pathways = unlist(strsplit(opt$selected_pathways, split=","))
  
  cat("selected_pathways_\n")
  cat(sprintf(as.character(selected_pathways)))
  cat("\n")
  
  
  ##### Subset ORA_ActivePathways_results ----
  
  ORA_ActivePathways_results_sel<-ORA_ActivePathways_results[which(ORA_ActivePathways_results$id%in%selected_pathways),]
  
  cat("ORA_ActivePathways_results_sel_0\n")
  cat(str(ORA_ActivePathways_results_sel))
  cat("\n")
  
  #### classification of pathways -----
  
  ORA_ActivePathways_results_sel$CLASS<-NA
  
  MEP_progenitor<-c('GSE15330_LYMPHOID_MULTIPOTENT_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_UP',
                    'GSE15330_LYMPHOID_MULTIPOTENT_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_IKAROS_KO_UP',
                    'GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_UP',
                    'GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_DN',
                    'ZHENG_CORD_BLOOD_C3_MEGAKARYOCYTE_ERYTHROID_PROGENITOR',
                    'ZHENG_CORD_BLOOD_C1_PUTATIVE_MEGAKARYOCYTE_PROGENITOR',
                    'GSE15330_MEGAKARYOCYTE_ERYTHROID_VS_GRANULOCYTE_MONOCYTE_PROGENITOR_IKAROS_KO_DN',
                    'GSE15330_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_VS_PRO_BCELL_UP',
                    'GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_IKAROS_KO_DN')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%MEP_progenitor)]<-"MEP_progenitor"
  
  Megakaryocyte<-c('DESCARTES_FETAL_SPLEEN_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_LUNG_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_KIDNEY_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_HEART_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_ADRENAL_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_MUSCLE_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_LIVER_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_CEREBRUM_MEGAKARYOCYTES',
                   'DESCARTES_MAIN_FETAL_MEGAKARYOCYTES',
                   'DESCARTES_FETAL_PLACENTA_MEGAKARYOCYTES',
                   'REACTOME_FACTORS_INVOLVED_IN_MEGAKARYOCYTE_DEVELOPMENT_AND_PLATELET_PRODUCTION')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Megakaryocyte)]<-"Megakaryocyte"
  
  Platelet<-c('GNATENKO_PLATELET_SIGNATURE',
              'HP_IMPAIRED_COLLAGEN_INDUCED_PLATELET_AGGREGATION',
              'HP_IMPAIRED_PLATELET_AGGREGATION',
              'HP_ABNORMAL_PLATELET_FUNCTION',
              'GAVISH_3CA_MALIGNANT_METAPROGRAM_34_PLATELET_ACTIVATION',
              'RAGHAVACHARI_PLATELET_SPECIFIC_GENES',
              'GOBP_PLATELET_ACTIVATION',
              'HP_ABNORMAL_PLATELET_VOLUME')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Platelet)]<-"Platelet"
  
  
  Erythrocyte<-c("HP_ABNORMAL_ERYTHROCYTE_PHYSIOLOGY")
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Erythrocyte)]<-"Erythrocyte"
  
  Platelet_volume<-c('HP_ABNORMAL_PLATELET_VOLUME')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Platelet_volume)]<-"Platelet_volume"
  
  Mitosis<-c('GOBP_MICROTUBULE_CYTOSKELETON_ORGANIZATION_INVOLVED_IN_MITOSIS',
             'REICHERT_MITOSIS_LIN9_TARGETS')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%Mitosis)]<-"Mitosis"
  
  CUX1_target_genes<-"CUX1_TARGET_GENES"
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%CUX1_target_genes)]<-"CUX1_target_genes"
  
  GATA_target_genes<-c('GSE40273_GATA1_KO_VS_WT_TREG_DN',
                       'GATA2_01',
                       'GATA1_02',
                       'HUANG_GATA2_TARGETS_DN',
                       'GSE40277_EOS_AND_LEF1_TRANSDUCED_VS_GATA1_AND_SATB1_TRANSDUCED_CD4_TCELL_DN',
                       'GSE40274_GATA1_VS_FOXP3_AND_GATA1_TRANSDUCED_ACTIVATED_CD4_TCELL_UP',
                       'GATA1_05',
                       'GATA1_03',
                       'GATA1_01')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%GATA_target_genes)]<-"GATA_target_genes"
  
  RUNX1_target_genes<-c('RACCACAR_AML_Q6','REACTOME_RUNX1_REGULATES_GENES_INVOLVED_IN_MEGAKARYOCYTE_DIFFERENTIATION_AND_PLATELET_FUNCTION',
                        'MORF_RUNX1','AML1_Q6','AML1_01','AML_Q6')
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%RUNX1_target_genes)]<-"RUNX1_target_genes"
  
  PU1_target_genes<-'PU1_Q6'
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%PU1_target_genes)]<-"PU1_target_genes"
  
  
  AML<-c('HP_ABNORMAL_MYELOID_LEUKOCYTE_MORPHOLOGY',
         'GSE23502_WT_VS_HDC_KO_MYELOID_DERIVED_SUPPRESSOR_CELL_COLON_TUMOR_UP',
         'GSE10325_CD4_TCELL_VS_MYELOID_UP',
         'RUBENSTEIN_SKELETAL_MUSCLE_MYELOID_CELLS',
         'NAKAYA_MYELOID_DENDRITIC_CELL_FLUMIST_AGE_18_50YO_7DY_UP',
         'GSE10325_LUPUS_CD4_TCELL_VS_LUPUS_MYELOID_UP',
         'GSE10325_LUPUS_CD4_TCELL_VS_LUPUS_MYELOID_DN',
         'VALK_AML_CLUSTER_15',
         'VALK_AML_CLUSTER_3',
         'ROSS_ACUTE_MYELOID_LEUKEMIA_CBF',
         'ALCALAY_AML_BY_NPM1_LOCALIZATION_DN')
 
 
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%AML)]<-"AML_and_myeloid"
  
  AKT_progenitor<-c('ALCALAY_AML_BY_NPM1_LOCALIZATION_DN',
                    'AKT_UP_MTOR_DN.V1_UP',
                    'AKT_UP.V1_DN',
                    'HALLMARK_PI3K_AKT_MTOR_SIGNALING',
                    'XU_HGF_SIGNALING_NOT_VIA_AKT1_6HR')
  
  
  ORA_ActivePathways_results_sel$CLASS[which(ORA_ActivePathways_results_sel$id%in%AKT_progenitor)]<-"AKT_progenitor"
  
  
  ORA_ActivePathways_results_sel$CLASS<-droplevels(factor(ORA_ActivePathways_results_sel$CLASS,
                                                          levels=c("MEP_progenitor","Megakaryocyte","AKT_progenitor","Platelet","Platelet_volume","Erythrocyte","Mitosis",
                                                                   "CUX1_target_genes","PU1_target_genes","RUNX1_target_genes","GATA_target_genes","AML_and_myeloid"),
                                                          ordered=T))
  
  cat("ORA_ActivePathways_results_sel_0\n")
  cat(str(ORA_ActivePathways_results_sel))
  cat("\n")
  cat(sprintf(as.character(names(summary(ORA_ActivePathways_results_sel$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(ORA_ActivePathways_results_sel$CLASS))))
  cat("\n")
  
  check.NA<-droplevels(ORA_ActivePathways_results_sel[is.na(ORA_ActivePathways_results_sel$CLASS),])
  
  cat("check.NA_0\n")
  cat(str(check.NA))
  cat("\n")
  cat(sprintf(as.character(names(summary(check.NA$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(check.NA$CLASS))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check.NA$id))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(check.NA$id)))))
  cat("\n")
  
  
  
  
  
  
  
  
  
  ####Split LONG ----
  
  ORA_ActivePathways_results_sel_long<-unique(as.data.frame(cSplit(ORA_ActivePathways_results_sel,sep = '|', direction = "long",
                                                                    splitCols = "overlap_symbol"),stringsAsFactors=F))
  
  colnames(ORA_ActivePathways_results_sel_long)[which(colnames(ORA_ActivePathways_results_sel_long) == 'overlap_symbol')]<-'Symbol'
  
  cat("ORA_ActivePathways_results_sel_long_0\n")
  cat(str(ORA_ActivePathways_results_sel_long))
  cat("\n")
  
  
  Subset_symbol<-unique(ORA_ActivePathways_results_sel_long[,c(which(colnames(ORA_ActivePathways_results_sel_long) == 'Symbol'),
                                                                which(colnames(ORA_ActivePathways_results_sel_long) == 'CLASS'))])
  
  Subset_symbol$CLASS<-as.character(Subset_symbol$CLASS)
  
  Subset_symbol$CLASS<-factor(Subset_symbol$CLASS,
                              levels=c("PU1_target_genes","Platelet_volume","Erythrocyte","Mitosis","AKT_progenitor","CUX1_target_genes","RUNX1_target_genes","Megakaryocyte","MEP_progenitor","GATA_target_genes","Platelet","AML_and_myeloid"),
                              ordered=T)
  
  cat("Subset_symbol_0\n")
  cat(str(Subset_symbol))
  cat("\n")
  cat(sprintf(as.character(names(summary(Subset_symbol$CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Subset_symbol$CLASS))))
  cat("\n")
  
  #### collapse CLASS ----
  
  Subset_symbol.dt<-data.table(Subset_symbol, key="Symbol")
  
  Subset_symbol_collapsed<-as.data.frame(Subset_symbol.dt[,.(string_CLASS=paste(unique(sort(CLASS)), collapse='|')), by=key(Subset_symbol.dt)] , stringsAsFactors=F)
  
  cat("Subset_symbol_collapsed_0\n")
  cat(str(Subset_symbol_collapsed))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Subset_symbol_collapsed$string_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Subset_symbol_collapsed$string_CLASS)))))
  cat("\n")
  
 
 
 
 
  
  
  CUX1_target_genes<-c('CUX1_target_genes|RUNX1_target_genes|AML_and_myeloid',
                       'CUX1_target_genes|AML_and_myeloid',
                       'CUX1_target_genes|GATA_target_genes|AML_and_myeloid',
                       'CUX1_target_genes|Megakaryocyte|AML_and_myeloid',
                       'CUX1_target_genes|MEP_progenitor|AML_and_myeloid',
                       'CUX1_target_genes|Platelet|AML_and_myeloid',
    'CUX1_target_genes|MEP_progenitor|GATA_target_genes',
                        'CUX1_target_genes|MEP_progenitor|Platelet',
    'CUX1_target_genes|RUNX1_target_genes|MEP_progenitor|AML_and_myeloid','CUX1_target_genes|RUNX1_target_genes|Platelet|AML_and_myeloid',
                        'CUX1_target_genes|Platelet',
                        'CUX1_target_genes|RUNX1_target_genes|MEP_progenitor',
                        'CUX1_target_genes|RUNX1_target_genes|Platelet',
                        'CUX1_target_genes','CUX1_target_genes|GATA_target_genes','CUX1_target_genes|Megakaryocyte','CUX1_target_genes|Megakaryocyte|GATA_target_genes','CUX1_target_genes|Megakaryocyte|Platelet','CUX1_target_genes|MEP_progenitor','CUX1_target_genes|RUNX1_target_genes','CUX1_target_genes|RUNX1_target_genes|Megakaryocyte|Platelet','Erythrocyte|CUX1_target_genes|MEP_progenitor')
  
  MEP_progenitor<-c('MEP_progenitor|AML_and_myeloid',"AKT_progenitor|MEP_progenitor","AKT_progenitor|MEP_progenitor|Platelet",
                    'Megakaryocyte|MEP_progenitor|Platelet|AML_and_myeloid',
                    'Megakaryocyte|MEP_progenitor|AML_and_myeloid',
                    'MEP_progenitor|GATA_target_genes|AML_and_myeloid',
                    'MEP_progenitor|Platelet|AML_and_myeloid',
                    'MEP_progenitor|GATA_target_genes|Platelet',
                    'Erythrocyte|MEP_progenitor','Erythrocyte|MEP_progenitor|GATA_target_genes','Erythrocyte|Megakaryocyte|MEP_progenitor|GATA_target_genes|Platelet','Erythrocyte|RUNX1_target_genes|Megakaryocyte|MEP_progenitor|GATA_target_genes','Erythrocyte|RUNX1_target_genes|MEP_progenitor','Megakaryocyte|MEP_progenitor','Megakaryocyte|MEP_progenitor|Platelet','MEP_progenitor','MEP_progenitor|GATA_target_genes','MEP_progenitor|Platelet','Platelet_volume|Erythrocyte|Megakaryocyte|MEP_progenitor|Platelet','Platelet_volume|Megakaryocyte|MEP_progenitor|GATA_target_genes|Platelet','Platelet_volume|Megakaryocyte|MEP_progenitor|Platelet','Platelet_volume|MEP_progenitor','Platelet_volume|MEP_progenitor|Platelet','RUNX1_target_genes|Megakaryocyte|MEP_progenitor|Platelet','RUNX1_target_genes|MEP_progenitor','RUNX1_target_genes|MEP_progenitor|Platelet')
  
  PU1_target_genes<-c('PU1_target_genes|RUNX1_target_genes','PU1_target_genes','PU1_target_genes|GATA_target_genes','PU1_target_genes|Megakaryocyte','PU1_target_genes|Platelet','PU1_target_genes|RUNX1_target_genes|Platelet',
                      'PU1_target_genes|AML_and_myeloid','PU1_target_genes|RUNX1_target_genes|Platelet|AML_and_myeloid')
  
  Platelet_volume<-c("Platelet_volume|RUNX1_target_genes|Megakaryocyte|MEP_progenitor|Platelet","Platelet_volume|GATA_target_genes|Platelet",'Platelet_volume','Platelet_volume|Megakaryocyte','Platelet_volume|Megakaryocyte|GATA_target_genes|Platelet','Platelet_volume|Megakaryocyte|Platelet','Platelet_volume|Platelet','Platelet_volume|RUNX1_target_genes|GATA_target_genes|Platelet')
  
  Megakaryocyte<-c('Megakaryocyte|AML_and_myeloid',
                   'Megakaryocyte|GATA_target_genes|AML_and_myeloid',
                   'Megakaryocyte|Platelet|AML_and_myeloid','Megakaryocyte','Megakaryocyte|GATA_target_genes','Megakaryocyte|Platelet','Erythrocyte|RUNX1_target_genes|Megakaryocyte|Platelet')
  
  Platelet<-c('GATA_target_genes|Platelet|AML_and_myeloid',
              'Platelet|AML_and_myeloid',
              'GATA_target_genes|Platelet','Platelet','RUNX1_target_genes|Platelet')
  
  Erythrocyte<-c('Erythrocyte|RUNX1_target_genes|Megakaryocyte|AML_and_myeloid','Erythrocyte','Erythrocyte|GATA_target_genes','Erythrocyte|RUNX1_target_genes','Erythrocyte|RUNX1_target_genes|Platelet','Erythrocyte|CUX1_target_genes',
                 'Erythrocyte|AML_and_myeloid','Erythrocyte|CUX1_target_genes|AML_and_myeloid','Erythrocyte|MEP_progenitor|AML_and_myeloid',
                 'Erythrocyte|Megakaryocyte')
  
  Mitosis<-c('Mitosis','Mitosis|AML_and_myeloid','Mitosis|GATA_target_genes|Platelet','Mitosis|MEP_progenitor')
  
  
  RUNX1_target_genes<-c("AKT_progenitor|RUNX1_target_genes",'RUNX1_target_genes','RUNX1_target_genes|GATA_target_genes','RUNX1_target_genes|AML_and_myeloid','RUNX1_target_genes|MEP_progenitor|AML_and_myeloid')
  
  AML_myeloid<-c('AML_and_myeloid',"AKT_progenitor|AML_and_myeloid")
  
  AKT_progenitor<-'AKT_progenitor'
  
  
  GATA_target_genes<-c('GATA_target_genes','GATA_target_genes|AML_and_myeloid')
  
  
  
  Subset_symbol_collapsed$GENE_CLASS<-NA
  
  
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%CUX1_target_genes)]<-'CUX1_target_genes'
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%MEP_progenitor)]<-'MEP_progenitor'
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%PU1_target_genes)]<-'PU1_target_genes'
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%Platelet_volume)]<-'Platelet_volume'
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%Megakaryocyte)]<-'Megakaryocyte'
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%Platelet)]<-'Platelet'
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%Erythrocyte)]<-'Erythrocyte'
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%RUNX1_target_genes)]<-'RUNX1_target_genes'
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%GATA_target_genes)]<-'GATA_target_genes'
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%AML_myeloid)]<-'AML_myeloid'
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%Mitosis)]<-'Mitosis'
  Subset_symbol_collapsed$GENE_CLASS[which(Subset_symbol_collapsed$string_CLASS%in%AKT_progenitor)]<-'AKT_progenitor'
  
  
  
  Subset_symbol_collapsed$GENE_CLASS<-droplevels(factor(Subset_symbol_collapsed$GENE_CLASS,
                                                levels=c("MEP_progenitor","Megakaryocyte","Platelet","Platelet_volume","Erythrocyte","Mitosis",'AKT_progenitor',
                                                         "CUX1_target_genes","PU1_target_genes","RUNX1_target_genes","GATA_target_genes",'AML_myeloid'),
                                                ordered=T))
  
  cat("Subset_symbol_collapsed_0\n")
  cat(str(Subset_symbol_collapsed))
  cat("\n")
  cat(sprintf(as.character(names(summary(Subset_symbol_collapsed$GENE_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Subset_symbol_collapsed$GENE_CLASS))))
  cat("\n")
  
  check.NA<-droplevels(Subset_symbol_collapsed[is.na(Subset_symbol_collapsed$GENE_CLASS),])
  
  cat("check.NA_0\n")
  cat(str(check.NA))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(check.NA$string_CLASS))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(check.NA$string_CLASS)))))
  cat("\n")
  
  
  #### SAVE RDS ----
  
  
  
  setwd(out)
  
  saveRDS(Subset_symbol_collapsed,file='Selected_genes_classified.rds')
  write.table(Subset_symbol_collapsed,file='Selected_genes_classified.tsv',sep="\t",quote=F,row.names = F)
 

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
    make_option(c("--ORA_ActivePathways_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selected_pathways"), type="character", default=NULL, 
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
  
  data_wrangling_upsetr(opt)
  data_wrangling_lolliplot(opt)
  data_wrangling_classification_of_genes(opt)

}


###########################################################################

system.time( main() )
  