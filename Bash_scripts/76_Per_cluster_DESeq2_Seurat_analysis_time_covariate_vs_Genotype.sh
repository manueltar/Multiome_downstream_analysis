#!/bin/bash

eval "$(conda shell.bash hook)"
 

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Master_path_analysis=$(echo "$output_dir")

#rm -rf $Master_path_analysis
#mkdir -p $Master_path_analysis



Log_files=$(echo "$Master_path_analysis""/""Log_files/")

#rm -rf $Log_files
#mkdir -p $Log_files



conda activate multiome_downstream


######################################################## LOOP #####################################################

declare -a arr

cluster_array=$(echo '1,2,3,4,5,6,7,8,9,10,11,12,13')
 
 a=($(echo "$cluster_array" | tr "," '\n'))

 declare -a arr

 for i  in "${a[@]}"
 do

     cluster_array_sel=$i
     echo "$cluster_array_sel"

     ### DE_test


     DE_route=$(echo "$Master_path_analysis""$cluster_array_sel""/")
     echo "$DE_route"

     rm -rf $DE_route
     mkdir -p $DE_route
     
     type=$(echo "$cluster_array_sel""_""DE_test")
     outfile_DE_test=$(echo "$Log_files""outfile_1_""$type"".log")
     touch $outfile_DE_test
     echo -n "" > $outfile_DE_test
     name_DE_test=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")



     Rscript_DE_test=$(echo "$Rscripts_path""241_Per_cluster_Seurat_to_DESeq2_multiome_GeneEXP_time_covariate.R")

     Seurat_object=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/merged_downstream_analysis_FOR_DANIELA.rds")
     
     myjobid_DE_test=$(sbatch --job-name $name_DE_test --output=$outfile_DE_test --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=6 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_DE_test --Seurat_object $Seurat_object --seurat_cluster $cluster_array_sel --type $type --out $DE_route")
     myjobid_seff_DE_test=$(sbatch --dependency=afterany:$myjobid_DE_test --open-mode=append --output=$outfile_DE_test --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_DE_test >> $outfile_DE_test")

     
     ### MySigDB_GSEA

  
     type=$(echo "$cluster_array_sel""_""MySigDB_GSEA")
     outfile_MySigDB_GSEA=$(echo "$Log_files""outfile_2_""$type"".log")
     touch $outfile_MySigDB_GSEA
     echo -n "" > $outfile_MySigDB_GSEA
     name_MySigDB_GSEA=$(echo "$type""_job")

 
     Rscript_MySigDB_GSEA=$(echo "$Rscripts_path""242_Per_cluster__MySigDb_DESeq2_Seurat.R")

     DE_route=$(echo "$Master_path_analysis""$cluster_array_sel""/")
     multiome_edgeR_results=$(echo "$DE_route""DE_genes.tsv")
     path_to_GMT=$(echo "/home/manuel.tardaguila/GMT_files/msigdb_v2023.1.Hs_files_to_download_locally_ENTREZ/")
     search_terms=$(echo "PLATELET,ERYTHROCYTE,CUX1,MEGAKARYOCYTE,GATA1,GATA2,TET2,RUNX1,RUNX2,MITOSIS,ANEUPLOIDY,CYTOKINESIS,MYELOID,AML,HEPATOCYTE,NEURON,LIPID,SPHINGOSINE,FOXM1,SPI1,PU1,WP_PI3K_AKT_SIGNALING_PATHWAY,PI3K,AKT")   # ADD HSC TERMS
     background_genes=$(echo "/home/manuel.tardaguila/GMT_files/msigdb_v2023.1.Hs_files_to_download_locally_ENTREZ/c5.hpo.v2023.2.Hs.entrez.gmt")
     maxGS_size=$(echo "500")
     minGS_size=$(echo "10")
     pval_threshold=$(echo "0.05")
     log2FC_threshold=$(echo "0.25")
     seff_name=$(echo "seff""_""$type")

     # --dependency=afterany:$myjobid_DE_test
     
     myjobid_MySigDB_GSEA=$(sbatch --dependency=afterany:$myjobid_DE_test --job-name=$name_MySigDB_GSEA --output=$outfile_MySigDB_GSEA --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="Rscript $Rscript_MySigDB_GSEA --multiome_edgeR_results $multiome_edgeR_results --path_to_GMT $path_to_GMT --search_terms $search_terms --background_genes $background_genes --maxGS_size $maxGS_size --minGS_size $minGS_size --seurat_cluster $cluster_array_sel --pval_threshold $pval_threshold --log2FC_threshold $log2FC_threshold --type $type --out $DE_route")
     myjobid_seff_MySigDB_GSEA=$(sbatch --dependency=afterany:$myjobid_MySigDB_GSEA --open-mode=append --output=$outfile_MySigDB_GSEA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_MySigDB_GSEA >> $outfile_MySigDB_GSEA")



     echo "->>>$myjobid_MySigDB_GSEA"
     arr[${#arr[@]}]="$myjobid_MySigDB_GSEA"

 done
      


conda deactivate

 done_string=$(echo "--dependency=afterany:"""""${arr[@]}"""")
 echo "$done_string"

 dependency_string=$(echo $done_string|sed -r 's/ /:/g')

 echo "$dependency_string"

#####################collect_DE

module load R/4.1.0

Rscript_collect_DE=$(echo "$Rscripts_path""243_Per_cluster_Collect_DE_results_DESeq2_Seurat.R")

type=$(echo "collect_DE")


outfile_collect_DE=$(echo "$Log_files""outfile_3_""$type"".log")
touch $outfile_collect_DE
echo -n "" > $outfile_collect_DE
 name_collect_DE=$(echo "$type""_job")

cluster_array=$(echo '1,2,3,4,5,6,7,8,9,10,11,12,13')
comparison_array=$(echo 'Genotype_A.G_vs_G.G,Genotype_A.A_vs_G.G,Genotype_Del16_vs_G.G,Genotype_Del80_vs_G.G')
seff_name=$(echo "seff""_""$type")

#$dependency_string

myjobid_collect_DE=$(sbatch $dependency_string --job-name=$name_collect_DE --output=$outfile_collect_DE --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_collect_DE --cluster_array $cluster_array --comparison_array $comparison_array --type $type --out $output_dir")
myjobid_seff_collect_DE=$(sbatch --dependency=afterany:$myjobid_collect_DE --open-mode=append --output=$outfile_collect_DE --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_DE >> $outfile_collect_DE")



#####################collect_ORA

module load R/4.1.0

Rscript_collect_ORA=$(echo "$Rscripts_path""244_Per_cluster_Collect_ORA_DESeq2_Seurat.R")

type=$(echo "collect_ORA")


outfile_collect_ORA=$(echo "$Log_files""outfile_4_""$type"".log")
touch $outfile_collect_ORA
echo -n "" > $outfile_collect_ORA
name_collect_ORA=$(echo "$type""_job")

cluster_array=$(echo '1,2,3,4,5,6,7,8,9,10,11,12,13')
comparison_array=$(echo 'Genotype_A.G_vs_G.G,Genotype_A.A_vs_G.G,Genotype_Del16_vs_G.G,Genotype_Del80_vs_G.G')
seff_name=$(echo "seff""_""$type")

#$dependency_string

myjobid_collect_ORA=$(sbatch $dependency_string --job-name=$name_collect_ORA --output=$outfile_collect_ORA --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_collect_ORA --cluster_array $cluster_array --comparison_array $comparison_array --type $type --out $output_dir")
myjobid_seff_collect_ORA=$(sbatch --dependency=afterany:$myjobid_collect_ORA --open-mode=append --output=$outfile_collect_ORA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_ORA >> $outfile_collect_ORA")

####Classify_ORA_pathway_genes

module load R/4.1.0

Rscript_Classify_ORA_pathway_genes=$(echo "$Rscripts_path""254_Per_cluster_ORA_selecting_and_classifying_genes_time_covariate_analysis.R")

type=$(echo "Classify_ORA_pathway_genes")


outfile_Classify_ORA_pathway_genes=$(echo "$Log_files""outfile_5_""$type""_""$analysis"".log")
touch $outfile_Classify_ORA_pathway_genes
echo -n "" > $outfile_Classify_ORA_pathway_genes
name_Classify_ORA_pathway_genes=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")


ORA_ActivePathways_results=$(echo "$output_dir""ORA_ActivePathways_results.rds")

selected_pathways=$(echo 'GSE15330_LYMPHOID_MULTIPOTENT_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_UP,GSE15330_LYMPHOID_MULTIPOTENT_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_IKAROS_KO_UP,GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_UP,GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_DN,ZHENG_CORD_BLOOD_C3_MEGAKARYOCYTE_ERYTHROID_PROGENITOR,ZHENG_CORD_BLOOD_C1_PUTATIVE_MEGAKARYOCYTE_PROGENITOR,GSE15330_MEGAKARYOCYTE_ERYTHROID_VS_GRANULOCYTE_MONOCYTE_PROGENITOR_IKAROS_KO_DN,GSE15330_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_VS_PRO_BCELL_UP,GSE15330_HSC_VS_MEGAKARYOCYTE_ERYTHROID_PROGENITOR_IKAROS_KO_DN,DESCARTES_FETAL_SPLEEN_MEGAKARYOCYTES,DESCARTES_FETAL_LUNG_MEGAKARYOCYTES,DESCARTES_FETAL_KIDNEY_MEGAKARYOCYTES,DESCARTES_FETAL_HEART_MEGAKARYOCYTES,DESCARTES_FETAL_ADRENAL_MEGAKARYOCYTES,DESCARTES_FETAL_MUSCLE_MEGAKARYOCYTES,DESCARTES_FETAL_LIVER_MEGAKARYOCYTES,DESCARTES_FETAL_CEREBRUM_MEGAKARYOCYTES,DESCARTES_MAIN_FETAL_MEGAKARYOCYTES,DESCARTES_FETAL_PLACENTA_MEGAKARYOCYTES,REACTOME_FACTORS_INVOLVED_IN_MEGAKARYOCYTE_DEVELOPMENT_AND_PLATELET_PRODUCTION,GNATENKO_PLATELET_SIGNATURE,HP_IMPAIRED_COLLAGEN_INDUCED_PLATELET_AGGREGATION,HP_IMPAIRED_PLATELET_AGGREGATION,HP_ABNORMAL_PLATELET_FUNCTION,GAVISH_3CA_MALIGNANT_METAPROGRAM_34_PLATELET_ACTIVATION,RAGHAVACHARI_PLATELET_SPECIFIC_GENES,GOBP_PLATELET_ACTIVATION,HP_ABNORMAL_PLATELET_VOLUME,HP_ABNORMAL_PLATELET_VOLUME,GOBP_MICROTUBULE_CYTOSKELETON_ORGANIZATION_INVOLVED_IN_MITOSIS,REICHERT_MITOSIS_LIN9_TARGETS,CUX1_TARGET_GENES,GSE40273_GATA1_KO_VS_WT_TREG_DN,GATA2_01,GATA1_02,HUANG_GATA2_TARGETS_DN,GSE40277_EOS_AND_LEF1_TRANSDUCED_VS_GATA1_AND_SATB1_TRANSDUCED_CD4_TCELL_DN,GSE40274_GATA1_VS_FOXP3_AND_GATA1_TRANSDUCED_ACTIVATED_CD4_TCELL_UP,GATA1_05,GATA1_03,GATA1_01,RACCACAR_AML_Q6,REACTOME_RUNX1_REGULATES_GENES_INVOLVED_IN_MEGAKARYOCYTE_DIFFERENTIATION_AND_PLATELET_FUNCTION,MORF_RUNX1,AML1_Q6,AML1_01,AML_Q6,HP_ABNORMAL_MYELOID_LEUKOCYTE_MORPHOLOGY,GSE23502_WT_VS_HDC_KO_MYELOID_DERIVED_SUPPRESSOR_CELL_COLON_TUMOR_UP,GSE10325_CD4_TCELL_VS_MYELOID_UP,RUBENSTEIN_SKELETAL_MUSCLE_MYELOID_CELLS,NAKAYA_MYELOID_DENDRITIC_CELL_FLUMIST_AGE_18_50YO_7DY_UP,GSE10325_LUPUS_CD4_TCELL_VS_LUPUS_MYELOID_UP,GSE10325_LUPUS_CD4_TCELL_VS_LUPUS_MYELOID_DN,VALK_AML_CLUSTER_15,VALK_AML_CLUSTER_3,ROSS_ACUTE_MYELOID_LEUKEMIA_CBF,ALCALAY_AML_BY_NPM1_LOCALIZATION_DN,AKT_UP_MTOR_DN.V1_UP,AKT_UP.V1_DN,HALLMARK_PI3K_AKT_MTOR_SIGNALING,XU_HGF_SIGNALING_NOT_VIA_AKT1_6HR')

# --dependency=afterany:$myjobid_collect_ORA

myjobid_Classify_ORA_pathway_genes=$(sbatch --dependency=afterany:$myjobid_collect_ORA --job-name=$name_Classify_ORA_pathway_genes --output=$outfile_Classify_ORA_pathway_genes --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Classify_ORA_pathway_genes --ORA_ActivePathways_results $ORA_ActivePathways_results --selected_pathways $selected_pathways --type $type --out $output_dir")
myjobid_seff_Classify_ORA_pathway_genes=$(sbatch --dependency=afterany:$myjobid_Classify_ORA_pathway_genes --open-mode=append --output=$outfile_Classify_ORA_pathway_genes --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Classify_ORA_pathway_genes >> $outfile_Classify_ORA_pathway_genes")

#####################collect_GeneEXP

module load R/4.1.0

Rscript_collect_GeneEXP=$(echo "$Rscripts_path""245_Per_cluster_Collect_GeneEXP_results_DESeq2_Seurat_time_covariate.R")

type=$(echo "collect_GeneEXP")


outfile_collect_GeneEXP=$(echo "$Log_files""outfile_6_""$type"".log")
touch $outfile_collect_GeneEXP
echo -n "" > $outfile_collect_GeneEXP
 name_collect_GeneEXP=$(echo "$type""_job")


cluster_array=$(echo '1,2,3,4,5,6,7,8,9,10,11,12,13')
Genotype_array=$(echo 'G.G,A.G,A.A,Del16,Del80')
time_point_array=$(echo 'h24,h48,h72,h96')
clone_line_array=$(echo 'chrGFP_WTA,chrGFP_WTB,chrGFP_WTC,chrGFP_HET,chrGFP_KI_13,chrGFP_KI_27,chrGFP_KI_29,chrGFP_Del_16bp,chrGFP_Del_233,chrGFP_Del_235,chrGFP_Del_287')
seff_name=$(echo "seff""_""$type")


#$dependency_string

myjobid_collect_GeneEXP=$(sbatch $dependency_string --job-name=$name_collect_GeneEXP --output=$outfile_collect_GeneEXP --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=8 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_collect_GeneEXP --cluster_array $cluster_array --time_point_array $time_point_array --Genotype_array $Genotype_array --clone_line_array $clone_line_array --type $type --out $output_dir")
myjobid_seff_collect_GeneEXP=$(sbatch --dependency=afterany:$myjobid_collect_GeneEXP --open-mode=append --output=$outfile_collect_GeneEXP --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_GeneEXP >> $outfile_collect_GeneEXP")


#####################Heatmaps

module load R/4.1.0

Rscript_Heatmaps=$(echo "$Rscripts_path""246_Per_cluster_pheatmap_DESeq2_Seurat.R")

type=$(echo "Heatmaps")


outfile_Heatmaps=$(echo "$Log_files""outfile_7_""$type"".log")
touch $outfile_Heatmaps
echo -n "" > $outfile_Heatmaps
 name_Heatmaps=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

ery_genes_array=$(echo 'HBA1,HBG1,HBG2,MYL4,ANKRD26')
myeloid_genes_array=$(echo 'SLC44A1,AXL,INHBA,IGF1R,KIF21A,IL33,IFI16,ZEB1,ALDH1A2,IL6ST,TJP2,LTBP1,BACH2,CCL3,EPB41')
megak_genes_array=$(echo 'GP9,KIF22,CMTM5,KLF1,TUBB1,ITGA2B,FGF13,SPTA1,GFI1B,MYH9,MYLK,ITGB3,KCNT2,ADAM10,MEF2C,LRP12,ATP2C1,SPTB,TBL1X,BCOR,TBXA2R,FCGR2A,RHAG,GP1BB,EGF,GNAQ')
marker_genes_array=$(echo 'ANAPC7,STIL,KIF15,HBQ1,HBZ,HBA2,GP6,FYB1,CCR7,IL1B,IL33,GYPA,ITGA2B')
CUX1=$(echo 'CUX1')
ANAPC_genes=$(echo 'ANAPC1,ANAPC10,ANAPC11,ANAPC13,ANAPC15,ANAPC16,ANAPC2,ANAPC4,ANAPC5')
original_selection_of_genes=$(echo 'TUBB4B,NUP107,AURKA,KIF22,KNL1,CENPE,KIF4A,CENPF,BUB1B,FAM111A,DHFR,XRCC2,KIF23,NUP155,KIF11,KIF15,RACGAP1,HSPD1,CCT5,HSPB1,RPS14,ACTB,ACTG1,GLA,CDC42,UROD,CYCS,FLNA,IL1RL1,IFT80,ANGPT1,GATA1,DIAPH1,STAT5A,DOCK8,RUNX1,LMBR1,WIPF1,STAT1,IKZF1,NT5C3A,CEP152,UBE3C,ABL1,ANKRD26,SLC4A1,SMARCC1,SCAPER,ANK1,TGFB1,PUS7,GPHN,ZEB2,STAT5B,AGK,CUX1,STAT3,UBA1,CD109,CCR7,CCL3,ARFGEF2,IL23A,SH2B3,IL23R,CRADD,KIF2A,ACTN4,EPB41,TNIK,SLC44A1,AXL,INHBA,IGF1R,KIF21A,IL33,IFI16,ZEB1,ALDH1A2,IL6ST,TJP2,LTBP1,BACH2,ZBTB20,PLCG2,FN1,SOCS1,IGF1,KIT,B2M,SLC16A1,BANK1,ECM1,APOE,INHA,ENO3,MLH3,PSAP,FOLR1,CLCN3,OPTN,ACTN1,LAMP2,IL15,INPP5D,IREB2,PRKCH,FOXP1,NEDD4L,SCN3A,SCN2A,IL7R,PRDX3,GTPBP2,SYK,CALR,SOD1,SAMD9,ADCY6,ITM2B,TFRC,HBQ1,KLF1,TUBB1,ITGA2B,FGF13,SPTA1,GFI1B,HBZ,HBA1,MYL4,HBA2,HBG1,HBG2,MYH9,MYLK,ITGB3,KCNT2,ADAM10,MEF2C,LRP12,ATP2C1,SPTB,FYB1,TBL1X,BCOR,TBXA2R,GP6,FCGR2A,RHAG,GP1BB,EGF,GNAQ')

Selected_genes_classified=$(echo "$output_dir""Selected_genes_classified.rds")

GeneEXP_matrix=$(echo "$output_dir""Normalised_count_matrix.tsv")
metadata=$(echo "$output_dir""Metadata.rds")


# --dependency=afterany:$myjobid_collect_GeneEXP:$myjobid_Classify_ORA_pathway_genes

myjobid_Heatmaps=$(sbatch  --dependency=afterany:$myjobid_collect_GeneEXP:$myjobid_Classify_ORA_pathway_genes --job-name=$name_Heatmaps --output=$outfile_Heatmaps --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Heatmaps --ery_genes_array $ery_genes_array --megak_genes_array $megak_genes_array --myeloid_genes_array $myeloid_genes_array --marker_genes_array $marker_genes_array --CUX1 $CUX1 --GeneEXP_matrix $GeneEXP_matrix --ANAPC_genes $ANAPC_genes --original_selection_of_genes $original_selection_of_genes --metadata $metadata --Selected_genes_classified $Selected_genes_classified --type $type --out $output_dir")
myjobid_seff_Heatmaps=$(sbatch --dependency=afterany:$myjobid_Heatmaps --open-mode=append --output=$outfile_Heatmaps --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Heatmaps >> $outfile_Heatmaps")


#####################volcano_plots

module load R/4.1.0

Rscript_volcano_plots=$(echo "$Rscripts_path""247_Per_cluster_Volcano_plot_Seurat_DESeq2.R")

type=$(echo "volcano_plots")


outfile_volcano_plots=$(echo "$Log_files""outfile_8_""$type"".log")
touch $outfile_volcano_plots
echo -n "" > $outfile_volcano_plots
name_volcano_plots=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")


#highlighted_genes=$(echo 'CUX1,RUNX1,UBE3C,TBL1X,ALDH1A2,BANK1,IKZF2,SPI1,GATA1,ITGA2B,FYB1,GFI1B,CENPE,FLNA,NUP62,KIF4A,GJA1,WDR62,RACGAP1,HNRNPU,KIF2A,KIF11,MECP2,PCNT,CENPJ,STIL,AURKA,AURKC,TPR,SMC1A,KIF23,CDC20,GP6,LCP2,EZH2,XRCC2,CDKN1A,CDKN2D,KIF15,CENPE,CENPF,AURKB,TUBB6,TUBB4B,TUBB1,TBXA2R,TBXAS1,PDGFA,CCR7,BCOR,IL7R,IGF1R,LRMDA,IL4R,IL1A,IL1B,BACH2,WIPF1,WASL,KMT2E,GRID2,ABCA1,IL6ST,ZEB1,AARS2,ABCB8,ABCC5-AS1,ABCD4,ABCG1,ABCG2,ABHD16A,ABHD2,ACBD4,ACBD5,ACSS3,ACTL6A,ACYP2,ADARB1,ADGRB3,ADGRB3-DT,ADGRF1,ADGRF4,AGPAT1,AGR2,AHCYL2,AHI1,AIP,AKAP1,AKR1E2,ALDH1A2,ALDH3A2,ALDH3B2,ALOXE3,AMMECR1,ANAPC10,ANKHD1,ANKHD1-DT,ANKHD1-EIF4EBP3,ANKRD13A,ANKRD34A,ANO6,ANXA7,AP1M1,AP1M2,AP3M2,AP3S2,AP5M1,APOC1,APOOL,ARAP1,ARFGEF2,ARHGAP32,ARHGEF28,ARHGEF37,ARID4A,ARID5B,ARL14EPP1,ARMH4,ARMT1,ARPC5,ARRDC3,ASAH2B,ASB3,ASB8,ASH2L,ASMTL,ASNSD1,ATF1,ATF7IP,ATF7IP2,ATG101,ATG12,ATG13,ATG4B,ATG5,ATMIN,ATP10B,ATP5MC1,ATP6V1G2,ATP6V1G2-DDX39B,ATPAF1,ATPSCKMT,ATXN7L3B,AURKA,BABAM1,BAHCC1,BAIAP2L1,BATF,BAZ2A,BBOX1-AS1,BCAS4,BCKDHA,BETALINC1,BICDL1,BLOC1S2,BLOC1S6,BLTP1,BMAL1,BMF,BNIP2,BOD1L1,BRD8,BRI3BP,BRINP3,BRINP3-DT,BTN3A2,BTNL8,BUB1B,C11orf21,C12orf76,C19orf38,C1QA,C2CD3,C2CD5,C4orf46,C4orf47,CA1,CAB39L,CABIN1,CACNA1A,CAPN8,CAPS2,CARD8-AS1,CARF,CASC11,CCDC124,CCDC137,CCDC159,CCDC186,CCDC88C,CCN2,CCNP,CCNQ,CCT5,CCT7P2,CCT8,CD36,CDC73,CDH11,CDK2AP1,CDV3,CEACAM19,CEBPG,CENPN-AS1,CEP112,CEP152,CEP290,CEP350,CEP78,CFLAR-AS1,CHCHD3P1,CHCHD5,CHD6,CHD9NB,CIAO2A,CIDEB,CISTR,CLCN3,CLIP1,CLIP4,CLK3,CLK4,CLN8,CLNK,CLPX,CLSTN3,CLTC,CMTM3,CNOT1,CNOT6L,CNPY1,COG3,COG7,COG8,COL17A1,COL4A2,COMMD2,COPB2,COPE,COPS4,COQ8B,COX15,COX16,COX17,COX7C,CPED1,CPNE2,CRADD,CREB1,CREB3L1,CRLF2,CRPPA-AS1,CRYGS,CSNK1A1,CSNK1G1,CSTF1,CUTC,CUX1,CXXC1,CYP2B7P,DACT3,DACT3-AS1,DAGLB,DARS1-AS1,DBTP1,DCC,DCTN1,DDX49,DEDD,DET1,DHRS4-AS1,DHRSX,DIXDC1,DLGAP1-AS2,DLX2,DMAP1,DNAH7,DNAJB2,DNAJC1,DNAJC16,DOLPP1,DOP1A,DRAIC,DRG1,DTNA,DTWD1,DTX4,DUSP6,DYNC2H1,EDEM2,EFCAB14-AS1,EFL1,EFTUD2,EGR2,EHD4,EIF2AK2,EIF2B1,EIF4E,EIF4G3,ELK2AP,ELP3,EML6,EN1,ENO3,ENSG00000200999,ENSG00000222095,ENSG00000230615,ENSG00000254363,ENSG00000255314,ENSG00000270571,EOGT,EPB41L4B,EPHA7,ERAP1,ERCC2,ERGIC1,ERLEC1,ESR1,ETFDH,EVI5L,EXOGP1,EXOSC5,EXOSC8,EXPH5,EXTL3,EYA3,FABP5,FAM111A,FAM114A2,FAM13A,FAM217B,FAM227B,FAM3B,FAM47E,FARP2,FAXDC2,FBRS,FBXL18,FBXO24,FBXO34,FBXO34-AS1,FBXO46,FCHSD1,FDX2,FEZ2,FMC1,FMC1-LUC7L2,FMO1,FNBP4,FNIP2,FOXJ3,FOXN3-AS1,FREM2,FREM2-AS1,FSTL4,FTCDNL1,FTSJ1,FXYD6-AS1,GALNT2,GARRE1,GATA3,GATA3-AS1,GBE1,GCDH,GCN1,GCNT3,GDE1,GFM2,GHET1,GIRGL,GLIPR1L2,GMDS,GNG4,GOLT1B,GPBP1L1,GPCPD1,GPX2,GREB1,GRPEL2,GSE1,GSTA4,GSTP1,GTF2B,GTF2H3,GUSBP1,GUSBP2,GXYLT1,H2BC15,HABP2,HACD2,HARBI1,HAUS5,HAUS5-DT,HAVCR2,HCG20,HDAC8,HEXIM2,HEXIM2-AS1,HM13,HMGB1,HNF4A,HNRNPA1P42,HNRNPA3,HOXA-AS3,HOXB7,HRG-AS1,HSD17B2,HSDL1,HSPD1,HSPH1,IBTK,ICE2,ID2,ID2-AS1,IFI16,IFRD1,IFT46,IGLV3-32,IKBKB-DT,IKZF2,IL1R1,ILF2,IMPDH1,INO80,INO80B,INO80B-WBP1,INSM2,INTS13,INTS14,INTS8,IQCG,ISCA2P1,ISY1,ISY1-RAB43,ITFG2,ITFG2-AS1,ITPKC,JPT1,JTB,JTB-DT,KALRN,KAT5,KATNB1,KBTBD2,KCNJ15,KCNS2,KCTD3,KDM5A,KDM8,KIAA0319,KIAA0513,KIAA1217,KIAA2026,KIF20A,KIF2A,KLHDC9,KLHL21,KLHL38,KNL1,KPNA6,KRT18P46,KRT8,KRTAP3-1,KYAT3,KYNU,LAPTM4B,LARP1B,LARP7,LARS1,LATS1,LEKR1,LENG1,LEPROTL1,LIMA1,LINC00466,LINC00992,LINC01132,LINC01354,LINC01579,LINC01701,LINC01732,LINC01734,LINC01823,LINC01972,LINC02015,LINC02243,LINC02252,LINC02343,LINC02533,LINC02707,LINC02842,LINC02934,LINC02950,LINC02960,LINC02984,LMO4,LNCATV,LNMICC,LOC124900473,LOC124900883,LPAR2,LPP,LPXN,LRBA,LRCH4,LRIG2,LRRC23,LRRC37A3,LSM10,LSM14A,LSM14B,LTB4R2,LURAP1L-AS1,LYSETP1,LYZ,MACC1,MAFA,MALAT1,MAN2A1,MAN2A1-DT,MAN2A2,MAP4,MARCHF10,MARF1,MAX,MBP,MCC,MCRIP1,MCTS1,MED16,MED7,MEF2C,MEF2C-AS1,MEIS1,MEIS1-AS3,MFAP3,MGST3,MIOS,MIR1265,MIR200CHG,MIR3646,MIR3926-1,MIR4659B,MIR4766,MIR5091,MIR5188,MIR548AQ,MIR6070,MIR644A,MIR7849,MKRN2OS,MLLT3,MLPH,MMP11,MON1B,MPND,MPV17L2,MROH8,MRPS23,MRPS27,MRPS33,MSI2,MTFR1,MTMR11,MTO1,MYADM,MYADM-AS1,MYC,MYCL,MYLK-AS1,MYO3A,MYT1,NADK2,NAGLU,NAIF1,NCOA7,NDC1,NFE2,NFKBIL1,NHEJ1,NHSL1,NINJ2,NIP7,NIPSNAP1,NME2P2,NOL6,NOSIP,NR2F1-AS1,NR2F2,NSA2,NSL1,NT5C3A,NUP107,NUP155,NUSAP1,OAT,ODAD3,OGDH,OIP5,OIP5-AS1,OPLAH,OR10J2P,OR1AA1P,OTX1,OXSM,P2RY6,P4HB,PANK2,PARD6B,PARP2,PATJ,PATZ1,PCBP1-AS1,PCBP2,PCDHB3,PCLAF,PCM1,PCTP,PDCD6,PEX7,PGP,PHACTR3,PIERCE2,PIGO,PIGO-AS1,PIH1D1,PIK3C2B,PIK3CA,PIP5K1C,PISD,PJA2,PKD1L2,PKM,PLA2G12A,PLA2G4E-AS1,PLAAT3,PLAC8,PLBD1,PLEC,PLEKHA5,PLEKHA8P1,PLEKHB1,PLEKHH1,PLEKHM1,PLIN5,PLSCR4,PM20D2,PMEL,PMVK,PNPLA7,POC1B,POC1B-GALNT4,POLR1G,POLR2D,POLR2G,POLR3G,POU2AF1,PPP1R13L,PPP1R3D,PPP2CA,PPP2R5B,PPP4R1L,PPP6R1,PPP6R3,PPRC1,PRDM2,PREP,PRKACA,PRKAR1A,PRKCI,PRKCSH,PRORSD1P,PROSER1,PRR13,PRR5,PSD4,PSKH2,PSMA3-AS1,PSMB1,PSMB3,PSMD10P2,PSMD9,PSME2P3,PSME3,PSMG3,PTCH1,PTPRF,QPCTL,R3HDM2,R3HDML-AS1,RAB34,RAB5B,RABGGTB,RAD51B,RAD9B,RANGRF,RARG,RBBP4,RBBP5,RBBP8,RBM27,RBMXL1,RECQL,REG4,RGS17P1,RGS6,RIF1,RIMKLB,RIN3,RMND1,RN7SK,RN7SKP192,RN7SL346P,RN7SL39P,RN7SL445P,RNA5SP323,RNF44,RNU1-132P,RNU4-1,RNU4-2,RNU5B-1,RNU5D-1,RNU6-1158P,RNU6-433P,RNU6-821P,RNU6-9,RNU7-195P,RNY1,RNY3,ROCK1P1,RPL19P14,RPL21P12,RPL21P131,RPL23A,RPL27A,RPL30P11,RPL35A,RPL36,RPL39P40,RPL3P4,RPL41,RPL7P30,RPL7P41,RPN2,RPP21,RPS11,RPS18,RPS29,RPS6KA5,RPS6KB1,RPS8,RSL24D1,RSRP1,RUNX1,SACM1L,SAE1,SAFB,SAFB2,SAMD1,SAMD4B,SAMD9L,SAPCD2P2,SAR1B,SART3,SCAMP1,SCAMP5,SCN2A,SEC24C,SEMA4B,SEMA4G,SENP2,SEPHS2,SERPINB9P1,SETDB2,SFR1,SFT2D2,SGPP1,SHLD1,SLAIN2,SLC11A2,SLC15A2,SLC1A3,SLC22A17,SLC22A5,SLC24A1,SLC25A25,SLC2A4RG,SLC33A1,SLC35E1,SLC35E3,SLC38A2,SLC4A1AP,SLC9A6,SLTM,SMAD1,SMARCD2,SMIM10L2B,SMIM13,SNAI3-AS1,SND1-DT,SNHG1,SNIP1,SNORD111B,SNORD113-9,SNORD118,SNORD13,SNORD22,SNORD28,SNORD30,SNORD55,SNRNP35,SNRNP70,SNRPD2,SNRPE,SNX1,SNX12,SNX8,SOX2-OT,SPAG9,SPAST,SPATA24,SPATS2L,SPECC1P1,SPICE1,SPINK5,SPRED2,SPRYD4,SPTA1,SREK1,SRI,SSR1,SSR4P1,STAP2,STARD10,STPG2,STRIP1,STT3A,STUM,STX16,STX16-NPEPL1,STYK1,SUMO2,SUPT7L,SYNE1,SYS1,SYTL3,TAB3,TAF4,TANGO6,TATDN1P1,TATDN3,TBC1D30,TBL1X,TBP,TBX3,TBX6,TCF4,TCF7L2,TECPR1,TEDC1,TENT2,TESK2,TGFBI,THBS1,THBS4-AS1,TIPARP,TLE4,TLE6,TMC1,TMCO6,TMED1,TMEM125,TMEM161A,TMEM161B,TMEM161B-DT,TMEM177,TMEM214,TMEM238L,TMEM243,TMEM62,TMF1,TMOD3,TMPRSS11F,TMTC3,TNPO1,TOGARAM2,TOMM22P6,TPD52,TPM4,TRAV33,TRIM54,TRIP4,TRMT12,TSHB,TSPAN31,TTC39C,TTI2,TTLL6,TUBB4B,TUBD1,TXN2,TXNP5,UBALD2,UBC,UBE2D2,UBE2O,UBE3C,UBQLN1,UBR3,UCA1,UFSP2,ULBP1,UMODL1,UQCC1,USPL1,UTP11,UTP3,UTS2B,VAC14,VAMP1,VCP,VCPIP1,VDAC2,VEPH1,VMP1,VN1R28P,VPS29,VPS36,VPS4A,VPS52,VWA5A,VWA8,WDPCP,WDR45B,WDR83,WDR83OS,WSB1,YAE1,YAP1,YY1,ZBTB8OS,ZDHHC1,ZFAND3,ZFAND3-DT,ZFAND6,ZFHX2,ZFHX3,ZGRF1,ZIC3,ZIM2-AS1,ZNF213,ZNF213-AS1,ZNF234,ZNF252P,ZNF263,ZNF282,ZNF398,ZNF410,ZNF423,ZNF446,ZNF585B,ZNF609,ZNF623,ZNF689,ZNF713,ZNF731P,ZNF74,ZNHIT3,ZSCAN31')


platelet_volume_genes=$(echo 'ITGA2B,FYB1,GFI1B,TUBB1,FLI1,FLNA,DIAPH1,ITGB3')
platelet_genes=$(echo 'GP6,TBXAS1,TBXA2R')
EZH2_signature=$(echo 'EZH2,XRCC2,CDKN1A,CDKN2D,ZEB1,ZEB2,TWIST')


AKT_signature=$(echo 'DHCR24,SLC44A1,GUSB,SLC6A6,EDARADD,SASH1,FUT8,ITPR2,MAMLD1,IL6ST,MOCOS,ZEB1,PDE10A,SLC44A1,MSR1,PLCB1,SMAD2')

#EZH2_signature=$(echo 'ABCA5,ABI3BP,ABRAXAS2,ACAD8,ACKR2,ADCY7,AHCTF1,AKAP13,ANO1,ANOS1,ANTXR2,ANXA10,AOX1,AP3B2,APBB2,APP,ARHGAP5,ARHGEF3,ATP8B1,ATP9A,AVL9,BBOF1,BCL2A1,C18orf25,C1D,C1S,C6orf89,CACTIN,CAP2,CCDC50,CCL20,CCM2L,CCN2,CD164,CD59,CFDP1,CMTR2,CNIH3,CNTNAP3,COA3,COL1A1,COL3A1,COL8A1,COMMD10,CP,CPE,CRISPLD1,CSGALNACT2,CXCL5,CYP19A1,DLX2,DNAJB8-AS1,DNER,DOCK4,DPF3,DSEL,DYNLT2,EDIL3,EDN1,EGFR,ELMO1,ELOVL2,EMP1,EPHA5,FABP4,FAM43A,FAR2,FHL1,FLI1,FOXA2,GABARAPL2,GAL3ST1,GALNS,GBA2,GJA1,GLP1R,GNA14,GPC6,GRK4,GULP1,H4C8,HAS2,HAS3,HEMGN,HLF,HOXC8,HPGD,HSD17B2,ID2,IFI6,IFIT2,IFIT3,IGFBP3,IGFBP5,IL13RA1,IL6R,IL9R,IQCA1,IRX5,ITGA6,ITGBL1,KAAG1,KCNMA1,KIAA1549L,KLF2,KLF3,KLF9,KRTAP4-12,LINC00550,LMAN1,LMCD1,LOXL4,LPAR6,LY9,MAB21L2,MAGEC2,MAGT1,MAP3K8,MAP4K3,MATN3,MCPH1,MEF2D,MGLL,MMP1,MS4A1,MSX2,MYEF2,MYLK,MYO16,MYO1D,MYRFL,NAV2,NCF4,NDST2,NEBL,NEGR1,NEK3,NIP7,NIPSNAP3A,NMNAT1,NNMT,NOL3,NPAS2,NRP1,NT5E,NTN4,OR12D2,OR7D2,OSR1,PCK1,PCLO,PIERCE1,PLA2G12B,PLAAT2,PPP1R3B,PRB3,PRIM2,PRKAR1B,PRR16,PRSS35,PYROXD1,RAB3B,RBMS3,RBMX,RBP4,RNASE4,RNF144B,RNF213-AS1,RNPC3,RRAGD,SAMD4A,SCARB1,SCN2A,SERPINB3,SERPINE1,SGK1,SHISA3,SIGLEC5,SLAMF6,SLC12A2,SLC2A10,SLC4A8,SNORC,SNTB2,SORL1,SPARC,SPART-AS1,SPOCK1,SPRY4,SRRM2,SSPN,ST3GAL1,ST3GAL6,ST6GALNAC4,SUSD4,SYT6,TACR1,TFPI,TMEM116,TMEM45A,TMSB15A,TMTC1,TNFAIP6,TNFRSF6B,TNP2,TNRC6B,TPM3,TRIM2,TSC1,TSPAN8,TTC28,TYRP1,UBXN10,UHRF2,USP34,UST,UTS2,VAMP4,VPS13C,VTI1A,WWTR1,YPEL5,ZEB1,ZFX,ZNF175,ZNF37BP,ZNF567,ZNF76,ZNF766')

#ZEB1_signature=$(echo 'CD24,CDH1,CDH11,CDH3,CDH4,CLDN7,CRB3,CXADR,DMKN,DSC2,DSP,EPCAM,EPPK1,F11R,GJB2,GJB3,MAL2,MARVELD2,MPZL2,MUC1,OCLN,PATJ,PCDH7,PKP3,PMEPA1,PPL,SCEL,SFN,SH3YL1,SHROOM3,SYTL1,TACSTD2,TMEM30B,TSPAN1,TSPAN15')
ZEB1_signature=$(echo 'CD24,CDH1')
#AKT_signature=$(echo 'CUX1,ADAMTSL4-AS1,AKT1,AKT2,AKT3,ANGPT1,ANGPT2,ANGPT4,ATF2,ATF4,ATF6B,BAD,BCL2,BCL2L1,BCL2L11,BDNF,BRCA1,CASP9,CCND1,CCND2,CCND3,CCNE1,CCNE2,CD19,CDC37,CDK2,CDK4,CDK6,CDKN1A,CDKN1B,CHAD,CHRM1,CHRM2,CHUK,COL1A1,COL1A2,COL2A1,COL4A1,COL4A2,COL4A3,COL4A4,COL4A5,COL4A6,COL6A1,COL6A2,COL6A3,COL6A5,COL6A6,COL9A1,COL9A2,COL9A3,COMP,CREB1,CREB3,CREB3L1,CREB3L2,CREB3L3,CREB3L4,CREB5,CSF1,CSF1R,CSF3,CSF3R,CSH1,CSH2,DDIT4,EFNA1,EFNA2,EFNA3,EFNA4,EFNA5,EGF,EGFR,EIF4B,EIF4E,EIF4E1B,EIF4E2,EIF4EBP1,EPHA2,EPO,EPOR,F2R,FASLG,FGF1,FGF10,FGF11,FGF12,FGF13,FGF14,FGF17,FGF18,FGF19,FGF2,FGF20,FGF21,FGF22,FGF23,FGF3,FGF4,FGF5,FGF6,FGF7,FGF8,FGF9,FGFR1,FGFR2,FGFR3,FGFR4,FLT1,FLT3,FLT3LG,FLT4,FN1,FOXO3,G6PC1,G6PC2,G6PC3,GH1,GH2,GHR,GNB1,GNB2,GNB3,GNB4,GNB5,GNG10,GNG11,GNG12,GNG13,GNG2,GNG3,GNG4,GNG5,GNG7,GNG8,GNGT1,GNGT2,GRB2,GSK3B,GYS1,GYS2,HGF,HRAS,HSP90AA1,HSP90AB1,HSP90B1,IBSP,IFNA1,IFNA10,IFNA13,IFNA14,IFNA16,IFNA17,IFNA2,IFNA21,IFNA4,IFNA5,IFNA6,IFNA7,IFNA8,IFNAR1,IFNAR2,IFNB1,IGF1,IGF1R,IGF2,IKBKB,IKBKG,IL2,IL2RA,IL2RB,IL2RG,IL3,IL3RA,IL4,IL4R,IL6,IL6R,IL7,IL7R,INS,INSR,IRS1,ITGA1,ITGA10,ITGA11,ITGA2,ITGA2B,ITGA3,ITGA4,ITGA5,ITGA6,ITGA7,ITGA8,ITGA9,ITGAV,ITGB1,ITGB3,ITGB4,ITGB5,ITGB6,ITGB7,ITGB8,JAK1,JAK2,JAK3,KDR,KIT,KITLG,KRAS,LAMA1,LAMA2,LAMA3,LAMA4,LAMA5,LAMB1,LAMB2,LAMB3,LAMB4,LAMC1,LAMC2,LAMC3,LPAR1,LPAR2,LPAR3,LPAR4,LPAR5,LPAR6,MAP2K1,MAP2K2,MAPK1,MAPK3,MCL1,MDM2,MET,MLST8,MTOR,MYB,MYC,NFKB1,NGF,NGFR,NOS3,NRAS,NTF3,NTF4,NTRK1,NTRK2,OSM,OSMR,PCK1,PCK2,PDGFA,PDGFB,PDGFC,PDGFD,PDGFRA,PDGFRB,PDPK1,PGF,PHLPP1,PHLPP2,PIK3AP1,PIK3CA,PIK3CB,PIK3CD,PIK3CG,PIK3R1,PIK3R2,PIK3R3,PIK3R5,PIK3R6,PKN1,PKN2,PKN3,PPP2CA,PPP2CB,PPP2R1A,PPP2R1B,PPP2R2A,PPP2R2B,PPP2R2C,PPP2R2D,PPP2R3A,PPP2R3B,PPP2R3C,PPP2R5A,PPP2R5B,PPP2R5C,PPP2R5D,PPP2R5E,PRKAA1,PRKAA2,PRKCA,PRL,PRLR,PTEN,PTK2,RAC1,RAF1,RBL2,RELA,RELN,RHEB,RPS6,RPS6KB1,RPS6KB2,RPTOR,SGK1,SGK2,SGK3,SOS1,SOS2,SPP1,STK11,SYK,TCL1A,TCL1B,TEK,TGFA,THBS1,THBS2,THBS3,THBS4,THEM4,TLR2,TLR4,TNC,TNN,TNR,TNXB,TP53,TSC1,TSC2,VEGFA,VEGFB,VEGFC,VEGFD,VTN,VWF')
selected_clusters=$(echo '3_1')
selected_comparisons=$(echo 'Genotype_A_G_vs_G_G,Genotype_A_A_vs_G_G,Genotype_Del16_vs_G_G,Genotype_Del80_vs_G_G')
DE_results=$(echo "$output_dir""DE_results.rds")



# --dependency=afterany:$myjobid_collect_GeneEXP:$myjobid_Classify_ORA_pathway_genes

myjobid_volcano_plots=$(sbatch --dependency=afterany:$myjobid_collect_GeneEXP:$myjobid_Classify_ORA_pathway_genes --job-name=$name_volcano_plots --output=$outfile_volcano_plots --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_volcano_plots --platelet_volume_genes $platelet_volume_genes --platelet_genes $platelet_genes --EZH2_signature $EZH2_signature --ZEB1_signature $ZEB1_signature --AKT_signature $AKT_signature --selected_clusters $selected_clusters --selected_comparisons $selected_comparisons --DE_results $DE_results --type $type --out $output_dir")
myjobid_seff_volcano_plots=$(sbatch --dependency=afterany:$myjobid_volcano_plots --open-mode=append --output=$outfile_volcano_plots --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_volcano_plots >> $outfile_volcano_plots")
