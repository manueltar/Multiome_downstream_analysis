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

#  rm -rf $Log_files
#  mkdir -p $Log_files

 

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

     ### DA_test

 
     DA_route=$(echo "$Master_path_analysis""$cluster_array_sel""/")
     echo "$DA_route"

     rm -rf $DA_route
     mkdir -p $DA_route
     
     type=$(echo "$cluster_array_sel""_""DA_test")
     outfile_DA_test=$(echo "$Log_files""outfile_1_""$type"".log")
     touch $outfile_DA_test
     echo -n "" > $outfile_DA_test
     name_DA_test=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")



     Rscript_DA_test=$(echo "$Rscripts_path""248_Per_cluster_Seurat_to_DESeq2_multiome_ATAC_time_covariate.R")

     Seurat_object=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/merged_downstream_analysis_FOR_DANIELA.rds")
     Master_peak_file_with_SNP=$(echo '/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Peak_extraction_and_classification_for_time_covariate_model/Master_peak_file_with_SNP.rds')
     
     myjobid_DA_test=$(sbatch --job-name $name_DA_test --output=$outfile_DA_test --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=6 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_DA_test --Seurat_object $Seurat_object --Master_peak_file_with_SNP $Master_peak_file_with_SNP --seurat_cluster $cluster_array_sel --type $type --out $DA_route")
     myjobid_seff_DA_test=$(sbatch --dependency=afterany:$myjobid_DA_test --open-mode=append --output=$outfile_DA_test --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_DA_test >> $outfile_DA_test")

    
     echo "->>>$myjobid_DA_test"
     arr[${#arr[@]}]="$myjobid_DA_test"


 done
      


conda deactivate

 done_string=$(echo "--dependency=afterany:"""""${arr[@]}"""")
 echo "$done_string"

 dependency_string=$(echo $done_string|sed -r 's/ /:/g')

 echo "$dependency_string"

#####################collect_DA

module load R/4.1.0

Rscript_collect_DA=$(echo "$Rscripts_path""249_Per_cluster_Collect_DA_results_DESeq2_Seurat.R")

type=$(echo "collect_DA")


outfile_collect_DA=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_collect_DA
echo -n "" > $outfile_collect_DA
name_collect_DA=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

cluster_array=$(echo '1,2,3,4,5,6,7,8,9,10,11,12,13')
comparison_array=$(echo 'Genotype_A.G_vs_G.G,Genotype_A.A_vs_G.G,Genotype_Del16_vs_G.G,Genotype_Del80_vs_G.G')
Master_peak_file_with_SNP_numbered_peaks=$(echo '/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Peak_extraction_and_classification_for_time_covariate_model/Master_peak_file_with_SNP_numbered_peaks.rds')
#$dependency_string

myjobid_collect_DA=$(sbatch $dependency_string --job-name=$name_collect_DA $dependency_string --output=$outfile_collect_DA --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_collect_DA --Master_peak_file_with_SNP_numbered_peaks $Master_peak_file_with_SNP_numbered_peaks --cluster_array $cluster_array --comparison_array $comparison_array --type $type --out $output_dir")
myjobid_seff_collect_DA=$(sbatch --dependency=afterany:$myjobid_collect_DA --open-mode=append --output=$outfile_collect_DA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_DA >> $outfile_collect_DA")



#####################collect_ATAC_counts

module load R/4.1.0

Rscript_collect_ATAC_counts=$(echo "$Rscripts_path""250_Per_cluster_Collect_ATACcounts_results_DESeq2_Seurat_time_covariate.R")

type=$(echo "collect_ATAC_counts")


outfile_collect_ATAC_counts=$(echo "$Log_files""outfile_3_""$type"".log")
touch $outfile_collect_ATAC_counts
echo -n "" > $outfile_collect_ATAC_counts
name_collect_ATAC_counts=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

cluster_array=$(echo '1,2,3,4,5,6,7,8,9,10,11,12,13')
Genotype_array=$(echo 'G.G,A.G,A.A,Del16,Del80')
time_point_array=$(echo 'h24,h48,h72,h96')
clone_line_array=$(echo 'chrGFP_WTA,chrGFP_WTB,chrGFP_WTC,chrGFP_HET,chrGFP_KI_13,chrGFP_KI_27,chrGFP_KI_29,chrGFP_Del_16bp,chrGFP_Del_233,chrGFP_Del_235,chrGFP_Del_287')

#$dependency_string

myjobid_collect_ATAC_counts=$(sbatch $dependency_string --job-name=$name_collect_ATAC_counts --output=$outfile_collect_ATAC_counts --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=8 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_collect_ATAC_counts --cluster_array $cluster_array --time_point_array $time_point_array --Genotype_array $Genotype_array --clone_line_array $clone_line_array --type $type --out $output_dir")
myjobid_seff_collect_ATAC_counts=$(sbatch --dependency=afterany:$myjobid_collect_ATAC_counts --open-mode=append --output=$outfile_collect_ATAC_counts --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_collect_ATAC_counts >> $outfile_collect_ATAC_counts")



#####################Heatmaps

module load R/4.1.0

Rscript_Heatmaps=$(echo "$Rscripts_path""251_Per_cluster_pheatmap_DESeq2_Seurat_ATAC.R")

type=$(echo "Heatmaps")


outfile_Heatmaps=$(echo "$Log_files""outfile_4_""$type"".log")
touch $outfile_Heatmaps
echo -n "" > $outfile_Heatmaps
name_Heatmaps=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")


ery_genes_array=$(echo 'HBQ1,HBZ,HBA1,HBA2,HBG1,HBG2,MYL4,ANKRD26,')
myeloid_genes_array=$(echo 'CCR7,IL33,SLC44A1,AXL,INHBA,IGF1R,KIF21A,IL33,IFI16,ZEB1,ALDH1A2,IL6ST,TJP2,LTBP1,BACH2,CCL3,EPB41')
megak_genes_array=$(echo 'GP9,KIF15,KIF22,GP6,CMTM5,KLF1,TUBB1,ITGA2B,FGF13,SPTA1,GFI1B,MYH9,MYLK,ITGB3,KCNT2,ADAM10,MEF2C,LRP12,ATP2C1,SPTB,FYB1,TBL1X,BCOR,TBXA2R,GP6,FCGR2A,RHAG,GP1BB,EGF,GNAQ')
marker_genes_array=$(echo 'ANAPC7,STIL,KIF15,HBQ1,HBZ,HBA2,GP6,FYB1,CCR7,IL1B,IL33,GYPA,ITGA2B')
ANAPC_genes=$(echo 'ANAPC1,ANAPC10,ANAPC11,ANAPC13,ANAPC15,ANAPC16,ANAPC2,ANAPC4,ANAPC5,ANAPC7')
CUX1=$(echo 'CUX1')


Selected_genes_classified=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/DESeq2_time_covariate_model/""Selected_genes_classified.rds")

ATAC_matrix=$(echo "$output_dir""Normalised_count_matrix.tsv")
metadata=$(echo "$output_dir""Metadata.rds")


Master_peak_file_with_SNP_numbered_peaks=$(echo '/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Peak_extraction_and_classification_for_time_covariate_model/Master_peak_file_with_SNP_numbered_peaks.rds')

# --dependency=afterany:$myjobid_collect_ATAC_counts

myjobid_Heatmaps=$(sbatch --dependency=afterany:$myjobid_collect_ATAC_counts --output=$outfile_Heatmaps --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Heatmaps --ery_genes_array $ery_genes_array --megak_genes_array $megak_genes_array --myeloid_genes_array $myeloid_genes_array --marker_genes_array $marker_genes_array --ATAC_matrix $ATAC_matrix --CUX1 $CUX1 --ANAPC_genes $ANAPC_genes --Master_peak_file_with_SNP_numbered_peaks $Master_peak_file_with_SNP_numbered_peaks --metadata $metadata --Selected_genes_classified $Selected_genes_classified --type $type --out $output_dir")
myjobid_seff_Heatmaps=$(sbatch --dependency=afterany:$myjobid_Heatmaps --open-mode=append --output=$outfile_Heatmaps --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Heatmaps >> $outfile_Heatmaps")

#####################volcano_plots

module load R/4.1.0

Rscript_volcano_plots=$(echo "$Rscripts_path""252_Per_cluster_Volcano_plot_Seurat_DESeq2_DA.R")

type=$(echo "volcano_plots")


outfile_volcano_plots=$(echo "$Log_files""outfile_5_""$type"".log")
touch $outfile_volcano_plots
echo -n "" > $outfile_volcano_plots
name_volcano_plots=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

#highlighted_genes=$(echo 'CUX1,RUNX1,TBL1X,ALDH1A2,IKZF2,SPI1,GATA1,ITGA2B,FYB1,GFI1B,CENPE,FLNA,NUP62,KIF4A,GJA1,WDR62,RACGAP1,HNRNPU,KIF2A,KIF11,MECP2,PCNT,CENPJ,STIL,AURKA,AURKC,TPR,SMC1A,KIF23,CDC20,GP6,LCP2,EZH2,KIF15,CENPE,CENPF,AURKB,TUBB6,TUBB4B,TUBB1,TBXA2R,TBXAS1,PDGFA,CCR7,BCOR,IL7R,IGF1R,FUT8,LRMDA,IL4R,IL1A,IL1B,BACH2,WIPF1,WASL,KMT2E,GRID2,ABCA1,IL6ST,ZEB1,UBE3C,BANK1')

DA_results=$(echo "$output_dir""DA_results.rds")
platelet_volume_genes=$(echo 'ITGA2B,FYB1,GFI1B,TUBB1,FLI1,FLNA,DIAPH1,ITGB3')
platelet_genes=$(echo 'GP6,TBXAS1,TBXA2R')
EZH2_signature=$(echo 'EZH2,XRCC2,CDKN1A,CDKN2D,ZEB1,ZEB2,TWIST')
AKT_signature=$(echo 'DHCR24,SLC44A1,GUSB,SLC6A6,EDARADD,SASH1,FUT8,ITPR2,MAMLD1,IL6ST,MOCOS,ZEB1,PDE10A,SLC44A1,MSR1,PLCB1,SMAD2')
ZEB1_signature=$(echo 'CD24,CDH1')
selected_clusters=$(echo '3_1')
selected_comparisons=$(echo 'Genotype_A_G_vs_G_G,Genotype_A_A_vs_G_G,Genotype_Del16_vs_G_G,Genotype_Del80_vs_G_G')


# --dependency=afterany:$myjobid_collect_DA

myjobid_volcano_plots=$(sbatch --dependency=afterany:$myjobid_collect_DA --job-name=$name_volcano_plots --output=$outfile_volcano_plots --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_volcano_plots --platelet_volume_genes $platelet_volume_genes --DA_results $DA_results --type $type --out $output_dir --selected_clusters $selected_clusters --selected_comparisons $selected_comparisons --platelet_genes $platelet_genes --EZH2_signature $EZH2_signature --AKT_signature $AKT_signature --ZEB1_signature $ZEB1_signature")
myjobid_seff_volcano_plots=$(sbatch --dependency=afterany:$myjobid_volcano_plots --open-mode=append --output=$outfile_volcano_plots --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_volcano_plots >> $outfile_volcano_plots")





