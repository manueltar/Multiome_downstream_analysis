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




#####################Build_IR_and_GF_files

module load R/4.1.0

Rscript_Build_IR_and_GF_files=$(echo "$Rscripts_path""212_Build_Input_regions_for_DA_DESeq2_time_covariate_ADD_PoIs.R")

type=$(echo "Build_IR_and_GF_files")


outfile_Build_IR_and_GF_files=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_Build_IR_and_GF_files
echo -n "" > $outfile_Build_IR_and_GF_files
name_Build_IR_and_GF_files=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Peaks_to_display=$(echo 'CUX1\>Linked_Peak__3__CUX1\|Linked_Peak__7__CUX1\|Linked_Peak__5__CUX1\|Linked_Peak__2__CUX1,GUSB\>TSS_Peak__1__GUSB,XRCC2\>TSS_Peak__1__XRCC2\|Linked_Peak__1__XRCC2,MOCOS\>TSS_Peak__1__COSMOC\;MOCOS\|Linked_Peak__2__MOCOS,PDE10A\>Linked_Peak__2__PDE10A,EZH2\>Linked_Peak__1__EZH2,SMAD2\>TSS_Peak__1__SMAD2')
Master_peak_file_with_SNP_numbered_peaks=$(echo '/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Peak_extraction_and_classification_for_time_covariate_model/Master_peak_file_with_SNP_numbered_peaks.rds')
Input_regions=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Peak_extraction_and_classification_for_time_covariate_model/Input_regions.rds")
Linked_peak_to_selected_genes=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Peak_extraction_and_classification_for_time_covariate_model/Linked_peak_to_selected_genes.rds")
K562_Regulatory_Build=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/reference_files/homo_sapiens.GRCh38.K562.Regulatory_Build.regulatory_activity.20221007.gff")
ensembl_gtf=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/reference_files/Homo_sapiens.GRCh38.111.gtf")


myjobid_Build_IR_and_GF_files=$(sbatch --job-name=$name_Build_IR_and_GF_files --output=$outfile_Build_IR_and_GF_files --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Build_IR_and_GF_files --Linked_peak_to_selected_genes $Linked_peak_to_selected_genes --ensembl_gtf $ensembl_gtf --K562_Regulatory_Build $K562_Regulatory_Build --Peaks_to_display $Peaks_to_display --Input_regions $Input_regions --Master_peak_file_with_SNP_numbered_peaks $Master_peak_file_with_SNP_numbered_peaks --type $type --out $output_dir")
myjobid_seff_Build_IR_and_GF_files=$(sbatch --dependency=afterany:$myjobid_Build_IR_and_GF_files --open-mode=append --output=$outfile_Build_IR_and_GF_files --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Build_IR_and_GF_files >> $outfile_Build_IR_and_GF_files")


#################### Print_situation_plots

conda activate multiome_downstream

type=$(echo "Situation_plot_printer")
outfile_Print_situation_plots=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_Print_situation_plots
echo -n "" > $outfile_Print_situation_plots
name_Print_situation_plots=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")



Rscript_Print_situation_plots=$(echo "$Rscripts_path""253_Situation_plot_printer.R")

Seurat_object=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/merged_downstream_analysis_FOR_DANIELA.rds")
Genomic_features=$(echo "$output_dir""ALL_regions_Genomic_features.rds")
Input_regions=$(echo "$output_dir""Input_Regions_DEF.rds")
Master_peak_file_with_SNP_numbered_peaks=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Peak_extraction_and_classification_for_time_covariate_model/Master_peak_file_with_SNP_numbered_peaks.rds")
metadata_GeneEXP=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Per_cluster_DESeq2_time_covariate_model/Metadata.rds")
GeneEXP=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Per_cluster_DESeq2_time_covariate_model/Normalised_count_matrix.tsv")
metadata_ATAC=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Per_cluster_DESeq2_time_covariate_model_ATAC/Metadata.rds")
ATAC=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Per_cluster_DESeq2_time_covariate_model_ATAC/Normalised_count_matrix.tsv")

#  --dependency=afterany:$myjobid_Build_IR_and_GF_files

myjobid_Print_situation_plots=$(sbatch --dependency=afterany:$myjobid_Build_IR_and_GF_files --job-name $name_Print_situation_plots --output=$outfile_Print_situation_plots --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=12 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_Print_situation_plots --Seurat_object $Seurat_object --Genomic_features $Genomic_features --Input_regions $Input_regions --Master_peak_file_with_SNP_numbered_peaks $Master_peak_file_with_SNP_numbered_peaks --metadata_GeneEXP $metadata_GeneEXP --GeneEXP $GeneEXP --metadata_ATAC $metadata_ATAC --ATAC $ATAC --type $type --out $output_dir")
myjobid_seff_Print_situation_plots=$(sbatch --dependency=afterany:$myjobid_Print_situation_plots --open-mode=append --output=$outfile_Print_situation_plots --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Print_situation_plots >> $outfile_Print_situation_plots")


conda deactivate

















