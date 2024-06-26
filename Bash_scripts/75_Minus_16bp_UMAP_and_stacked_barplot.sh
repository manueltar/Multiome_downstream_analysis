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

### multiome_minus_16bp


type=$(echo "multiome_minus_16bp")
outfile_multiome_minus_16bp=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_multiome_minus_16bp
echo -n "" > $outfile_multiome_minus_16bp
name_multiome_minus_16bp=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")


Rscript_multiome_minus_16bp=$(echo "$Rscripts_path""239_Minus_16bp_UMAP_and_stacked_barplot_v2.R")

Seurat_object=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/merged_downstream_analysis_FOR_DANIELA.rds")
#marker_genes=$(echo 'ITGA2B,GYPA,SLC1A4,GP1BB,MEF2C,ALDH1A2,ZFHX2,CREB3L1,PIP5K1C,MECOM,PRDX1,ACBD5,AARS2,GALNT2,NADK2')
marker_genes=$(echo 'ITGA2B,GYPA,IL33,IL1B,CCR7,FYB1,GP6,HBA2,HBZ,HBQ1,KIF15,STIL,ANAPC7')

myjobid_multiome_minus_16bp=$(sbatch --output=$outfile_multiome_minus_16bp --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=6 --mem-per-cpu=4096 --parsable --job-name $name_multiome_minus_16bp --wrap="Rscript $Rscript_multiome_minus_16bp --Seurat_object $Seurat_object --marker_genes $marker_genes --type $type --out $output_dir")
myjobid_seff_multiome_minus_16bp=$(sbatch --dependency=afterany:$myjobid_multiome_minus_16bp --open-mode=append --output=$outfile_multiome_minus_16bp --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_multiome_minus_16bp >> $outfile_multiome_minus_16bp")



conda deactivate



