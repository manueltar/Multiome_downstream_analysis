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

# rm -rf $Log_files
# mkdir -p $Log_files



conda activate multiome_downstream

### Peak_extraction

type=$(echo "Peak_extraction")
outfile_Peak_extraction=$(echo "$Log_files""outfile_0_""$type"".log")
touch $outfile_Peak_extraction
echo -n "" > $outfile_Peak_extraction
name_Peak_extraction=$(echo "$type""_job")
seff_name=$(echo "seff_""$type")


Rscript_Peak_extraction=$(echo "$Rscripts_path""205_Extract_Peaks_for_DA_DESeq2_Seurat_time_covariate_v2.R")


Seurat_object=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/merged_downstream_analysis_FOR_DANIELA.rds")
ery_genes_array=$(echo 'HBQ1,HBZ,HBA1,HBA2,HBG1,HBG2,MYL4,ANKRD26')
myeloid_genes_array=$(echo 'CCR7,IL33,SLC44A1,AXL,INHBA,IGF1R,KIF21A,IL33,IFI16,ZEB1,ALDH1A2,IL6ST,TJP2,LTBP1,BACH2,CCL3,EPB41')
megak_genes_array=$(echo 'GP9,KIF15,KIF22,GP6,CMTM5,KLF1,TUBB1,ITGA2B,FGF13,SPTA1,GFI1B,MYH9,MYLK,ITGB3,KCNT2,ADAM10,MEF2C,LRP12,ATP2C1,SPTB,FYB1,TBL1X,BCOR,TBXA2R,GP6,FCGR2A,RHAG,GP1BB,EGF,GNAQ')
marker_genes_array=$(echo 'ANAPC7,STIL,KIF15,HBQ1,HBZ,HBA2,GP6,FYB1,CCR7,IL1B,IL33,GYPA,ITGA2B')
ANAPC_genes=$(echo 'ANAPC1,ANAPC10,ANAPC11,ANAPC13,ANAPC15,ANAPC16,ANAPC2,ANAPC4,ANAPC5,ANAPC7')
EZH2_signature=$(echo 'EZH2,XRCC2,CDKN1A,CDKN2D,ZEB1,ZEB2,TWIST')
CUX1=$(echo "CUX1")



Selected_genes_classified=$(echo '/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Per_cluster_DESeq2_time_covariate_model/Selected_genes_classified.rds')



myjobid_Peak_extraction=$(sbatch --job-name=$name_Peak_extraction --output=$outfile_Peak_extraction --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=6 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_Peak_extraction --Seurat_object $Seurat_object --ery_genes_array $ery_genes_array --myeloid_genes_array $myeloid_genes_array --megak_genes_array $megak_genes_array --marker_genes_array $marker_genes_array --ANAPC_genes $ANAPC_genes --EZH2_signature $EZH2_signature --CUX1 $CUX1 --Selected_genes_classified $Selected_genes_classified --type $type --out $output_dir")
myjobid_seff_Peak_extraction=$(sbatch --dependency=afterany:$myjobid_Peak_extraction --open-mode=append --output=$outfile_Peak_extraction --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Peak_extraction >> $outfile_Peak_extraction")

conda deactivate


#### Assign_K562_promoters_to_Ensembl_gene_TSS #################

module load R/4.1.0

Rscript_Assign_K562_promoters_to_Ensembl_gene_TSS=$(echo "$Rscripts_path""177_Classify_DA_peaks_with_promoters_and_other_features_v2_function_1_Assign_K562_promoters_to_Ensembl_gene_TSS.R")

type=$(echo "Assign_K562_promoters_to_Ensembl_gene_TSS")


outfile_Assign_K562_promoters_to_Ensembl_gene_TSS=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_Assign_K562_promoters_to_Ensembl_gene_TSS
echo -n "" > $outfile_Assign_K562_promoters_to_Ensembl_gene_TSS
name_Assign_K562_promoters_to_Ensembl_gene_TSS=$(echo "$type""_job")
seff_name=$(echo "seff_""$type")


ensembl_gtf=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/reference_files/Homo_sapiens.GRCh38.111.gtf")
Promoter_distance_to_TSS=$(echo "2500")
tracking_genes=$(echo 'CUX1,GFI1B,RUNX1')
K562_Regulatory_Build=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/reference_files/homo_sapiens.GRCh38.K562.Regulatory_Build.regulatory_activity.20221007.gff")


myjobid_Assign_K562_promoters_to_Ensembl_gene_TSS=$(sbatch --job-name=$name_Assign_K562_promoters_to_Ensembl_gene_TSS --output=$outfile_Assign_K562_promoters_to_Ensembl_gene_TSS --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Assign_K562_promoters_to_Ensembl_gene_TSS --Promoter_distance_to_TSS $Promoter_distance_to_TSS --K562_Regulatory_Build $K562_Regulatory_Build --tracking_genes $tracking_genes --ensembl_gtf $ensembl_gtf --type $type --out $output_dir")
myjobid_seff_Assign_K562_promoters_to_Ensembl_gene_TSS=$(sbatch --dependency=afterany:$myjobid_Assign_K562_promoters_to_Ensembl_gene_TSS --open-mode=append --output=$outfile_Assign_K562_promoters_to_Ensembl_gene_TSS --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Assign_K562_promoters_to_Ensembl_gene_TSS >> $outfile_Assign_K562_promoters_to_Ensembl_gene_TSS")


#### Assign_Assay_peaks_to_Ensembl_gene_TSS #################

module load R/4.1.0
 
Rscript_Assign_Assay_peaks_to_Ensembl_gene_TSS=$(echo "$Rscripts_path""177_Classify_DA_peaks_with_promoters_and_other_features_v2_function_2_Assign_Assay_peaks_to_Ensembl_gene_TSS.R")

type=$(echo "Assign_Assay_peaks_to_Ensembl_gene_TSS")


outfile_Assign_Assay_peaks_to_Ensembl_gene_TSS=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_Assign_Assay_peaks_to_Ensembl_gene_TSS
echo -n "" > $outfile_Assign_Assay_peaks_to_Ensembl_gene_TSS
name_Assign_Assay_peaks_to_Ensembl_gene_TSS=$(echo "$type""_job")
seff_name=$(echo "seff_""$type")


ensembl_gtf=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/reference_files/Homo_sapiens.GRCh38.111.gtf")
distance_to_TSS=$(echo "2500")
ALL_PoI=$(echo "$output_dir""ALL_PoI.rds")
tracking_genes=$(echo 'CUX1,GFI1B,RUNX1')

# --dependency=afterany:$myjobid_Assign_K562_promoters_to_Ensembl_gene_TSS:$myjobid_Peak_extraction

myjobid_Assign_Assay_peaks_to_Ensembl_gene_TSS=$(sbatch --job-name $name_Assign_Assay_peaks_to_Ensembl_gene_TSS --dependency=afterany:$myjobid_Assign_K562_promoters_to_Ensembl_gene_TSS:$myjobid_Peak_extraction --output=$outfile_Assign_Assay_peaks_to_Ensembl_gene_TSS --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Assign_Assay_peaks_to_Ensembl_gene_TSS --distance_to_TSS $distance_to_TSS --tracking_genes $tracking_genes --ensembl_gtf $ensembl_gtf --ALL_PoI $ALL_PoI --type $type --out $output_dir")
myjobid_seff_Assign_Assay_peaks_to_Ensembl_gene_TSS=$(sbatch --dependency=afterany:$myjobid_Assign_Assay_peaks_to_Ensembl_gene_TSS --open-mode=append --output=$outfile_Assign_Assay_peaks_to_Ensembl_gene_TSS --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Assign_Assay_peaks_to_Ensembl_gene_TSS >> $outfile_Assign_Assay_peaks_to_Ensembl_gene_TSS")


#### classify_TSS_and_linked_peaks_to_K562_regulatory_features #################

module load R/4.1.0

Rscript_classify_TSS_and_linked_peaks_to_K562_regulatory_features=$(echo "$Rscripts_path""177_Classify_DA_peaks_with_promoters_and_other_features_v2_function_3_1_Assign_TSS_and_linked_Assay_peaks_to_K562_regulatory_features.R")

type=$(echo "classify_TSS_and_linked_peaks_to_K562_regulatory_features")


outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features=$(echo "$Log_files""outfile_3_""$type"".log")
touch $outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features
echo -n "" > $outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features
name_classify_TSS_and_linked_peaks_to_K562_regulatory_features=$(echo "$type""_job")
seff_name=$(echo "seff_""$type")


TSS_PoI=$(echo "$output_dir""ALL_PoI_in_TSS.rds")
Linked_PoI=$(echo "$output_dir""Linked_peak_to_selected_genes.rds")
tracking_genes=$(echo 'CUX1,GFI1B,RUNX1')
K562_Regulatory_Build=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/reference_files/homo_sapiens.GRCh38.K562.Regulatory_Build.regulatory_activity.20221007.gff")
K562_Promoters_in_TSS=$(echo "$output_dir""reference_files/Ensembl_promoters_K562_linked_to_genes.rds")

# --dependency=afterany:$myjobid_Assign_Assay_peaks_to_Ensembl_gene_TSS:$myjobid_Assign_K562_promoters_to_Ensembl_gene_TSS

myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features=$(sbatch --dependency=afterany:$myjobid_Assign_Assay_peaks_to_Ensembl_gene_TSS:$myjobid_Assign_K562_promoters_to_Ensembl_gene_TSS  --job-name=$name_classify_TSS_and_linked_peaks_to_K562_regulatory_features --output=$outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_classify_TSS_and_linked_peaks_to_K562_regulatory_features --K562_Regulatory_Build $K562_Regulatory_Build --tracking_genes $tracking_genes --TSS_PoI $TSS_PoI --Linked_PoI $Linked_PoI --K562_Promoters_in_TSS $K562_Promoters_in_TSS --type $type --out $output_dir")
myjobid_seff_classify_TSS_and_linked_peaks_to_K562_regulatory_features=$(sbatch --dependency=afterany:$myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features --open-mode=append --output=$outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features >> $outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features")

#### classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS #################

module load R/4.1.0

Rscript_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS=$(echo "$Rscripts_path""177_Classify_DA_peaks_with_promoters_and_other_features_v2_function_3_2_Assign_TSS_and_linked_Assay_peaks_to_K562_regulatory_features_NON_PROMOTERS.R")

type=$(echo "classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS")


outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS=$(echo "$Log_files""outfile_4_""$type"".log")
touch $outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS
echo -n "" > $outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS
name_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS=$(echo "$type""_job")
seff_name=$(echo "seff_""$type")



TSS_PoI=$(echo "$output_dir""ALL_PoI_in_TSS.rds")
Linked_PoI=$(echo "$output_dir""Linked_peak_to_selected_genes.rds")
tracking_genes=$(echo 'CUX1,GFI1B,RUNX1')
K562_Regulatory_Build=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/reference_files/homo_sapiens.GRCh38.K562.Regulatory_Build.regulatory_activity.20221007.gff")
PoI_promoters=$(echo "$output_dir""PoI_concordant_promoters.rds")
ensembl_gtf=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/reference_files/Homo_sapiens.GRCh38.111.gtf")

# --dependency=afterany:$myjobid_Assign_Assay_peaks_to_Ensembl_gene_TSS:$myjobid_Assign_K562_promoters_to_Ensembl_gene_TSS

myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS=$(sbatch --dependency=afterany:$myjobid_Assign_Assay_peaks_to_Ensembl_gene_TSS:$myjobid_Assign_K562_promoters_to_Ensembl_gene_TSS --job-name=$name_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS --output=$outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS --K562_Regulatory_Build $K562_Regulatory_Build --tracking_genes $tracking_genes --TSS_PoI $TSS_PoI --Linked_PoI $Linked_PoI --PoI_promoters $PoI_promoters --ensembl_gtf $ensembl_gtf --type $type --out $output_dir")
myjobid_seff_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS=$(sbatch --dependency=afterany:$myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS --open-mode=append --output=$outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS >> $outfile_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS")


#### Classify_THE_REST #################

module load R/4.1.0

Rscript_Classify_THE_REST=$(echo "$Rscripts_path""177_Classify_DA_peaks_with_promoters_and_other_features_v2_function_3_3_Assign_Rest_of_ALL_PoI.R")

type=$(echo "Classify_THE_REST")


outfile_Classify_THE_REST=$(echo "$Log_files""outfile_5_""$type"".log")
touch $outfile_Classify_THE_REST
echo -n "" > $outfile_Classify_THE_REST
name_Classify_THE_REST=$(echo "$type""_job")
seff_name=$(echo "seff_""$type")

ALL_PoI=$(echo "$output_dir""ALL_PoI.rds")
Non_promoters_PoI=$(echo "$output_dir""PoI_non_promoters.rds")
PoI_promoters=$(echo "$output_dir""PoI_concordant_promoters.rds")
tracking_genes=$(echo 'CUX1,GFI1B,RUNX1')
K562_Regulatory_Build=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/reference_files/homo_sapiens.GRCh38.K562.Regulatory_Build.regulatory_activity.20221007.gff")


# --dependency=afterany:$myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features:$myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS

myjobid_Classify_THE_REST=$(sbatch --dependency=afterany:$myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features:$myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS --job-name=$name_Classify_THE_REST --dependency=afterany:$myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features:$myjobid_classify_TSS_and_linked_peaks_to_K562_regulatory_features_NON_PROMOTERS  --output=$outfile_Classify_THE_REST --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Classify_THE_REST --K562_Regulatory_Build $K562_Regulatory_Build --tracking_genes $tracking_genes --ALL_PoI $ALL_PoI --Non_promoters_PoI $Non_promoters_PoI --PoI_promoters $PoI_promoters --type $type --out $output_dir")
myjobid_seff_Classify_THE_REST=$(sbatch --dependency=afterany:$myjobid_Classify_THE_REST --open-mode=append --output=$outfile_Classify_THE_REST --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Classify_THE_REST >> $outfile_Classify_THE_REST")


#### Intersect_SNP #################

module load R/4.1.0

Rscript_Intersect_SNP=$(echo "$Rscripts_path""177_Classify_DA_peaks_with_promoters_and_other_features_v2_function_4_Intersect_SNPs.R")

type=$(echo "Intersect_SNP")


outfile_Intersect_SNP=$(echo "$Log_files""outfile_6_""$type"".log")
touch $outfile_Intersect_SNP
echo -n "" > $outfile_Intersect_SNP
name_Intersect_SNP=$(echo "$type""_job")
seff_name=$(echo "seff_""$type")

Master_peak_file=$(echo "$output_dir""Master_peak_file.rds")
selected_variants=$(echo 'rs139141690__chr7_101856650_G_A')

# --dependency=afterany:$myjobid_Classify_THE_REST

myjobid_Intersect_SNP=$(sbatch --dependency=afterany:$myjobid_Classify_THE_REST --job-name=$name_Intersect_SNP --output=$outfile_Intersect_SNP --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Intersect_SNP --selected_variants $selected_variants --Master_peak_file $Master_peak_file --type $type --out $output_dir")
myjobid_seff_Intersect_SNP=$(sbatch --dependency=afterany:$myjobid_Intersect_SNP --open-mode=append --output=$outfile_Intersect_SNP --job-name=seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Intersect_SNP >> $outfile_Intersect_SNP")

#### Number_the_peaks #################

module load R/4.1.0

Rscript_Number_the_peaks=$(echo "$Rscripts_path""177_Classify_DA_peaks_with_promoters_and_other_features_v2_function_5_Number_the_peaks.R")

type=$(echo "Number_the_peaks")


outfile_Number_the_peaks=$(echo "$Log_files""outfile_7_""$type"".log")
touch $outfile_Number_the_peaks
echo -n "" > $outfile_Number_the_peaks
name_Number_the_peaks=$(echo "$type""_job")
seff_name=$(echo "seff_""$type")

Master_peak_file=$(echo "$output_dir""Master_peak_file_with_SNP.rds")

# --dependency=afterany:$myjobid_Intersect_SNP

myjobid_Number_the_peaks=$(sbatch --dependency=afterany:$myjobid_Intersect_SNP --job-name=$name_Number_the_peaks --output=$outfile_Number_the_peaks --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Number_the_peaks --Master_peak_file $Master_peak_file --type $type --out $output_dir")
myjobid_seff_Number_the_peaks=$(sbatch --dependency=afterany:$myjobid_Number_the_peaks --open-mode=append --output=$outfile_Number_the_peaks --job-name=seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Number_the_peaks >> $outfile_Number_the_peaks")



#### Build_Input_regions #################

module load R/4.1.0

Rscript_Build_Input_regions=$(echo "$Rscripts_path""206_Build_Input_regions_for_DA_DESeq2_time_covariate.R")

type=$(echo "Build_Input_regions")


outfile_Build_Input_regions=$(echo "$Log_files""outfile_8_""$type"".log")
touch $outfile_Build_Input_regions
echo -n "" > $outfile_Build_Input_regions
name_Build_Input_regions=$(echo "$type""_job")
seff_name=$(echo "seff_""$type")

Master_peak_file=$(echo "$output_dir""Master_peak_file_with_SNP_numbered_peaks.rds")
ery_genes_array=$(echo 'HBQ1,HBZ,HBA1,HBA2,HBG1,HBG2,MYL4,ANKRD26')
myeloid_genes_array=$(echo 'CCR7,IL33,SLC44A1,AXL,INHBA,IGF1R,KIF21A,IL33,IFI16,ZEB1,ALDH1A2,IL6ST,TJP2,LTBP1,BACH2,CCL3,EPB41')
megak_genes_array=$(echo 'GP9,KIF15,KIF22,GP6,CMTM5,KLF1,TUBB1,ITGA2B,FGF13,SPTA1,GFI1B,MYH9,MYLK,ITGB3,KCNT2,ADAM10,MEF2C,LRP12,ATP2C1,SPTB,FYB1,TBL1X,BCOR,TBXA2R,GP6,FCGR2A,RHAG,GP1BB,EGF,GNAQ')
marker_genes_array=$(echo 'GYPA,ITGA2B,GP1BB,MECOM,PRDX1,ACSL4')
ANAPC_genes=$(echo 'ANAPC1,ANAPC10,ANAPC11,ANAPC13,ANAPC15,ANAPC16,ANAPC2,ANAPC4,ANAPC5,ANAPC7')
ANAPC_genes=$(echo 'ANAPC1,ANAPC10,ANAPC11,ANAPC13,ANAPC15,ANAPC16,ANAPC2,ANAPC4,ANAPC5,ANAPC7,XRCC2,EZH2,CUX1')

Selected_genes_classified=$(echo '/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/Per_cluster_DESeq2_time_covariate_model/Selected_genes_classified.rds')

ensembl_gtf=$(echo "/group/soranzo/manuel.tardaguila/SC_RNA_seq/k562_multiome/NEW_object_output/reference_files/Homo_sapiens.GRCh38.111.gtf")

# --dependency=afterany:$myjobid_Number_the_peaks

myjobid_Build_Input_regions=$(sbatch --dependency=afterany:$myjobid_Number_the_peaks --job-name=$name_Build_Input_regions --output=$outfile_Build_Input_regions --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Build_Input_regions --Master_peak_file $Master_peak_file --ensembl_gtf $ensembl_gtf --ery_genes_array $ery_genes_array --myeloid_genes_array $myeloid_genes_array --megak_genes_array $megak_genes_array --marker_genes_array $marker_genes_array --ANAPC_genes $ANAPC_genes --Selected_genes_classified $Selected_genes_classified --type $type --out $output_dir")
myjobid_seff_Build_Input_regions=$(sbatch --dependency=afterany:$myjobid_Build_Input_regions --open-mode=append --output=$outfile_Build_Input_regions --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Build_Input_regions >> $outfile_Build_Input_regions")




