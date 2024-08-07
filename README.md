
## See multiome_downstream_dependencies.txt for the list of the dependencies of the conda environment needed to run Seurat and Signac


## Perform initial UMAP plot, Stacked barplot and Dotplot

\$ bash ~/Scripts/Wraper_scripts/75_Minus_16bp_UMAP_and_stacked_barplot.sh \<path_to_analysis\> \<analysis_name\>


## Perform the Differential Expression (DE) analysis comparing per cluster the gene expression of every genotype against the wt genotype

\$ bash ~/Scripts/Wraper_scripts/76_Per_cluster_DESeq2_Seurat_analysis_time_covariate_vs_Genotype.sh \<path_to_analysis\> \<analysis_name\>

## Extract linked peaks to DE genes, classify them into TSS and non TSS, overlap against the ENSEMBL K562_Regulatory_Build (https://ftp.ensembl.org/pub/release-111/regulation/homo_sapiens/RegulatoryFeatureActivity/K562/), intersect with SNPs and number them. Also it build the initial input table for locus representation

\$ bash ~/Scripts/Wraper_scripts/66_Peaks_extraction_and_classification_from_DESeq2_Seurat_analysis_time_covariate_vs_Genotype.sh \<path_to_analysis\> \<analysis_name\>

## Perform Differential Accessibility (DA) analysis comparing per cluster the peak accessibility of every genotype against the wt genotype

\$ bash ~/Scripts/Wraper_scripts/77_Per_cluster_DESeq2_Seurat_analysis_time_covariate_vs_Genotype_ATAC.sh \<path_to_analysis\> \<analysis_name\>

## Print the situation plots with all the annotations for every locus plus the CPM values for the gene EXP and the chromatin accessibility

\$ bash ~/Scripts/Wraper_scripts/68_Per_locus_representation_DESeq2_Seurat_time_as_a_covariate.sh \<path_to_analysis\> \<analysis_name\>

## Unique DE genes

\$ awk -F"\t" '{if($8 != "NA") print $0}' DE_results.tsv|awk -F"\t" '{if($8 >= 1.3) print $0}'|cut -f1|sort|uniq -c|sort -rg|wc -l
1862

## Unique DA peaks

\$ awk -F"\t" '{if($8 != "NA") print $0}' DA_results.tsv|awk -F"\t" '{if($8 >= 1.3) print $0}'|cut -f1|sort|uniq -c|sort -rg|wc -l
485

## List of gene sets assayed is in List_of_tested_gene_sets.txt (n=803)

## Classification of peaks

\$ cut -f2 Master_peak_file_with_SNP_numbered_peaks.tsv|sort|uniq -c|sort -rg
   6495 Linked_Peak
   4284 Non_gene_associated_peak
   4072 TSS_Peak


