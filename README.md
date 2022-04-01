# Peripheral Blood Nuclear Cells (PBMCs) cell-type classification using scRNA-sequencing data 

Single-cell RNA sequencing (scRNA-seq) allows for the profiling of gene expression at single-cell resolution. It has the potential to uncover rare cell populations, reveal regulatory relationships between genes, and track the trajectories of distinct cell lineages across development. Cell-type classification is a critical step in analyzing the scRNA-sequencing data. For instance, understanding the cells that make up the tumor microenvironment can inform clinicians how a patient’s immune system may react to various therapeutic strategies. In this study our goal is to develop a machine learning-based classifier for single-cell cell-type classification of peripheral blood mononuclear cells (PBMCs) and conduct various ablation experiments. 

## Data ##
In order to run these analyses you will need to download the following data

Training Set: scRNA-seq data of bead-enriched sub-populations of PBMCs from Donor A 

Test Set: scRNA-seq data of 68k PBMCs from Donor A

All of the data were processed by the Cell Ranger 1.1 pipeline. The processed data used for the single cell RNA-seq secondary analysis are stored in three R data files:

* [pbmc68k_data.rds](https://cf.10xgenomics.com/samples/cell/pbmc68k_rds/pbmc68k_data.rds) (77MB): consists of the gene expression profiles of ~68k PBMCs
* [all_pure_pbmc_data.rds](https://cf.10xgenomics.com/samples/cell/pbmc68k_rds/all_pure_pbmc_data.rds) (1.5GB): consists of the gene expression profiles of 10 bead-enriched PBMC samples
* [all_pure_select_11types.rds](https://cf.10xgenomics.com/samples/cell/pbmc68k_rds/all_pure_select_11types.rds) (687KB): consists of the gene expression and meta-information of the 11 sub-populations of PBMC identified from the 10 samples
* [68k_pbmc_barcodes_annotation.tsv](https://raw.githubusercontent.com/10XGenomics/single-cell-3prime-paper/master/pbmc68k_analysis/68k_pbmc_barcodes_annotation.tsv) (4.58MB): consists of the tsv file with TSNE coordinates, cell barcodes and annotations for the ~68k PBMCs

## Pre-processing ##

~~~
    Rscript preprocess.R [-h] [-v] \
        [--id] file_prefix \
        [--data_dir] /path/to/data_dir \
        [--figure_dir] /path/to/figure_dir
~~~

## Models ##
 
~~~
    Rscript build_model.R [-h] [-v] \
        [--genes] e.g. "100", "diff1k", "all" \
        [--size] e.g.  "1k", "2500", "all"\
        [--model] e.g. "rf_default", "rf_tuned", "glm", "svmLinear", "svmRadial" \
        [--data_dir] /path/to/data_dir \
        [--results_dir] /path/to/results_dir
~~~

## Evaluation & Visualization ##

~~~
    Rscript compute_metrics.R [-h] [-v] \
        [--file] filename_of_model \
        [--data_dir] /path/to/data_dir \
        [--results_dir] /path/to/results_dir
        
    Rscript plot_tsne.R [-h] [-v] \
        [--file] filename_of_model \
        [--data_dir] /path/to/data_dir
~~~ 
