# Peripheral Blood Nuclear Cells (PBMCs) Cell-Type Classification

Single-cell RNA sequencing technology quantifies transcriptomic expressions at single-cell level, leading to discoveries rooted in cell-to-cell heterogeneity. Cell-type classification, or assignment of cell type labels to single cells, is an important step in analyzing the RNA sequencing data. For example, correctly identifying cells within the tumor microenvironment can provide oncologists with a snapshot of how a patient’s immune system reacts to the tumor.

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
    preprocess.R [-h] [-v] [--id] file_prefix [--data_dir] /path/to/data_dir [--figure_dir] /path/to/figure_dir
~~~

## Models ##
 
~~~
    build_model.R [-h] [-v] [--genes] genes [--size] size [--model] rf_default [--data_dir] /path/to/data_dir [--results_dir] /path/to/results_dir
~~~

## Evaluation & Visualization ##

~~~
    compute_metrics.R [-h] [-v] [--file] filename_of_model [--data_dir] /path/to/data_dir [--results_dir] /path/to/results_dir
    plot_tsne.R [-h] [-v] [--file] filename_of_model [--data_dir] /path/to/data_dir
~~~ 
