# load relevant libraries
# ----------------------------
suppressMessages(require(Matrix))
suppressMessages(require(ggplot2))
suppressMessages(require(Rtsne))
suppressMessages(require(svd))
suppressMessages(require(dplyr))
suppressMessages(require(plyr))
suppressMessages(require(data.table))
suppressMessages(require(pheatmap))
suppressMessages(require(ggpubr))
suppressMessages(require(argparse))
source('util.R')
source('select_pure_pbmc.R')

parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
    dest="verbose", help="Print little output")
parser$add_argument("--f", "--file", type="character",  
    help = "File name of model")
parser$add_argument("--d", "--data_dir", type="character",  
    help = "Directory that contains data")
args <- parser$parse_args()

# -------------------------------------
# parse variables
# -------------------------------------

DATA_DIR <- args$data_dir
FILE <- args$file

annot <- read.csv(file.path(DATA_DIR,'68k_pbmc_barcodes_annotation.tsv'),sep='\t')

# ------------------------------------------------------------
# load 68k PBMC data, 11 purified PBMC data and meta-data
# ------------------------------------------------------------
pbmc_68k <- readRDS(file.path(DATA_DIR,'pbmc68k_data.rds'))
pure_11 <- readRDS(file.path(DATA_DIR,'all_pure_select_11types.rds'))
all_data <- pbmc_68k$all_data
purified_ref_11 <- load_purified_pbmc_types(pure_11,pbmc_68k$ens_genes)
# --------------------------------------------------------------------------------------
# normalize by RNA content (umi counts) and select the top 1000 most variable genes
# --------------------------------------------------------------------------------------
m<-all_data[[1]]$hg19$mat
l<-.normalize_by_umi(m)   
m_n<-l$m
df<-.get_variable_gene(m_n) 
disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[1000]
df$used<-df$dispersion_norm >= disp_cut_off
# --------------------------------------------
# use top 1000 variable genes for PCA 
# --------------------------------------------
set.seed(0)
m_n_1000<-m_n[,head(order(-df$dispersion_norm),1000)]
pca_n_1000<-.do_propack(m_n_1000,50)
# --------------------------------------------
# generate 2-D tSNE embedding
# this step may take a long time
# --------------------------------------------
tsne_n_1000<-Rtsne(pca_n_1000$pca,pca=F)
tdf_n_1000<-data.frame(tsne_n_1000$Y)
# --------------------------------------------------
# use k-means clustering to specify populations
# --------------------------------------------------
set.seed(0)
k_n_1000<-kmeans(pca_n_1000$pca,10,iter.max=150,algorithm="MacQueen")
tdf_n_1000$k<-k_n_1000$cluster
plt_kmeans <- ggplot(tdf_n_1000,aes(X1,X2,col=as.factor(k)))+geom_point(size=0,alpha=0.6)+theme_classic()+
   scale_color_manual(values=c("#FB9A99","#FF7F00","yellow","orchid","grey",
                               "red","dodgerblue2","tan4","green4","#99c9fb"))
# ---------------------------------------------------------------------------------------------------------------------------
# assign IDs by comparing the transcriptome profile of each cell to the reference profile from purified PBMC populations
# ---------------------------------------------------------------------------------------------------------------------------
m_filt<-m_n_1000
use_genes_n<-order(-df$dispersion_norm)
use_genes_n_id<-all_data[[1]]$hg19$gene_symbols[l$use_genes][order(-df$dispersion_norm)]
use_genes_n_ens<-all_data[[1]]$hg19$genes[l$use_genes][order(-df$dispersion_norm)]
z_1000_11<-.compare_by_cor(m_filt,use_genes_n_ens[1:1000],purified_ref_11) 
# reassign IDs, as there're some overlaps in the purified pbmc populations
test<-.reassign_pbmc_11(z_1000_11)
cls_id<-factor(colnames(z_1000_11)[test])
tdf_n_1000$cls_id<-cls_id
# adjust ordering of cells for plotting aesthetics
tdf_mod <- tdf_n_1000[tdf_n_1000$cls_id!='CD4+/CD45RA+/CD25- Naive T',]
tdf_mod <- rbind(tdf_mod,tdf_n_1000[tdf_n_1000$cls_id=='CD4+/CD45RA+/CD25- Naive T',])
tdf_mod_2 <- tdf_mod[tdf_mod$cls_id!='CD56+ NK',]
tdf_mod_2 <- rbind(tdf_mod_2,tdf_mod[tdf_mod$cls_id=='CD56+ NK',])
plt_marker <- ggplot(tdf_mod_2,aes(X1,X2,col=cls_id))+geom_point(size=0,alpha=1)+theme_classic()+.set_pbmc_color_11()

# ---------------------------------------------------------------------------------------------------------------------------
# assign IDs by using output from machine learning model
# ---------------------------------------------------------------------------------------------------------------------------
res <- readRDS(file.path(DATA_DIR,FILE))
tdf_n_1000$cls_id <- res[[3]]
tdf_mod <- tdf_n_1000[tdf_n_1000$cls_id!='CD4+/CD45RA+/CD25- Naive T',]
tdf_mod <- rbind(tdf_mod,tdf_n_1000[tdf_n_1000$cls_id=='CD4+/CD45RA+/CD25- Naive T',])
tdf_mod_2 <- tdf_mod[tdf_mod$cls_id!='CD56+ NK',]
tdf_mod_2 <- rbind(tdf_mod_2,tdf_mod[tdf_mod$cls_id=='CD56+ NK',])
plt_model <- ggplot(tdf_mod_2,aes(X1,X2,col=cls_id))+geom_point(size=0,alpha=1)+theme_classic()+.set_pbmc_color_11()

ggarrange(plt_kmeans,plt_marker,plt_marker)
