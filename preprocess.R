rm(list=ls()) # clear workspace
# ----------------------------
# load relevant libraries
# ----------------------------
suppressMessages(require(Matrix))
suppressMessages(require(ggplot2))
suppressMessages(require(Rtsne))
suppressMessages(require(svd))
suppressMessages(require(plyr))
suppressMessages(require(dplyr))
suppressMessages(require(data.table))
suppressMessages(require(pheatmap))
suppressMessages(require(argparse))
source('util.R')
source('select_pure_pbmc.R')

parser <- ArgumentParser()

parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
    dest="verbose", help="Print little output")
parser$add_argument("--i", "--id", type="character",  
    help = "Prefix for filename")
parser$add_argument("--d", "--data_dir", type="character",  
    help = "Directory that contains data")
parser$add_argument("--r", "--figure_dir", type="character",  
    help = "Output directory for figures")
args <- parser$parse_args()

# -------------------------------------
# parse variables
# -------------------------------------
DATA_DIR <- args$data_dir
FIG_DIR <- args$figure_dir
id <- args$id

# -----------------------
# load purified data
# -----------------------
pure_pbmcs <- readRDS(file.path(DATA_DIR,'all_pure_pbmc_data.rds'))
all_data <- pure_pbmcs$all_data
all_json <- pure_pbmcs$all_json
all_metrics <- pure_pbmcs$all_metrics
all_mol_info <- pure_pbmcs$all_mol_info
genes <- all_data[[1]]$hg19$gene_symbols

# -------------------------------------------------------------------------
# downsample mapped reads/cell
# so that all samples have the same # of confidently mapped reads/cell
# -------------------------------------------------------------------------
set.seed(1)
rpc <- all_metrics %>% 
  mutate(conf_mapped_rpc=raw_rpc*conf_mapped_frac*good_bc_frac*good_umi_frac) %>%
  select(sample_id, description, conf_mapped_rpc)
tgt_rpc <- floor(min(rpc$conf_mapped_rpc)) # 13995

subsampled_purified_mats <- lapply(1:length(all_data), function(i) { # subsample the matrix to match tgt_rpc
    cat(sprintf("%d...\n", i))
    .downsample_gene_bc_mtx(all_json[[i]], all_data[[i]], all_mol_info[[i]], tgt_rpc, 'conf_mapped_reads')[[1]]
} )

# -------------------------------------------
# run PCA and tSNE on down-sampled data
# the steps take a long time to complete
# -------------------------------------------
# using 10 PCs provides enough information to filter out undesired populations from purified populations
set.seed(1)
all_pure_pca<-lapply(1:length(subsampled_purified_mats),function(i) pure_pca_i<-.do_propack(subsampled_purified_mats[[i]],100))
all_pure_tsne<-lapply(1:length(all_pure_pca),function(i) pure_tsne_i<-Rtsne(all_pure_pca[[i]]$pca,pca=F))
options(bitmapType='cairo')

# ------------------------------------
# curate the purified populations
# ------------------------------------
pure_id<-c("CD34+","CD56+ NK","CD4+/CD45RA+/CD25- Naive T", "CD4+/CD25 T Reg","CD8+/CD45RA+ Naive Cytotoxic",
           "CD4+/CD45RO+ Memory","CD8+ Cytotoxic T","CD19+ B","CD4+ T Helper2","CD14+ Monocyte","Dendritic")

sub_idx <-list(data.frame(sample=1, use=(get_pure_pop_idx(genes,pure_id[1],all_pure_pca[[1]],all_pure_tsne[[1]],FIG_DIR))), 
                 data.frame(sample=2, use=(get_pure_pop_idx(genes,pure_id[2],all_pure_pca[[2]],all_pure_tsne[[2]],FIG_DIR))),
                 data.frame(sample=3, use=(get_pure_pop_idx(genes,pure_id[3],all_pure_pca[[3]],all_pure_tsne[[3]],FIG_DIR))),
                 data.frame(sample=4, use=(get_pure_pop_idx(genes,pure_id[4],all_pure_pca[[4]],all_pure_tsne[[4]],FIG_DIR))),
                 data.frame(sample=5, use=(get_pure_pop_idx(genes,pure_id[5],all_pure_pca[[5]],all_pure_tsne[[5]],FIG_DIR))),
                 data.frame(sample=6, use=(get_pure_pop_idx(genes,pure_id[6],all_pure_pca[[6]],all_pure_tsne[[6]],FIG_DIR))),
                 data.frame(sample=7, use=(get_pure_pop_idx(genes,pure_id[7],all_pure_pca[[7]],all_pure_tsne[[7]],FIG_DIR))),
                 data.frame(sample=8, use=(get_pure_pop_idx(genes,pure_id[8],all_pure_pca[[8]],all_pure_tsne[[8]],FIG_DIR))),
                 data.frame(sample=9, use=(get_pure_pop_idx(genes,pure_id[9],all_pure_pca[[9]],all_pure_tsne[[9]],FIG_DIR))),
                 data.frame(sample=10,use=(get_pure_pop_idx(genes,pure_id[10],all_pure_pca[[10]],all_pure_tsne[[10]],FIG_DIR))),
                 data.frame(sample=10,use=(get_pure_pop_idx(genes,pure_id[11],all_pure_pca[[11]],all_pure_tsne[[11]],FIG_DIR))))
pure_select_11<-lapply(1:length(sub_idx),function(i) {subsampled_purified_mats[[sub_idx[[i]]$sample[1]]][sub_idx[[i]]$use,]})

# -------------------------------------------------
# Get training data all cells
# -------------------------------------------------
genes = "all" ; size = "all"
pure_use_genes<-which(colSums(do.call(rbind,lapply(pure_select_11,function(x) x)))>1)
pure_select_use_genes<-lapply(1:length(pure_select_11),function(i) pure_select_11[[i]][,pure_use_genes])
x_train<- do.call("rbind", pure_select_use_genes)
y_train <- unlist(lapply(1:length(pure_select_use_genes),function(i) rep(pure_id[[i]],dim(pure_select_use_genes[[i]])[1])))
saveRDS(list(x_train,y_train),file.path(DATA_DIR,paste(id,genes,size'train.rds',collapse="_")))

# -------------------------------------------------
# Get training data 1k cells
# -------------------------------------------------
genes = "all" ; size = "1k"
pure_select_use_genes_1k<-lapply(1:(length(pure_select_11)-2),function(i) pure_select_11[[i]][1:1000,pure_use_genes])
pure_select_use_genes_1k[[10]] <- pure_select_11[[10]][1:dim(pure_select_use_genes[[10]])[1],pure_use_genes]
pure_select_use_genes_1k[[11]] <- pure_select_11[[11]][1:dim(pure_select_use_genes[[11]])[1],pure_use_genes]
x_train_1k<- do.call("rbind", pure_select_use_genes_1k)
y_train_1k <- unlist(lapply(1:length(pure_select_use_genes_1k),function(i) rep(pure_id[[i]],dim(pure_select_use_genes_1k[[i]])[1])))
saveRDS(list(x_train,y_train),file.path(DATA_DIR,paste(id,genes,size,'train.rds',collapse="_")))

# ------------------------------------------------------------
# load 68k PBMC data, 11 purified PBMC data and meta-data
# ------------------------------------------------------------
pbmc_68k <- readRDS(file.path(DATA_DIR,'pbmc68k_data.rds'))
pure_11 <- readRDS(file.path(DATA_DIR,'all_pure_select_11types.rds'))
all_data <- pbmc_68k$all_data
annot <- read.csv('/hpf/largeprojects/davidm/vsubasri/single-cell-3prime-paper/pbmc68k_analysis/68k_pbmc_barcodes_annotation.tsv',sep='\t')

genes = "all" ; size = "all"
m<-all_data[[1]]$hg19$mat
x_test <- m[,pure_use_genes]
y_test <- annot$celltype
saveRDS(list(x_test,y_test),file.path(DATA_DIR,paste(id,genes,'test.rds',collapse="_")))

# --------------------------------------------------------------------------------------
# normalize testing data by RNA content (umi counts)
# --------------------------------------------------------------------------------------
genes = "all" ; size = "all"
test_l<-.normalize_by_umi(x_test)
test_m_n <- test_l$m
saveRDS(list(test_m_n,y_test),file.path(DATA_DIR,paste(id,'normalized',genes,'test.rds',collapse="_")))

# --------------------------------------------------------------------------------------
# normalize training data by RNA content (umi counts)
# --------------------------------------------------------------------------------------
train_l<-.normalize_by_umi(x_train)
train_m_n <- train_l$m
train_m_n <- train_m_n[,test_l$use_genes]
saveRDS(list(train_m_n,y_train),file.path(DATA_DIR,paste(id,'normalized', genes,size,'train.rds',collapse="_")))

# --------------------------------------------------------------------------------------
# select the top 1000 most variable genes
# --------------------------------------------------------------------------------------
df<-.get_variable_gene(train_m_n)
disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[1000]
df$used<-df$dispersion_norm >= disp_cut_off
train_m_n_1000<-train_m_n[,head(order(-df$dispersion_norm),1000)]
saveRDS(list(train_m_n_1000,y_train),file.path(DATA_DIR,paste(id,'normalized',genes,size,'train.rds',collapse="_")))

test_m_n_1000<-test_m_n[,head(order(-df$dispersion_norm),1000)]
saveRDS(list(test_m_n_1000,y_test),file.path(DATA_DIR,paste(id,'normalized',genes,'test.rds',collapse="_")))

#x_train_1k <- do.call("rbind",list(train_m_n_1000[1:1000,], train_m_n_1000[3793:4792,], train_m_n_1000[12178:13177,], train_m_n_1000[22657:23656,], train_m_n_1000[32920:33919,], train_m_n_1000[44873:45872,], train_m_n_1000[55097:56096,], train_m_n_1000[65306:66305,], train_m_n_1000[75391:76390,], train_m_n_1000[86604:87181,], train_m_n_1000[87182:87277,]))
#y_train_1k <- c(unlist(lapply(1:(length(pure_select_use_genes)-2),function(i) rep(pure_id[[i]],1000))),rep(pure_id[[10]],578),rep(pure_id[[11]],96))

x_train_diff1k_1k <- do.call("rbind",list(train_m_n_1000[1:1000,], train_m_n_1000[3787:4786,], train_m_n_1000[12172:13171,], train_m_n_1000[22651:23650,], train_m_n_1000[32914:33913,], train_m_n_1000[44867:45866,], train_m_n_1000[55091:56090,], train_m_n_1000[65300:66299,], train_m_n_1000[75385:76384,], train_m_n_1000[86598:87597,], train_m_n_1000[88331:88423,]))
y_train_diff1k_1k <- c(unlist(lapply(1:(length(pure_select_use_genes)-1),function(i) rep(pure_id[[i]],1000))),rep(pure_id[[11]],93))
saveRDS(list(x_train_diff1k_1k,y_train_diff1k_1k),file.path(DATA_DIR,paste(id,'normalized', genes,size,'train.rds',collapse="_")))

#--------------------------------------------------
# plot dispersion vs. mean for the genes
# --------------------------------------------------
ggplot(df,aes(mean,dispersion,col=used))+geom_point(size=0.5)+scale_x_log10()+scale_y_log10()+
  scale_color_manual(values=c("grey","black"))+theme_classic()

# --------------------------------------------------------------------------------------
# select the top 100 most variable genes
# --------------------------------------------------------------------------------------

disp_cut_off<-sort(df$dispersion_norm,decreasing=T)[100]
df$used<-df$dispersion_norm >= disp_cut_off
train_m_n_100<-train_m_n[,head(order(-df$dispersion_norm),100)]
saveRDS(list(train_m_n_100,y_train),file.path(DATA_DIR,paste(id,'normalized', genes,size,'train.rds',collapse="_")))

test_m_n_100<-test_m_n[,head(order(-df$dispersion_norm),100)]
saveRDS(list(test_m_n_100,y_test),file.path(DATA_DIR,paste(id,'normalized', genes,'test.rds',collapse="_")))

x_train_diff100_1k <- do.call("rbind",list(train_m_n_100[1:1000,], train_m_n_100[3793:4792,], train_m_n_100[12178:13177,], train_m_n_100[22657:23656,], train_m_n_100[32920:33919,], train_m_n_100[44873:45872,], train_m_n_100[55097:56096,], train_m_n_100[65306:66305,], train_m_n_100[75391:76390,], train_m_n_100[86604:87181,], train_m_n_100[87182:87277,]))
#x_train_diff1k_1k <- do.call("rbind",list(train_m_n_100[1:1000,], train_m_n_100[3787:4786,], train_m_n_100[12172:13171,], train_m_n_100[22651:23650,], train_m_n_100[32914:33913,], train_m_n_100[44867:45866,], train_m_n_100[55091:56090,], train_m_n_100[65300:66299,], train_m_n_100[75385:76384,], train_m_n_100[86598:87597,], train_m_n_100[88331:88423,]))
#x_train_1k <- do.call("rbind",list(train_m_n_1000[1:100,], train_m_n_1000[3793:3892,], train_m_n_1000[12178:12277,], train_m_n_1000[22657:22756,], train_m_n_1000[32920:33019,], train_m_n_1000[44873:44972,], train_m_n_1000[55097:55196,], train_m_n_1000[65306:65405,], train_m_n_1000[75391:75490,], train_m_n_1000[86604:86703,], train_m_n_1000[87182:87277,]))
#x_train_1k <- do.call("rbind",list(train_m_n_100[1:100,], train_m_n_100[3787:3886,], train_m_n_100[12172:12271,], train_m_n_100[22651:22750,], train_m_n_100[32914:33013,], train_m_n_100[44867:44966,], train_m_n_100[55091:55190,], train_m_n_100[65300:65399,], train_m_n_1000[75385:75484,], train_m_n_100[86598:86697,], train_m_n_100[88331:88423,]))
#y_train_1k <- c(unlist(lapply(1:(length(pure_select_use_genes)-1),function(i) rep(pure_id[[i]],100))),rep(pure_id[[11]],93))
saveRDS(list(x_train_diff100_1k,y_train_1k),file.path(DATA_DIR,paste(id,'normalized', genes,size,'train.rds',collapse="_")))
saveRDS(list(test_m_n_100,y_test),file.path(DATA_DIR,paste(id,'normalized', genes,'test.rds',collapse="_")))

#--------------------------------------------------
# plot dispersion vs. mean for the genes
# --------------------------------------------------
ggplot(df,aes(mean,dispersion,col=used))+geom_point(size=0.5)+scale_x_log10()+scale_y_log10()+
  scale_color_manual(values=c("grey","black"))+theme_classic()

