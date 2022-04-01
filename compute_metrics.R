suppressMessages(require(caret))
suppressMessages(require(pROC))
suppressMessages(require(stringr))
suppressMessages(require(RColorBrewer))
suppressMessages(require(dendextend))
suppressMessages(require(gplots))
suppressMessages(require(argparse))

parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
    dest="verbose", help="Print little output")
parser$add_argument("--f", "--file", type="character",  
    help = "File name of model")
parser$add_argument("--d", "--data_dir", type="character",  
    help = "Directory that contains data")
parser$add_argument("--r", "--results_dir", type="character",  
    help = "Output directory for models")
args <- parser$parse_args()

# ----------------------
# parse variables
# ----------------------

FILE <- args$file
RESULTS_DIR = args$results_dir
DATA_DIR = args$data_dir

# ----------------------
# get multiclass roc
# ----------------------
get_roc <- function(dir,file,model) {
  if (model %in% c("svmLinear","svmRadial")) {
    res <- readRDS(file.path(dir,file))
    pred <- res[[2]]
    pred_prob <- as.matrix(data.frame(attr(pred, "probabilities")))
    colnames(pred_prob) <- c("CD34+", "CD56+ NK" , "CD4+/CD45RA+/CD25- Naive T","CD4+/CD25 T Reg",
                             "CD8+/CD45RA+ Naive Cytotoxic","CD4+/CD45RO+ Memory" ,"CD8+ Cytotoxic T",
                             "CD19+ B","CD4+ T Helper2","CD14+ Monocyte","Dendritic")
    labels <- factor(annot$celltype)
    roc.multi <- multiclass.roc(labels,pred_prob)
  } else if (model == "cvglmnet") {
    res <- readRDS(file.path(dir,file))
    pred_prob <- as.matrix(data.frame(res[[2]]))
    colnames(pred_prob) <- levels(labels)
    labels <- factor(annot$celltype)
    roc.multi <- multiclass.roc(labels,pred_prob)
  } else if (model == "rf") {
    res <- readRDS(file.path(dir,file))
    pred_prob <- as.matrix(data.frame(res[[2]]))
    colnames(pred_prob) <- levels(labels)
    labels <- factor(annot$celltype)
    roc.multi <- multiclass.roc(labels,pred_prob)
  }
  return(roc.multi)
}

# ----------------------
# get auc for all models
# ----------------------
get_model_auc <- function(dir) {
  files = list.files(dir)
  model_auc <- list()
  for (file in files) {
    print(file)
    id <- str_replace_all(file,c("pure_use_genes_" = "",".rds" = "","_default"="","_scaled"="","_normalized"="","normalized_"=""))
    model = rev(unlist(str_split(id,"_")))[1]
    roc.multi <- get_roc(file,model)
    model_auc[[id]] <- unlist(str_split(roc.multi$auc,':'))
  }

  model_auc <- data.frame(auc = do.call(rbind,model_auc))
  model_auc$dataset <- rownames(model_auc)
  return(model_auc)
}

# ----------------------
# plot auc for all pairs of classes 
# ----------------------
plt_pair_auc <- function(file,model) {
  roc.multi <- get_roc(file,model)
  srocs <- roc.multi[['rocs']]
  for (j in 1:length(rocs)) {
    rs <- rocs[[j]]
    print(paste0(rs[[1]]$levels[1]," vs ",rs[[1]]$levels[2],": "))
    print(auc(rs[[1]]))
    plot.roc(rs[[1]])
    sapply(2:length(rs),function(j) lines.roc(rs[[j]],col=j))  
  }
  aucs <- data.frame(do.call("rbind",lapply(1:length(rocs),function(i) c(as.character(rocs[[i]][[1]]$levels),unlist(str_split(auc(rocs[[i]][[1]]),':'))))),stringsAsFactors = F)
  colnames(aucs) <- c("Cell Type 1","Cell Type 2","AUC")
  aucs$AUC <- as.numeric(aucs$AUC)
  aucs$AUC <- ifelse(aucs$AUC < 0.5, 1-aucs$AUC, aucs$AUC)
  #write.csv(aucs,file.path(RESULTS_DIR,'bestmodel_auc.csv'))
  pl <- ggplot(aucs, aes(x = `Cell Type 1`, y = `Cell Type 2`, fill = AUC)) +
    geom_tile() + 
    scale_fill_gradient(low="green", high="red", limits=c(0.5, 1)) + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(aucs)
}


# ----------------------
# plot random forest importance 
# ----------------------
plt_rf_imp <- function(dir,file,labels) {
  rf <- readRDS(file.path(dir,file))[[1]]
  imp <- data.frame(rf$importance)[1:11]
  top10_imp <- unlist(data.frame(sapply(imp, function(x) head(row.names(imp)[order(x, decreasing = TRUE)], 10)),stringsAsFactors = F))
  imp <- imp[row.names(imp) %in% top10_imp,]
  colnames(imp) <- levels(labels)

  colors <- brewer.pal(3, "Set2")

  dend <- dist(t(imp), method = 'euclidean') %>%  
    hclust(method = "complete") %>% 
    as.dendrogram %>%
    hang.dendrogram(hang_height=0.1) %>%
    set("labels_cex", 0.5) %>%
    color_branches(k=3, col=colors)

  par(mar=c(8,8,4,14))
  # plot heatmap 

  hmcols <- colorpanel(3000, "white","black")

  heatmap.2(as.matrix(t(imp)),  
            srtCol = 90,
            key = TRUE,
            trace="none",
            dendrogram = "row",
            Rowv = dend,
            Colv = "NA", 
            col=hmcols,
            margins=c(8,14),
            cexCol=1,
            cexRow=1)     
  return(imp)
}

# ----------------------
# compute and plot metrics
# ----------------------

model <- rev(unlist(str_split(str_replace(FILE,".rds",""),"_")))[1]
model_aucs <- get_model_auc(RESULTS_DIR)
pair_aucs <- plt_pair_auc(FILE,model)

annot <- read.csv(file.path(DATA_DIR,'68k_pbmc_barcodes_annotation.tsv'),sep='\t')
labels <- annot$celltype
plt_rf_imp(RESULTS_DIR,FILE,labels)
