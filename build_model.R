suppressMessages(require(stringr))
suppressMessages(require(parallel))
suppressMessages(require(doParallel))
suppressMessages(require(dplyr))
suppressMessages(require(magrittr))
suppressMessages(require(caret))
suppressMessages(require(stringr))
suppressMessages(require(e1071))
suppressMessages(require(randomForest))
suppressMessages(require(glmnet))
suppressMessages(require(Matrix))
suppressMessages(require(argparse))
source('models.R')
source('util.R')

cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

parser <- ArgumentParser()

parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
    dest="verbose", help="Print little output")
parser$add_argument("-g", "--genes", type="integer", default=1000, 
    help="Number of DEGs")
parser$add_argument("-s", "--size", type="integer", default=1000,
    help="Number of each cells for each cell type to include in the training set")
parser$add_argument("-m", "--model", type="character", default="rf",
    help="Model to train for cell type classficiation")
parser$add_argument("--i", "--id", type="character",  
    help = "Prefix for filename")
parser$add_argument("--d", "--data_dir", type="character",  
    help = "Directory that contains preprocessed data")
parser$add_argument("--r", "--results_dir", type="character",  
    help = "Output directory for models")

args <- parser$parse_args()

# -------------------------------------
# parse variables
# -------------------------------------

genes <- args$genes
size <- args$size
id <- args$id
data_dir <- args$data_dir 
results_dir <- args$results_dir 
model <- args$model

# -------------------------------------
# load data 
# -------------------------------------

data <- load_data(data_dir,id,genes,size)
x_train <- data[[1]] ; y_train <- data[[2]]
x_test <- data[[3]] ; y_test <- data[[4]]

# -------------------------------------
# fit model
# -------------------------------------

if (model == "rf_default") {
    rf_default(results_dir,id,genes,size,model,x_train,y_train,x_test)
}

if (model == "rf_tuned") {
    rf_tuned(results_dir,id,genes,size,model,x_train,y_train,x_test) 
}

if (model == "svm_linear") {
    svm_linear(results_dir,id,genes,size,x_train,y_train,x_test)
}

if (model == "svm_radial") {
    svm_radial(results_dir,id,genes,size,x_train,y_train,x_test)
}

if (model == "glm") {
    glm(results_dir,id,genes,size,x_train,y_train,x_test)
}
