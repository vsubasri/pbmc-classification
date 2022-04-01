# ----------------------
# run random forest with default parameters
# ----------------------
rf_default <- function(dir,id,genes,size,model,x_train,y_train,x_test) {
	id <- paste0(id,genes,size,model,collapse="_")
	file <- file.path(dir,paste0(id,'_default.rds'))
	if (file.exists(file)) {
  		rf_res <- readRDS(file)
		rf <- rf_res[[1]]  
	} else {
	rf <- randomForest(x=x_train,
                     y=y_train,
                     type = "class",
                     importance=TRUE,
                     do.trace=100)
	}

	y_pred = predict(rf, newdata = x_test, type="prob")
	y_class = predict(rf, newdata = x_test, type="response")
	saveRDS(list(rf,y_pred,y_class),file)
	return(list(rf,y_pred,y_class))
}

# ----------------------
# run random forest with parameter tuning
# ---------------------- 
rf_tuned <- function(dir,id,genes,size,model,x_train,y_train,x_test) {
	id <- paste0(id,genes,size,model,collapse="_")
	res <- tuneRF(x_train, y_train,
              stepFactor = 1.5,
              plot = TRUE,
              ntreeStart = 3,
              ntreeTry = 500,
              trace = TRUE,
              improve = 0.05)

	      # fit tuned rf
	file <- file.path(dir,paste0(id,'_tuned.rds'))
	tuned_rf <- randomForest(x=x_train,
                         y=y_train,
                         type = "class",
                         mtry = res[1],
                         importance=TRUE,
                         do.trace=100 )
	y_pred = predict(tuned_rf, newdata = x_test, type = "prob")
	y_class = predict(tuned_rf, newdata = x_test, type = "response")
	saveRDS(list(tuned_rf,y_pred,y_class),file)
	return(list(tuned_rf,y_pred,y_class))

}

# ----------------------
# run svm with linear kernel
# ----------------------
svm_linear <- function(dir,id,genes,size,x_train,y_train,x_test) {
	costs = c(0.01,0.1,1,10)
	for (c in costs) {
  		cat(paste0("Tuning c: ",c))
		file = file.path(dir,paste(id,genes,size,c,'svmLinear.rds',collapse="_"))
  		if (file.exists(file)) {
    			model = readRDS(file)
	} else { 
    		model <- svm(x_train, 
                     y_train,
                     kernal="linear",
                     cost = c,
                     scale = TRUE,
                     probability = TRUE) 
    		y_pred <- predict(model, x_test,decision.values = TRUE, probability = TRUE)
    		y_class <- predict(model, x_test)
    		saveRDS(list(model,y_pred,y_class),file)
		return(list(model,y_pred,y_class))
  	}  
}

# ----------------------
# run svm with rbf kernel 
# ----------------------
svm_radial <- function(dir,id,genes,size,x_train,y_train,x_test) {
	costs = c(0.01,0.1,1,10)
	gammas = c(0.5,1,1.5)
	for (g in gammas) {
  		for (c in costs) {
    			cat(paste0("Tuning c: ",c," g: ",g))
    			file = file.path(dir,paste(id,genes,size,c,g,'svmLinear.rds',collapse="_"))
    			if (file.exists(file)) {
      				model = readRDS(file)
    			} else { 
      				model <- svm(x_train, 
                       			y_train,
                       			kernal="radial",
                       			scale = TRUE,
                       			probability = TRUE,
                       			cost = c,
                       			gamma = g) 
      				y_pred <- predict(model, x_test,decision.values = TRUE, probability = TRUE)
      				y_class <- predict(model, x_test)
      				saveRDS(list(model,y_pred,y_class),file)
    			}
  		}
	}
	return(list(model,y_pred,y_class))
}

# ----------------------
# run glm
# ----------------------
glm <- function(dir,id,genes,size,x_train,y_train,x_test) {
	file = file.path(dir,paste(id,genes,size,'cvglmnet.rds',collapse="_"))
	if (file.exists(file)) {
  		cvglmnet_res <- readRDS(file)
  		cvglmnet <- cvglmnet_res[[1]]  
	} else {
  	cvglmnet <- cv.glmnet(x=x_train,
                        y=y_train,
                        family = "multinomial",
                        type.measure = "class",
                        nfolds = 5,
                        trace.it=1,
                        keep=TRUE,
                        parallel = TRUE)
	}
	y_pred = predict(cvglmnet, newx = x_test, s="lambda.min",type = "response")
	y_class = predict(cvglmnet, newx = x_test, s="lambda.min")
	saveRDS(list(cvglmnet,y_pred,y_class),file)
	return(list(cvglmnet,y_pred,y_class))
}
