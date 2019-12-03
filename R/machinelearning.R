#!/usr/bin/env R

# Machine learning functions

# Run machine learning
#' @export
runLasso <- function(seset, seed=2019){
  # runLasso
  # Fit a model using penalized regression with lasso
  # Arguments:
  # * sese: Valid summarized experiment object
  # * seed: (int) set seed for randomization
  # Returns:
  # * resultslist (list) : Results of lasso fit
  set.seed(seed)

  gene.names = as.character(rownames(rowData(seset)))
  var.classifier = seset$deg.risk
  df = t(assay(seset))
  train.names = colnames(assay(seset[,seset$exptset.seahack=="train"]))
  test.names = colnames(assay(seset[,seset$exptset.seahack=="test"]))
  response <- var.classifier
  predictors <- gene.names
  y <- factor(response); names(y) <- colnames(assay(seset)) # response var obj
  x = df[,colnames(df) %in% predictors] # genes of interest
  contrast <- contrasts(y)
  grid <- 10^ seq(10,-2, length=100)
  standardize = FALSE
  fit <- glmnet::glmnet(x[train.names,], y[train.names], family = "binomial", alpha=1,
                standardize = standardize, lambda = grid, intercept = FALSE)

  # use cross-validation on the training model.CV only for lambda
  cv.fit <- glmnet::cv.glmnet(x[train.names,], y[train.names], family = "binomial",
                      type.logistic="modified.Newton", standardize = standardize,
                      lambda = grid, alpha=1, nfolds = length(train.names), #LOOCV
                      type.measure = "class", intercept = FALSE)
  #Select lambda min.
  lambda.min <- cv.fit$lambda.min
  #predict the classes
  pred.class <- predict(fit, newx = x[test.names,], type="class", s=lambda.min)
  #find the test error
  tab <- table(pred.class,y[test.names])
  testError <- mean(pred.class != y[test.names]) #how many predicted classes were incorrect
  #Fit the full dataset.
  final <- glmnet::glmnet(x, y,family = "binomial", standardize = standardize,
                  lambda = grid, alpha = 1, intercept = FALSE)
  #Extract the coefficients
  coef <- predict(final, type="coefficients", s=lambda.min)
  idx <- which(!as.numeric(coef)==0)
  nonZero <- coef[idx,]
  #Results
  resultslist <- list(train.names, test.names, contrast, fit,
                      cv.fit, tab, testError, final, nonZero, seed)
  names(resultslist) <- c("training.set", "testing.set","contrast", "train.fit",
                          "cv.fit", "confusionMatrix","test.error", "final.model",
                          "nonzero.coef", "seed")
  return(resultslist)
}

#' @export
runSVM <- function(seed,kerneltype="linear",seset, weightfilt=FALSE){
  # credit base code: Sean Maden
  # run SVM optimization
  # Arguments
  #   * seed : set seed (int) for randomization
  #   * kerneltype : (str) valid kernel type class for SVM (e.g. 'linear', 'radial', etc.)
  #   * seset : summarized expirment object with both test and training set data.
  #       * ndtr : training dataset (excluding sample classes)
  #       * ndtr.classes : classes for training sampels (vector) with 1:1 correspondence
  #       with trainset rows
  #       * ndte : test data (data frame or matrix), excluding classes
  #       * ndte.classes : classes for test samples (vector), with 1:1 row:pos correspondence
  #   * weightfilt : (FALSE or numeric float) top fraction weights to use in model
  #       (if FALSE, then all weights used)
  # Returns
  #   * rl (list) : list containing model fitted, predictions, and performacne metrics
  rl <- list(); str.options <- ""
  set.seed(seed)

  #training/testing sets
  ndtr <- t(assay(seset[,seset$exptset.seahack=="train"]))
  ndtr.classes <- seset[,seset$exptset.seahack=="train"]$deg.risk

  ndte <- t(assay(seset[,seset$exptset.seahack=="test"]))
  ndte.classes <- seset[,seset$exptset.seahack=="test"]$deg.risk
  # train svm model
  svm_model <- e1071::svm(as.factor(ndtr.classes)~.,
                   data=ndtr,
                   method="C-classification",
                   kernel=kerneltype,
                   probability=TRUE)
  weightsvect <- ndtr.weights <- t(svm_model$coefs) %*% svm_model$SV
  if(weightfilt){
    str.options <- c(str.options,paste0("weight filt = ",weightfilt))
    # order training data on relative weights
    ndtr.weightsort <- ndtr[,rev(order(abs(ndtr.weights)))]
    # select only top proportion weights
    nweight.col = round(ncol(ndtr.weightsort)*weightfilt,0)
    ndtr.weightfilt <- ndtr.weightsort[,c(1:nweight.col)]
    str.options <- c(str.options,paste("cols_retained:",colnames(ndtr.weightfilt),collapse=";"))
    # redefine training set, rerun SVM optimization
    ndtr <- ndtr.weightfilt
    svm_model <- e1071::svm(as.factor(ndtr.classes)~.,
                     data=ndtr,
                     method="C-classification",
                     kernel=kerneltype,
                     probability=TRUE)
  } else{
    str.options <- c(str.options,"no weight filt")
  }

  #training and test set predictions.
  pred_train <- predict(svm_model, ndtr, decision.valuesq = TRUE, probability=TRUE)
  pred_test <- predict(svm_model, ndte, decision.values = TRUE, probability=FALSE)
  pred_test2 <- predict(svm_model, ndte, decision.values = TRUE, probability=TRUE)
  #find the test errors
  tab <- table(pred = pred_test, true = ndte.classes)
  testError <- mean(pred_test != ndte.classes) #how many predicted classes were incorrect
  log.loss <- MLmetrics::LogLoss(y_pred = attr(pred_test2, which="probabilities")[,1], y_true = as.numeric(ndte.classes))
  # get performance metrics
  pred <- prediction(as.numeric(attr(pred_test,"decision.values")),ndte.classes)
  perf <- ROCR::performance(pred,"tpr","fpr")
  ppred <- pred_test[pred_test==1]
  tppred <- ndte.classes[pred_test==1]
  ppred <- as.numeric(as.character(ppred))
  testprec <- length(ppred[ppred==tppred])/length(ppred) # test precision
  rposi <- ndte.classes==1
  rtpred <- ndte.classes[rposi]
  rppred <- pred_test[rposi]
  rppred <- as.numeric(as.character(rppred))
  testrec <- length(rppred[rppred==1])/length(rppred) # test recall

  # return model, pred's, and performance metrics
  rl <- list(str.options,
             svm_model,
             weightsvect,
             pred_train,
             pred_test,
             pred_test2,
             perf,
             tppred,
             testprec,
             testrec,
             tab,
             testError,
             log.loss)
  names(rl) <- c("options_string",
                 "svm_model",
                 "weightsvect",
                 "predictions_train",
                 "predictions_test",
                 "predictions_test2",
                 "performance_test",
                 "TPR_test",
                 "precision_test",
                 "recall_test",
                 "confusionMatrix",
                 "testError",
                 "log.loss")
  return(rl)

}

# Importance Functions
#' @export
impLasso <- function(df, classes, trainindices, seed=2019){
  # df : data frame to parse, rownames = classifier groupings, colnames = feature ids
  set.seed(seed)
  var.classifier <- response <- as.character(classes)
  y <- factor(response); names(y) <- rownames(df) # response var obj
  x = df # genes of interest
  contrast <- contrasts(y)
  grid <- 10^ seq(10,-2, length=100)

  # use cross-validation on the training model.CV only for lambda
  message("performing cross-validation...")
  cv.fit <- glmnet::cv.glmnet(x[trainindices,], y[trainindices], family = "binomial",
                      type.logistic="modified.Newton", standardize = FALSE,
                      lambda = grid, alpha=1, nfolds = length(trainindices), #LOOCV
                      type.measure = "class", intercept = FALSE)
  #Select lambda min.
  message("selecting lambda min...")
  lambda.min <- cv.fit$lambda.min

  #Fit the full dataset.
  message("fitting whole dataset")
  lassofit <- glmnet::glmnet(x, y, family = "binomial", standardize = FALSE,
                     lambda = grid, alpha = 1, intercept = FALSE)

  #Extract the coefficients
  lassoimp <- predict(lassofit, type="coefficients", s=lambda.min)
  lassoimp <- lassoimp[2:nrow(lassoimp),1]
  return(lassoimp)
}

#' @export
impRF <- function(df, classes, ntrees=100, seed=2019){
  set.seed(seed)
  class <- as.numeric(classes)
  rffit <- randomForest::randomForest(class ~ ., data = as.matrix(df), ntree = ntrees,proximity = TRUE)
  rfimp <- randomForest::importance(rffit)[,1]
  return(rfimp)
}

#' @export
impXGB <- function(df, classes, seed=2019){
  set.seed(2019)
  message("fitting xgboost model...")
  xgbfit <- xgboost::xgboost(data = df, label = classes, max_depth = 2,
                    eta = 1, nthread = 2, nrounds = 2, objective = "binary:logistic")
  message("getting xgboost importances...")
  xgbimp <- xgboost::xgb.importance(feature_names = colnames(df), model = xgbfit)
  message("reformatting xgb importances...")

  xgbimp.format <- c(rep(0, ncol(df)))
  names(xgbimp.format) <- colnames(df)
  for(f in 1:nrow(xgbimp)){
    xgbimp.format[which(colnames(df)==as.character(xgbimp[f,1]))] <- as.numeric(xgbimp[f,2])
  }
  xgbimp.format <- as.numeric(xgbimp.format)
  names(xgbimp.format) <- colnames(df)

  return(xgbimp.format)
}

#' @export
impSVM <- function(df, classes, seed=2019){
  set.seed(2019)
  message("fitting svm model...")
  svmfit <- e1071::svm(as.factor(classes)~.,
                data=df,
                method="C-classification",
                kernel="linear")
  svmimp <-  t(svmfit$coefs) %*% svmfit$SV
  message("reformatting importances for output...")
  svmimp.format <- svmimp[1,]
  names(svmimp.format) <- colnames(svmimp)
  return(svmimp.format)
}

#' @export
impCML <- function(df, classes, trainindices, seed, algo.opt=c("lasso","rf","svm","xgb"),
                   standtable=FALSE, ranksummary="median"){
  # impCML
  # Get consensus importance ranks from disparate algorithms.
  # Arguments
  # * df (matrix) : data table (rows = instances, columns = features)
  # * classes (character): categorizations of instances
  # * algo.opt (list): List of valid algorithms to use for consensus
  # * seed (int): Set the seed for reproducibility
  # * trainindices (numeric, optional): Vector of dataset (rows) corresponding to training sample subsert. Only used for Lasso cross validation step.
  # * ranksummary (string): Either "median" or "mean", the operation used to calculate consensus rank.
  # Returns:
  # * Consensus rank, optionally a standard output table of importances for selected algorithms (if standtable==TRUE)
  print("getting importances...")
  implist <- list()
  implabellist <- c()
  if("lasso" %in% algo.opt){
    imp.lasso <- impLasso(df=as.matrix(df),trainindices=trainindices,classes=classes)
    implist[["lasso"]] <- imp.lasso
    implabellist <- c(implabellist, "lasso")
  }
  if("rf" %in% algo.opt){
    imp.rf <- impRF(df=df, classes=classes)
    implist[["rf"]] <- imp.rf
    implabellist <- c(implabellist, "rf")
  }
  if("xgb" %in% algo.opt){
    imp.xgb <- impXGB(df=df, classes=classes)
    implist[["xgb"]] <- imp.xgb
    implabellist <- c(implabellist, "xgb")
  }
  if("svm" %in% algo.opt){
    imp.svm <- impSVM(df=df, classes=classes)
    implist[["svm"]] <- imp.svm
    implabellist <- c(implabellist, "svm")
  }
  # match importance feature label ordering
  if(length(implist)>1){
    for(i in 1:(length(implist)-1)){
      implist[[i]] <- implist[[i]][order(match(names(implist[[i]]),names(implist[[i+1]])))]
    }
  }

  # get importance ranks and rankvars
  message("computing rank importances...")
  impranklist <- list()
  for(i in 1:length(implist)){
    # always compute absolute ranks
    impranklist[[names(implist)[i]]] <- rank(abs(implist[[i]]))
  }

  medrank <- c()
  meanrank <- c()
  for(r in 1:length(impranklist[[1]])){
    rvalr <- c()
    for(i in 1:length(impranklist)){
      rvalr <- c(rvalr, impranklist[[i]][r])
    }
    medrank <- c(medrank, median(rvalr))
    meanrank <- c(meanrank, mean(rvalr))
  }
  if(ranksummary=="median"){
    lr <- medrank
  } else{
    lr <- meanrank
  }

  if(standtable==TRUE){
    # form the stdout table
    message("forming standard table...")
    st <- matrix(nrow=length(impranklist[1]), ncol=0)
    for(r in 1:length(implist)){
      st <- cbind(st, matrix(implist[r], ncol=1))
      colnames(st)[r] <- implabellist[r]
    }
    st <- cbind(st, matrix(medrank, ncol=1)); colnames(st)[ncol(st)] <- "medrank"
    st <- cbind(st, matrix(meanrank, ncol=1)); colnames(st)[ncol(st)] <- "meanrank"

    lr <- list("ranksummary"=lr,
               "standtable"=st)
  }
  message("completed all tasks! Returning...")
  return(lr)
}

#' @export
impBorutaCML <- function(x=x, y=y, seed, algo.opt=c("lasso","rf","svm","xgb"),
                   ranksummary="median"){
  # impCML
  # Get consensus importance ranks from disparate algorithms.
  # Arguments
  # * df (matrix) : data table (rows = instances, columns = features)
  # * classes (character): categorizations of instances
  # * algo.opt (list): List of valid algorithms to use for consensus
  # * seed (int): Set the seed for reproducibility
  # * trainindices (numeric, optional): Vector of dataset (rows) corresponding to training sample subsert. Only used for Lasso cross validation step.
  # * ranksummary (string): Either "median" or "mean", the operation used to calculate consensus rank.
  # Returns:
  # * Consensus rank, optionally a standard output table of importances for selected algorithms (if standtable==TRUE)

  df <- x
  classes <- y

  print("getting importances...")
  implist <- list()
  implabellist <- c()
  if("lasso" %in% algo.opt){
    imp.lasso <- impLasso(df=as.matrix(df), trainindices=trainindices, classes=classes)
    implist[["lasso"]] <- imp.lasso
    implabellist <- c(implabellist, "lasso")
  }
  if("rf" %in% algo.opt){
    imp.rf <- impRF(df=as.matrix(df), classes=classes)
    implist[["rf"]] <- imp.rf
    implabellist <- c(implabellist, "rf")
  }
  if("xgb" %in% algo.opt){
    imp.xgb <- impXGB(df=as.matrix(df), classes=classes)
    implist[["xgb"]] <- imp.xgb
    implabellist <- c(implabellist, "xgb")
  }
  if("svm" %in% algo.opt){
    imp.svm <- impSVM(df=as.matrix(df), classes=classes)
    implist[["svm"]] <- imp.svm
    implabellist <- c(implabellist, "svm")
  }
  # match importance feature label ordering
  if(length(implist)>1){
    for(i in 1:(length(implist)-1)){
      implist[[i]] <- implist[[i]][order(match(names(implist[[i]]),names(implist[[i+1]])))]
    }
  }

  # get importance ranks and rankvars
  message("computing rank importances...")
  impranklist <- list()
  for(i in 1:length(implist)){
    # always compute absolute ranks
    impranklist[[names(implist)[i]]] <- rank(abs(implist[[i]]))
  }

  medrank <- c()
  meanrank <- c()
  for(r in 1:length(impranklist[[1]])){
    rvalr <- c()
    for(i in 1:length(impranklist)){
      rvalr <- c(rvalr, impranklist[[i]][r])
    }
    medrank <- c(medrank, median(rvalr))
    meanrank <- c(meanrank, mean(rvalr))
  }
  if(ranksummary=="median"){
    lr <- medrank
  } else{
    lr <- meanrank
  }
  message("completed all tasks! Returning...")
  return(lr)
}


