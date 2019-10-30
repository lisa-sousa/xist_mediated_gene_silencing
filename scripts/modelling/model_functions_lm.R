########################################
#this script contains all functions used in linear_model.R
#author: Lisa Barros de Andarde e Sousa
#email: lisa.barros.andrade.sousa@gmail.com
########################################
#Data: Load feature matrix and halftimes
########################################

load_data <- function(feature_matrix_file, thr_silencing_lower, thr_silencing_middle, thr_silencing_upper, high_confidence_escapees_file){
  
  print("load data...")
  
  load(feature_matrix_file)
  
  if(thr_silencing_middle == '-'){
    #XCI/escappe model
    data_set = data_set[halftime < thr_silencing_lower | halftime > thr_silencing_upper,]
    halftime = halftime[halftime < thr_silencing_lower | halftime > thr_silencing_upper]
    target = halftime
    target[target < thr_silencing_lower] = 0
    target[target > thr_silencing_upper] = 1
    target = as.factor(target)
  }else{
    #silencing dynamics model
    data_set = data_set[halftime < thr_silencing_lower | (halftime > thr_silencing_middle & halftime < thr_silencing_upper),]
    halftime = halftime[halftime < thr_silencing_lower | (halftime > thr_silencing_middle & halftime < thr_silencing_upper)]
    target = halftime
    target[target < thr_silencing_lower] = 0 #fastly silenced
    target[target > thr_silencing_middle] = 1 #slowly silenced
    target = as.factor(target)
    
    #exclude high confidence escapees 
    load(high_confidence_escapees_file)
    halftime = halftime[!rownames(data_set) %in% escapees]
    target = target[!rownames(data_set) %in% escapees]
    data_set = data_set[!rownames(data_set) %in% escapees,]
  }
  return(list(data_set, target, halftime))
}

#####################################################
#plot coefficients
#####################################################

sort_features_by_coefficients <- function(coefficients_df){
  features = row.names(coefficients_df)
  medians = rowMedians(as.matrix(coefficients_df))
  names(medians) = features
  medians = sort(medians, decreasing=TRUE)
  coefficients_df = coefficients_df[c(names(sort(medians, decreasing=TRUE))),]
  return(coefficients_df)
}

select_features_by_coefficients <- function(coefficients_df){
  medians = rowMedians(as.matrix(coefficients_df))
  names(medians) = row.names(coefficients_df)
  medians = sort(medians, decreasing=TRUE)
  selected_features = names(medians[abs(medians) > 0])
  return(selected_features)
}

plot_coefficients <- function(coefficients_df,runs){
  
  coefficients_df = coefficients_df[,complete.cases(t(coefficients_df))]
  rownames(coefficients_df) = gsub("HiC ","Hi-C ",gsub("number interactions","number",gsub("mean interaction ","",
                                  gsub("_"," ",gsub("8WG16","unphosphorylated",gsub('_GLIB_*-?[0-9]*_-?[0-9]*','',
                                  gsub('_ENCODE_[A-z]{4}_*-?[0-9]*_-?[0-9]*','',gsub('_GSE[0-9]*_-?[0-9,A-z]*_-?[0-9,A-z]*','',rownames(coefficients_df),perl=T))))))))
  
  #plot all coefficients
  sorted_coefficients = sort_features_by_coefficients(coefficients_df)
  gg_box = ggplot(melt(sorted_coefficients), aes(x=Var1,y=value)) + 
    geom_boxplot(alpha = 0.7,outlier.size=0.1,lwd=0.4) +
    scale_y_continuous(name = "coefficient)") + 
    scale_x_discrete(name = "feature") +
    theme_minimal(base_family = "Source Sans Pro") + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, angle = 90, hjust=1), axis.text.y = element_text(size=8), 
          axis.title=element_text(size=8),legend.position = c(0.85, 0.9),legend.text = element_text(size=8),legend.title = element_text(size=8)) 
  print(gg_box)
  
  #plot only features with median coefficients != 0
  selected_features = select_features_by_coefficients(coefficients_df)
  
  filtered_coefficients = sorted_coefficients[row.names(sorted_coefficients) %in% selected_features,]
  gg_box = ggplot(melt(filtered_coefficients), aes(x=Var1,y=value)) + 
    geom_boxplot(alpha = 0.7,outlier.size=0.1,lwd=0.4) +
    scale_y_continuous(name = "coefficient") + 
    scale_x_discrete(name = "feature") +
    theme_minimal(base_family = "Source Sans Pro") + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, angle = 90, hjust=1), axis.text.y = element_text(size=8), 
          axis.title=element_text(size=8),legend.position = c(0.85, 0.9),legend.text = element_text(size=8),legend.title = element_text(size=8)) 
  
  print(gg_box)
  
  #plot coefficients over runs
  colnames(filtered_coefficients) = paste("run",1:runs)
  matrix = melt(filtered_coefficients)
  
  ggplot = ggplot(matrix, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-5,5)) + # Change gradient color
    theme_minimal(base_family = "Source Sans Pro",base_size=8) +
    theme(axis.text.x = element_blank(),axis.text.y = element_text(size = 8),legend.text = element_text(size=8), 
          legend.title = element_text(size=8), axis.title=element_text(size=8),plot.background=element_blank(),panel.border=element_blank()) +
    coord_fixed(ratio=4) + 
    scale_y_discrete(name="feature",limits = rev(levels(matrix$Var1))) +
    scale_x_discrete(name="runs") 
  
  print(ggplot)
}

plot_error_rates <- function(error,class0_label,class1_label){
  
  error = as.data.frame(error[complete.cases(error),])
  error_rate = melt(error*100)
  
  gg_box = ggplot(error_rate, aes(x=variable,y=value)) + 
    geom_boxplot(alpha = 0.7,outlier.size=0.1,lwd=0.4) +
    scale_y_continuous(name = "error rate (%)",limits = c(0,70)) + 
    scale_x_discrete(name = "",labels=c("total",class0_label,class1_label)) +
    theme_minimal(base_family = "Source Sans Pro") + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, angle = 45, hjust=1), axis.text.y = element_text(size=8), 
          axis.title=element_text(size=8),legend.position = c(0.85, 0.9),legend.text = element_text(size=8),legend.title = element_text(size=8)) 
  print(gg_box)
  
}

#####################################################
#linear regression
#####################################################

linear_regression <- function(data_set,halftime,runs){
  
  data = cbind(halftime,data_set)
  
  #z-score transformation of data set
  for(i in 2:(ncol(data))){
    if(!is.factor(data[,i])){
      data[,i] = (data[,i] - mean(data[,i])) / sd(data[,i])
    }
  }
  #center halftime
  #data$halftime = data$halftime - mean(data$halftime)
  
  correlations = c()
  coefficients_df = NULL
  
  for(run in 1:runs){
    training_idx = sample(1:nrow(data),round(nrow(data)*0.8), replace=F)
    data_training = data[sort(training_idx),]
    data_test = data[-training_idx,]
    
    linear_regression_model = lm(halftime~.,data=data_training)
    coefficients = coef(linear_regression_model)
    coefficients = coefficients[complete.cases(coefficients)]
    coefficients = coefficients[2:length(coefficients)] #do not consider intercept for plot
    
    predictions = predict(linear_regression_model,newdata = data_test[,2:ncol(data_test)],type = 'response')
    
    cor = cor(data_test$halftime,predictions)
    correlations = append(correlations,cor)
    coefficients_df = cbind(coefficients_df,coefficients)
  }
  
  plot_coefficients(coefficients_df,runs)
  boxplot(correlations,ylab = 'correlations') 
  print(mean(correlations))
}

#####################################################
#logistic regression
#####################################################

logistic_regression <- function(data,training,runs){
  
  data_set = data[[1]]
  target = data[[2]]
  halftime = data[[3]]
  
  idx_y_0 = grep(0,target)
  idx_y_1 = grep(1,target)
  
  #z-score transformation of data set
  for(i in 1:(ncol(data_set))){
    if(!is.factor(data_set[,i])){
      data_set[,i] = (data_set[,i] - mean(data_set[,i])) / sd(data_set[,i])
    }
  }
  data = cbind(target,data_set)
  
  error = NULL
  coefficients_df = NULL
  
  for(run in 1:runs){
    training_idx = sort(c(sample(x=idx_y_0, size=training, replace=F),sample(x=idx_y_1, size=training, replace=F)))
    data_training = data[training_idx,]
    data_test = data[-training_idx,]
    
    remove_factor = c()
    for(f in 1:ncol(data_training)){
      if(is.factor(data_training[,f])){
        if(sum(as.numeric(data_training[,f])-1) == 0){
          remove_factor = c(remove_factor,f)
        }
      }
    }
    if(!is.null(remove_factor)){
      data_training = data_training[,-remove_factor]
      data_test = data_test[,-remove_factor]
    }

    logistic_regression_model = glm(target~.,family = binomial(link='logit'),data=data_training)
    coefficients = coef(logistic_regression_model)
    coefficients = coefficients[complete.cases(coefficients)]
    coefficients = coefficients[2:length(coefficients)] #do not consider intercept for plot
    coefficients_df = cbind(coefficients_df,coefficients)
    
    predictions = predict(logistic_regression_model,newdata = data_test,type = 'response')
    predictions = as.factor(ifelse(predictions > 0.5,1,0))
    
    misClasificError <- mean(predictions != data_test$target)
    error = append(error,misClasificError)
  }
  plot_coefficients(coefficients_df,runs) 
  boxplot(error,ylab = 'error') 
  return(error)
}

####################################
#regularized logistic regression
####################################

regularized_logistic_regression <- function(data,training,alpha,runs,runs_lambda,class0_label,class1_label){
  
  data_set = data[[1]]
  target = data[[2]]
  halftime = data[[3]]
  
  idx_y_0 = grep(0,target)
  idx_y_1 = grep(1,target)
  
  #z-score transformation of data set
  for(i in 1:(ncol(data_set))){
    if(!is.factor(data_set[,i])){
      data_set[,i] = (data_set[,i] - mean(data_set[,i])) / sd(data_set[,i])
    }
  }
  
  error = NULL
  coefficients_df = NULL
  
  for(run in 1:runs){
    
    best_lambdas = c()
    
    training_idx = sort(c(sample(x=idx_y_0, size=training, replace=F),sample(x=idx_y_1, size=training, replace=F)))
    x_training = data_set[sort(training_idx),]
    y_training = target[sort(training_idx)]
    x_test = data_set[-sort(training_idx),]
    y_test = target[-sort(training_idx)]
    
    remove_factor = c()
    for(f in 1:ncol(x_training)){
      if(is.factor(x_training[,f])){
        if(sum(as.numeric(x_training[,f])-1) == 0){
          remove_factor = c(remove_factor,f)
        }
      }
    }
    if(!is.null(remove_factor)){
      x_training = x_training[,-remove_factor]
      x_test = x_test[,-remove_factor]
    }
    
    #run the model "runs_lambda" times to get a sequence of "runs_lambda" best lambdas
    best_lambdas = foreach(i=1:runs_lambda,.combine=c,.packages='glmnet') %dopar% {
      lmodel_binom.cv <- cv.glmnet(x = data.matrix(x_training), y = y_training, family='binomial',standardize=T,alpha=alpha, nfolds=7, type.measure="class")
      lambda.min <- lmodel_binom.cv$lambda.min
      lambda.min
    }
    
    #build model with sequence of best lambdas and get the best among them
    best_lambdas = sort(best_lambdas)
    if(!length(unique(best_lambdas))==1){
      lmodel_binom.cv <- cv.glmnet(x = data.matrix(x_training), y = y_training, family='binomial',lambda=best_lambdas,standardize=T,alpha=alpha, nfolds=7, type.measure="class")
    }
    lambda.min <- lmodel_binom.cv$lambda.min

    #print(paste('optimal lambda:',round(log(lambda.min),2)))
    #plot(density(best_lambdas),col='darkgray')
    #abline(v=lambda.min,col='darkorange')
    
    #get model coefficients
    coefficients = coef(lmodel_binom.cv, s=lambda.min)
    names = rownames(coefficients)[2:nrow(coefficients)]
    coefficients = as.vector(coefficients)
    coefficients = coefficients[2:length(coefficients)]
    names(coefficients) = names
    coefficients_df = cbind(coefficients_df,coefficients)
    
    #calcuate error rate via test set
    predictions = predict(lmodel_binom.cv, newx = data.matrix(x_test), s = 'lambda.min', type = 'class')
    classification_df = data.frame(label = y_test, prediction = c(predictions), correct_predicted = c(predictions == as.character(y_test)))
    class_0_error = 1 - mean(classification_df$correct_predicted[classification_df$label == 0])
    class_1_error = 1 - mean(classification_df$correct_predicted[classification_df$label == 1])
    total_error = 1 - mean(classification_df$correct_predicted)
    
    error = rbind(error,c(total_error,class_0_error,class_1_error))
  }
  
  plot_error_rates(error,class0_label,class1_label)
  plot_coefficients(coefficients_df,runs)
  
  return(error)
}
