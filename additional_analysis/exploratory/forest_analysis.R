library(randomForest)

###################################
#load data
###################################

load('./work_stuff/projects/xist_epigenetics/computing/feature_matrix/hdac3/promoter_matrix_original_normRFixed.RData')

thr_silencing_lower = 0.9
thr_silencing_middle = "-"
thr_silencing_upper = 1.6
x = 10
mtry = 8

if(thr_silencing_middle == '-'){
  #silenced vs not silenced
  data_set = data_set[halftime < thr_silencing_lower | halftime > thr_silencing_upper,]
  halftime = halftime[halftime < thr_silencing_lower | halftime > thr_silencing_upper]
  target = halftime
  target[target < thr_silencing_lower] = 0 #fastly silenced
  target[target > thr_silencing_upper] = 1 #slowly silenced
  target = as.factor(target)
}else{
  #early vs late
  data_set = data_set[halftime < thr_silencing_lower | (halftime > thr_silencing_middle & halftime < thr_silencing_upper),]
  halftime = halftime[halftime < thr_silencing_lower | (halftime > thr_silencing_middle & halftime < thr_silencing_upper)]
  target = halftime
  target[target < thr_silencing_lower] = 0 #fastly silenced
  target[target > thr_silencing_middle] = 1 #slowly silenced
  target = as.factor(target)
  
  #exclude high confidence escapees and marks escapees
  load('/project/lncrna/Xist/data_lisa/annotation_files/escapees/escapees.RData')
  escapees_marks = read.table('/project/lncrna/Xist/data_lisa/annotation_files/escapees/metadata/escapees_2015_marks.txt')
  escapees_marks = as.character(escapees_marks$V1)
  escapees = unique(c(escapees,escapees_marks))
  
  halftime = halftime[!rownames(data_set) %in% escapees]
  target = target[!rownames(data_set) %in% escapees]
  data_set = data_set[!rownames(data_set) %in% escapees,]
}

load('./work_stuff/projects/xist_epigenetics/computing/silenced_vs_not_silenced/features_promoter_thr_0.9_-_1.6_selected_features.RData')
selected_features = unique(c(selected_feature_imp1[1:x],selected_feature_imp2[1:x]))
data_set = data_set[colnames(data_set)%in%selected_features]
head(data_set)

###################################
#Random Forest model - get best tree
###################################

#training parameters
ntree = 1000 #500 trees suggested as a good estimation for ntree (in paper)
runs = 500

#set sampsize
smaller_set = round(min(c(sum(target==0),sum(target==1))))
test = 10
training = round(min(c(smaller_set/100*80,smaller_set-test)))
sampsize=c(training,training)

#results that should be saved
proximity = NULL
predictions_class = NULL
predictions_vote = NULL

best_error_rate = 1
best_tree = NULL

for(i in 1:runs){
  class_forest = randomForest(data_set, target, type="classification", ntree=ntree, mtry=mtry, replace=FALSE, sampsize=sampsize, importance=TRUE, proximity = TRUE)
  
  error_rate = class_forest$err.rate[ntree,1]
  class_0_error = class_forest$err.rate[ntree,2]
  class_1_error = class_forest$err.rate[ntree,3]
  difference = abs(class_0_error - class_1_error)
  
  if(error_rate < best_error_rate & difference < 0.02){
    print("new best tree:")
    print(error_rate)
    print(class_0_error)
    print(class_1_error)
    best_error_rate = error_rate
    proximity = class_forest$proximity
    predictions_class = as.numeric(predict(class_forest,data_set,type = 'class'))-1
    predictions_vote = predict(class_forest,data_set,type = 'vote')[,1]
    best_tree = class_forest
  }
  
}


###################################
#count all cooccurences of features
###################################

cooccurences = as.data.frame(matrix(0,nrow = 15, ncol = 15))
colnames(cooccurences) = colnames(data_set)
rownames(cooccurences) = colnames(data_set)

for(forest in 1:1000){
  variables = as.vector(sort(unique(best_tree$forest$bestvar[,forest])))
  variables = variables[!variables==0]
  
  for(i in 1:length(variables)){
    cooccurences[variables[i],variables] = cooccurences[variables[i],variables] + 1
  }
}

pdf("coccurences.pdf",width = 15,height = 15)
corrplot(as.matrix(cooccurences), order='hclust', hclust.method='ward.D2', col=colorRampPalette(c('blue','white','red'))(200), diag=F, cl.cex=1.3, tl.cex=1.3, tl.col = "black",is.corr = FALSE,bg='white', addgrid.col=NA)
dev.off()

###################################
#plot distribution of features at each node level
###################################


pdf("feature_distribution_at_nodes.pdf",width = 10,height = 10)
for(i in 1:10){
  variables = as.vector(best_tree$forest$bestvar[i,])
  variables = tabulate(variables[!variables==0])
  names(variables) = colnames(data_set)
  par(mar=c(14,2,2,2),oma=c(1,1,1,1))
  barplot(variables, las=2, ylim = c(0,500))
}
dev.off()




