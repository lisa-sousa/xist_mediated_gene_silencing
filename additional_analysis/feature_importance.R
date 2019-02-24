library(matrixStats)
library(gplots)

load("./work_stuff/projects/xist_epigenetics/computing/data/modelling/model/silencing_dynamics_model/results_thr_0.5_0.9_1.3/random_forest_model_all_features.RData")

meanImpClass0 = data.frame(feature=row.names(random_forest_model_all_features[[2]]),
                           meanImpClass0 = rowMedians(random_forest_model_all_features[[2]]))
meanImpClass0$meanImpClass0[meanImpClass0$meanImpClass0 < 0] = 0

meanImpClass1 = data.frame(feature=row.names(random_forest_model_all_features[[3]]),
                           meanImpClass1 = rowMedians(random_forest_model_all_features[[3]]))
meanImpClass1$meanImpClass1[meanImpClass1$meanImpClass1 < 0] = 0

importance = merge(meanImpClass0,meanImpClass1,by="feature")
importance$mean = rowMeans(importance[,2:3])
importance = importance[order(importance$mean,decreasing = T),]
row.names(importance) = importance$feature
importance = as.matrix(importance[,2:3])

breaks=c(min(importance),0.01,0.2,0.3,0.4,0.5,1,2,3,4,5,10,max(importance))
mycol = c("white","#E2E2E2","#DBDCE0","#CFD0DA","#BEC1D4","#ABB0CC","#959CC3","#7D87B9","#6371AF","#4359A7","#023FA5","darkblue")

pdf(file = "./work_stuff/projects/xist_epigenetics/computing/plots/additional_analysis/feature_importance.pdf",height = 20, width = 10)
par(mfrow=c(1,1),mar=c(1,1,1,1),oma=c(7,3,3,5))
heatmap.2(importance,col=mycol,breaks=breaks,Colv=NA,Rowv=NA,scale='none',cexCol=1.8,main='feature importance',trace='none',dendrogram='none')
dev.off()
