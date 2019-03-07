library(gplots)


thrshold_table = data.frame(thr_lower = c(0.9,0.9,0.9,0.9,0.9,1.0,1.0,1.0,1.0,1.1,1.2,1.3,1.4), 
                            thr_middle = c("-","-","-","-","-","-","-","-","-","-","-","-","-"),
                            thr_upper = c(1.5,1.6,1.7,1.8,2.0,1.5,1.6,1.7,1.8,1.6,1.6,1.6,1.6))
thrshold_table = data.frame(thr_lower = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.7,0.7,0.7), 
                            thr_middle = c(0.7,0.7,0.7,0.7,0.8,0.8,0.8,0.8,0.9,0.7,0.8,0.8,0.8,0.8,0.9,0.8,0.8,0.9),
                            thr_upper = c(1.1,1.2,1.3,1.4,1.1,1.2,1.3,1.4,1.4,1.3,1.1,1.2,1.3,1.4,1.4,1.3,1.4,1.4))
thrshold_table = data.frame(thr_lower = c(0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.7), 
                            thr_middle = c(0.7,0.8,0.8,0.8,0.9,0.8,0.8,0.8,0.8,0.9,0.9),
                            thr_upper = c(1.4,1.2,1.3,1.4,1.4,1.1,1.2,1.3,1.4,1.4,1.4))


features_class1_top10 = c()
features_class2_top10 = c()

for(i in 1:nrow(thrshold_table)){
  #file = paste("/project/lncrna/Xist/data_lisa/classification/silenced_vs_not_silenced_narrow/features_promoter_thr_",thrshold_table$thr_lower[i],"_",thrshold_table$thr_middle[i],"_", thrshold_table$thr_upper[i],"_selected_features.RData",sep="")
  file = paste("/project/lncrna/Xist/data_lisa/classification/early_vs_late_narrow/features_promoter_thr_",thrshold_table$thr_lower[i],"_",thrshold_table$thr_middle[i],"_", thrshold_table$thr_upper[i],"_selected_features.RData",sep="")
  
  load(file)

  features_class1_top10 = c(features_class1_top10,selected_feature_imp1[1:10])
  features_class2_top10 = c(features_class2_top10,selected_feature_imp2[1:10])
}

sort(table(features_class1_top10))
sort(table(features_class2_top10))

features_class1 = as.data.frame(matrix(0,nrow = 11,ncol = length(unique(features_class1_top10))))
colnames(features_class1) = sort(unique(features_class1_top10))

features_class2 = as.data.frame(matrix(0,nrow = 11,ncol = length(unique(features_class2_top10))))
colnames(features_class2) = sort(unique(features_class2_top10))

for(i in 1:nrow(thrshold_table)){
  #file = paste("/project/lncrna/Xist/data_lisa/classification/silenced_vs_not_silenced_narrow/features_promoter_thr_",thrshold_table$thr_lower[i],"_",thrshold_table$thr_middle[i],"_", thrshold_table$thr_upper[i],"_selected_features.RData",sep="")
  file = paste("/project/lncrna/Xist/data_lisa/classification/early_vs_late_narrow/features_promoter_thr_",thrshold_table$thr_lower[i],"_",thrshold_table$thr_middle[i],"_", thrshold_table$thr_upper[i],"_selected_features.RData",sep="")
  load(file)
  
  class1 = selected_feature_imp1[1:10]
  features_class1[i,colnames(features_class1) %in% class1] = 1
  
  class2 = selected_feature_imp2[1:10]
  features_class2[i,colnames(features_class2) %in% class2] = 1 
}

pdf("/project/lncrna/Xist/plots_lisa/additional_analysis/feature_consistency_early_vs_late.pdf",width = 15,height = 10)
par(mfrow=c(1,1),mar=c(1,1,1,1),oma=c(3,3,3,10))
heatmap.2(as.matrix(t(features_class1)),Colv=NA,Rowv=NA,scale='none',cexCol=1.8,main='top 10 features early silenced',trace='none',dendrogram='none',col=rev(bluered(2)))
heatmap.2(as.matrix(t(features_class2)),Colv=NA,Rowv=NA,scale='none',cexCol=1.8,main='top 10 features late silenced',trace='none',dendrogram='none',col=rev(bluered(2)))
dev.off()
