
dir = "/project/lncrna/Xist/data/annotation_files/enhancers/"
file_HiCap_enhancers = "Promoter_Enhancer_Interactions.txt"
file_output = "gene_enhancers.bed"

hicap_enhancers = read.table(paste(dir,file_HiCap_enhancers,sep=""),header = T, sep = "\t")
hicap_enhancers_chrX = hicap_enhancers[hicap_enhancers$Promoter.chr == "chrX" & hicap_enhancers$Fragment.chromosome == "chrX",]
enhancer_bed = data.frame(chr = hicap_enhancers_chrX$Fragment.chromosome, start = hicap_enhancers_chrX$Fragment.start.coordinate, end = hicap_enhancers_chrX$Fragment.end.coordinate,
                          Genes = hicap_enhancers_chrX$Gene.Name, score = rowMeans(cbind(hicap_enhancers_chrX$Read.pairs.in.HiCap.replicate.1,hicap_enhancers_chrX$Read.pairs.in.HiCap.replicate.2)),
                          strand = hicap_enhancers_chrX$Promoter.Strand)

enhancer_bed$enhancer_length = abs(enhancer_bed$start-enhancer_bed$end)
enhancer_bed$enhancer_center = rowMeans(enhancer_bed[,2:3])
enhancer_bed$start[enhancer_bed$enhancer_length < 1000] = enhancer_bed$enhancer_center[enhancer_bed$enhancer_length < 1000] - 500
enhancer_bed$end[enhancer_bed$enhancer_length < 1000] = enhancer_bed$enhancer_center[enhancer_bed$enhancer_length < 1000] + 500
enhancer_bed$enhancer_length = abs(enhancer_bed$start-enhancer_bed$end)

write.table(enhancer_bed[1:6],paste(dir,file_output,sep=""),col.names = F, row.names = F, sep="\t",quote = F)


#######
#-> map chip-seq data
#-> create data matrix
#######

##############
#function
##############

plot_data_boxpots <- function(output_directory_plots_thr, data_set_all, data_set_strong){
  
  CairoPDF(file = paste(output_directory_plots_thr,'data_boxplots.pdf',sep=''), width = 15, height = 15)
  hic_inteactions = c("mean_interaction_strength_HiC_all","mean_interaction_strength_HiC_promoter","mean_interaction_strength_HiC_xist")
  
  data_set_all$enhancer_set = "all_enhancer"
  data_set_strong$enhancer_set = "strongest_enhancer"
  data_set_plot = rbind(data_set_all,data_set_strong)
  plot_df = melt(data_set_plot,id.vars = c("enhancer_set"))
  
  for(i in 1:ncol(data_set_all)){
    
    feature = colnames(data_set_all)[i]
    column_all = cbind.data.frame(feature=data_set_all[,i],target=data_set_all$target, halftime=data_set_all$halftime)
    
    
    
    
    if(is.factor(column)){
      
      column_target = cbind.data.frame(column,halftime)
      feature0 = column_target[column_target[,1]==0,] # the feature
      feature1 = column_target[column_target[,1]==1,] # the feature
      
      wilcox = wilcox.test(feature0[,2],feature1[,2])$p.value
      
      par(mar=c(15,7,7,7))
      title_box_plot = paste(feature,"\nwilcox p-value:",signif(wilcox,4))
      gg_box = ggplot(column_target, aes(x=column,y=halftime)) + geom_boxplot(notch=FALSE,fill = 'lightgrey', colour = 'black',alpha = 0.7,outlier.shape = 16) + ggtitle(title_box_plot) + theme_bw() + 
        theme(axis.text.x=element_text(hjust = 1,size=20),axis.text.y=element_text(size=20),axis.title = element_text(face="bold", size=20),
              plot.title = element_text(hjust = 0.5,size=35,face='bold'),plot.margin = unit(c(3,3,3,3), "cm"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
        scale_x_discrete(name = "overlap") + scale_y_continuous(name = "halftime")
      print(gg_box)
    }else{
      
      column_target = cbind.data.frame(column,halftime,target)
      class0 = column_target[column_target$target==0,]
      class1 = column_target[column_target$target==1,]
      
      pears_corr = cor(column_target[,1],column_target[,2])
      pears_pValue = cor.test(column_target[,1],column_target[,2])$p.value
      
      wilcox = wilcox.test(class0[,1],class1[,1])$p.value
      
      par(mar=c(15,7,7,7))
      title_box_plot = paste(feature,"\npearson cor",signif(pears_corr,4),"with p-value:",signif(pears_pValue,4),"\n wilcox test:",signif(wilcox,4))
      gg_box = ggplot(column_target, aes(x=target,y=column)) + geom_boxplot(notch=FALSE,fill = 'lightgrey', colour = 'black',alpha = 0.7,outlier.shape = 16) + ggtitle(title_box_plot) + theme_bw() + 
        theme(axis.text.x=element_text(hjust = 1,size=20),axis.text.y=element_text(size=20),axis.title = element_text(face="bold", size=20),
              plot.title = element_text(hjust = 0.5,size=35,face='bold'),plot.margin = unit(c(3,3,3,3), "cm"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
        scale_x_discrete(name = "class") + scale_y_continuous(name = "feature") 
      print(gg_box)
    }
    if(feature %in% hic_inteactions){
      
      column_target = cbind.data.frame(column,halftime)
      column_target$column[column_target$column > 0] = 1
      column_target$column = as.factor(column_target$column)
      
      feature0 = column_target[column_target[,1]==0,] # the feature
      feature1 = column_target[column_target[,1]==1,] # the feature
      
      wilcox = wilcox.test(feature0[,2],feature1[,2])$p.value
      
      par(mar=c(15,7,7,7))
      title_box_plot = paste(feature,"\nwilcox p-value:",signif(wilcox,4))
      gg_box = ggplot(column_target, aes(x=column,y=halftime)) + geom_boxplot(notch=FALSE,fill = 'lightgrey', colour = 'black',alpha = 0.7,outlier.shape = 16) + ggtitle(title_box_plot) + theme_bw() + 
        theme(axis.text.x=element_text(hjust = 1,size=20),axis.text.y=element_text(size=20),axis.title = element_text(face="bold", size=20),
              plot.title = element_text(hjust = 0.5,size=35,face='bold'),plot.margin = unit(c(3,3,3,3), "cm"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
        scale_x_discrete(name = "overlap") + scale_y_continuous(name = "halftime")
      print(gg_box)
    }
  }
  dev.off()
}


plot_error_rates <- function(output_directory_plots_thr, random_forest_model_all_features, random_forest_model_top_features){
  
  CairoPDF(file = paste(output_directory_plots_thr,'error_rate.pdf',sep=''), width = 10, height = 10)
  
  
  
  par(mfrow=c(1,1),mar=c(15,10,10,5),oma=c(5,5,5,5))
  
  ggbox = ggplot(plot_df, aes(x = variable, y = value, fill = enhancer_set)) + geom_boxplot(alpha=0.7,notch = T) +
    scale_y_continuous(name = "error rate (%)", limits=c(0, 50)) + scale_x_discrete(name = "") + ggtitle("Error rate for Random Forest models") +
    theme_bw() + theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold",hjust = 0.5),text = element_text(size = 12, family = "Tahoma"),
                       axis.title = element_text(face="bold"),axis.text.x=element_text(angle = 60, hjust = 1,size=11),plot.margin = unit(c(2,2,2,2), "cm"),
                       axis.line = element_line(colour = "black")) + 
    scale_fill_brewer(palette = "Greys")
  
  print(ggbox)
  
  dev.off()
}

source("/project/lncrna/Xist/xist_mediated_gene_silencing/modelling/model_functions.R")
output_dir = "/project/lncrna/Xist/plots/additional_analysis/"
file_halftimes = "/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_pro_seq_mm10_RSS_initial_ratio.txt"
file_feature_matrix = "/project/lncrna/Xist/data/modelling/feature_matrix/promoter_matrix_normRAdjusted_enhancer.RData"


thr_silencing_lower = 0.9
thr_silencing_middle = "-"
thr_silencing_upper = 2


halftimes = read.table(file_halftimes,header=T)

load(file_feature_matrix)
colnames(data_set)[1] = "Genes"
data_set$enhancer_strength = halftime
data_set = cbind(data_set,enhancer_bed[,c(2,3,8)])
data_set = merge(halftimes,data_set,by="Genes")
data_set = data_set[data_set$halftime < thr_silencing_lower | data_set$halftime > thr_silencing_upper,]
data_set$target = 0
data_set$target[data_set$halftime > 2] = 1

data_set_all = cbind(data_set[,10:86],target=data_set$target,halftime=data_set$halftime)

#all enhancers per gene
data_all=list()
data_all[[1]] = data_set[,10:86]
data_all[[2]] = as.factor(data_set$target)
data_all[[3]] = data_set$halftime

plot_data_boxpots(output_dir, data)


#choose strongest enhancer
data_set_strong = data_set
genes = unique(as.character(data_set_strong$Genes[duplicated(data_set_strong$Genes)]))
for(i in 1:length(genes)){
  entry = data_set_strong[data_set_strong$Genes==genes[i],]
  max = max(entry$enhancer_strength)
  data_set_strong = data_set_strong[!(data_set_strong$Genes==genes[i] & data_set_strong$enhancer_strength != max),]
}
data_set_strong = data_set_strong[!duplicated(data_set_strong),]  
data_set_strong = cbind(data_set_strong[,10:86],target=data_set_strong$target,halftime=data_set_strong$halftime)
  
data_strong=list()
data_strong[[1]] = data_set_strong[,10:86]
data_strong[[2]] = as.factor(data_set_strong$target)
data_strong[[3]] = data_set_strong$halftime

plot_data_boxpots(output_dir, data)
  
#choose closest enhancer per gene
genes = unique(as.character(data_set$Genes[duplicated(data_set$Genes)]))
data_set$TSS = 0
data_set$TSS[data_set$Strand == "+"] = data_set$Start[data_set$Strand == "+"]
data_set$TSS[data_set$Strand == "-"] = data_set$End[data_set$Strand == "-"]
for(i in 1:length(genes)){
  entry = data_set[data_set$Genes==genes[i],]
  enhancer_start = entry$start[which.min(abs(entry$TSS-entry$enhancer_center))]
  data_set = data_set[!(data_set$Genes==genes[i] & data_set$start != enhancer_start),]
}
data_set = data_set[!duplicated(data_set),]  



#choose only active enhancer
