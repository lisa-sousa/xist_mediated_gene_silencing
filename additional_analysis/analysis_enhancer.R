library(reshape2)

dir_data = "/Users/lisa/Desktop/data/"
dir_enhancers = paste(dir_data,"annotation_files/enhancers/",sep="")
file_HiCap_enhancers = "Promoter_Enhancer_Interactions.txt"
file_output = "gene_enhancers.bed"

hicap_enhancers = read.table(paste(dir_enhancers,file_HiCap_enhancers,sep=""),header = T, sep = "\t")
hicap_enhancers_chrX = hicap_enhancers[hicap_enhancers$Promoter.chr == "chrX" & hicap_enhancers$Fragment.chromosome == "chrX",]
enhancer_bed = data.frame(chr = hicap_enhancers_chrX$Fragment.chromosome, start = hicap_enhancers_chrX$Fragment.start.coordinate, end = hicap_enhancers_chrX$Fragment.end.coordinate,
                          Genes = hicap_enhancers_chrX$Gene.Name, score = rowMeans(cbind(hicap_enhancers_chrX$Read.pairs.in.HiCap.replicate.1,hicap_enhancers_chrX$Read.pairs.in.HiCap.replicate.2)),
                          strand = hicap_enhancers_chrX$Promoter.Strand)

enhancer_bed$enhancer_length = abs(enhancer_bed$start-enhancer_bed$end)
enhancer_bed$enhancer_center = rowMeans(enhancer_bed[,2:3])
enhancer_bed$start[enhancer_bed$enhancer_length < 1000] = enhancer_bed$enhancer_center[enhancer_bed$enhancer_length < 1000] - 500
enhancer_bed$end[enhancer_bed$enhancer_length < 1000] = enhancer_bed$enhancer_center[enhancer_bed$enhancer_length < 1000] + 500
enhancer_bed$enhancer_length = abs(enhancer_bed$start-enhancer_bed$end)

write.table(enhancer_bed[1:6],paste(dir_enhancers,file_output,sep=""),col.names = F, row.names = F, sep="\t",quote = F)


#######
#-> map chip-seq data
#-> create data matrix
#######

##############
#function
##############

get_wilcox_p_value <- function(column){
  feature0 = column[column$feature==0,] # the feature
  feature1 = column[column$feature==1,] # the feature
  wilcox = wilcox.test(feature0$halftime,feature1$halftime)$p.value
  return(wilcox)
}

get_wilcox_and_personcor <- function(column) {
  class0 = column[column$target==0,]
  class1 = column[column$target==1,]
  wilcox = wilcox.test(class0$feature,class1$feature)$p.value
  
  pears_corr = cor(column$feature,column$halftime)
  pears_pValue = cor.test(column$feature,column$halftime)$p.value
  
  return(list(wilcox,pears_corr,pears_pValue))
}

plot_binary_feature <- function(feature,column_all,column_strongest,column_close,n){
  wilcox_all = get_wilcox_p_value(column_all)
  #wilcox_all = p.adjust(get_wilcox_p_value(column_all), method = "BH", n=n)
  wilcox_strongest = get_wilcox_p_value(column_strongest)
  #wilcox_strongest = p.adjust(get_wilcox_p_value(wilcox_strongest), method = "BH", n=n)
  
  wilcox_close = get_wilcox_p_value(column_close)
  #wilcox_close = p.adjust(get_wilcox_p_value(wilcox_close), method = "BH", n=n)
  
  data_set_plot = rbind(column_all,column_strongest,column_close)[,c(4,1,3)]
  colnames(data_set_plot) = c("enhancer_set","overlap","value")
  
  par(mar=c(15,7,7,7))
  title_box_plot = paste(feature,
                         "\nwilcox test p-value for all enhancer: ",signif(wilcox_all,3),
                         "\nwilcox test p-value for strongest enhancer:",signif(wilcox_strongest,3),
                         "\nwilcox test p-value for closest enhancer:",signif(wilcox_close,3))
  
  ggbox = ggplot(data_set_plot, aes(x = enhancer_set, y = value, fill = overlap)) + geom_boxplot(alpha=0.7,notch = F) +
    scale_y_continuous(name = "half-time[days]") + scale_x_discrete(name = "enhancer set") + ggtitle(title_box_plot) +
    theme_bw() + theme(plot.title = element_text(size = 16, face = "bold",hjust = 0.5),
                       axis.title = element_text(size=15, face="bold"),
                       axis.text.x = element_text(size=15, angle = 60, hjust = 1),
                       axis.text.y = element_text(size=15),
                       axis.line = element_line(colour = "black"),
                       plot.margin = unit(c(2,2,2,2), "cm")) + 
    scale_fill_brewer(palette = "Greys")
  print(ggbox)
}

plot_continuous_feature <- function(feature,column_all,column_strongest,column_close,n){
  statistics_all = get_wilcox_and_personcor(column_all)
  wilcox_all = statistics_all[[1]]
  #wilcox_all = p.adjust(statistics_all[[1]], method = "BH", n=n)
  
  statistics_strongest = get_wilcox_and_personcor(column_strongest)    
  wilcox_strongest = statistics_strongest[[1]]
  #wilcox_strongest = p.adjust(statistics_strongest[[1]], method = "BH", n=n)
  
  statistics_close = get_wilcox_and_personcor(column_close)    
  wilcox_close = statistics_close[[1]]
  #wilcox_strongest = p.adjust(statistics_close[[1]], method = "BH", n=n)
  
  data_set_plot = rbind(column_all,column_strongest,column_close)[,c(4,2,1)]
  colnames(data_set_plot) = c("enhancer_set","silencing_class","value")
  data_set_plot$silencing_class = as.factor(data_set_plot$silencing_class)
  
  par(mar=c(15,7,7,7))
  title_box_plot = paste(feature,
                         "\nwilcox test p-value for all enhancer: ",signif(wilcox_all,3),
                         "\nwilcox test p-value for strongest enhancer:",signif(wilcox_strongest,3),
                         "\nwilcox test p-value for closest enhancer:",signif(wilcox_close,3))
  
  ggbox = ggplot(data_set_plot, aes(x = enhancer_set, y = value, fill = silencing_class)) + geom_boxplot(alpha=0.7,notch = F) +
    scale_y_continuous(name = "half-time[days]") + scale_x_discrete(name = "enhancer set") + ggtitle(title_box_plot) +
    theme_bw() + theme(plot.title = element_text(size = 16, face = "bold",hjust = 0.5),
                       axis.title = element_text(size=15, face="bold"),
                       axis.text.x = element_text(size=15, angle = 60, hjust = 1),
                       axis.text.y = element_text(size=15),
                       axis.line = element_line(colour = "black"),
                       plot.margin = unit(c(2,2,2,2), "cm")) + 
    scale_fill_brewer(palette = "Greys")
  print(ggbox)
}

###############
#enhancer sets
###############

source("/project/lncrna/Xist/xist_mediated_gene_silencing/modelling/model_functions.R")
output_dir = "/Users/lisa/work_stuff/projects/xist_epigenetics/computing/plots/additional_analysis/"
file_halftimes = paste(dir_data,"silencing_halftimes/fitted_data/halftimes_pro_seq_mm10_RSS_initial_ratio.txt",sep="")
file_feature_matrix = paste(dir_data,"modelling/feature_matrix/promoter_matrix_normRAdjusted_enhancer.RData",sep="")

hic_inteactions = c("mean_interaction_strength_HiC_all","mean_interaction_strength_HiC_promoter","mean_interaction_strength_HiC_xist")

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
n = ncol(data_set_all)-2

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
  
#choose closest enhancer per gene
data_set_close = data_set
genes = unique(as.character(data_set_close$Genes[duplicated(data_set_close$Genes)]))
data_set_close$TSS = 0
data_set_close$TSS[data_set_close$Strand == "+"] = data_set_close$Start[data_set_close$Strand == "+"]
data_set_close$TSS[data_set_close$Strand == "-"] = data_set_close$End[data_set_close$Strand == "-"]
for(i in 1:length(genes)){
  entry = data_set_close[data_set_close$Genes==genes[i],]
  enhancer_start = entry$start[which.min(abs(entry$TSS-entry$enhancer_center))]
  data_set_close = data_set_close[!(data_set_close$Genes==genes[i] & data_set_close$start != enhancer_start),]
}
data_set_close = data_set_close[!duplicated(data_set_close),]  
data_set_close = cbind(data_set_close[,10:86],target=data_set_close$target,halftime=data_set_close$halftime)


##plots boxplots
pdf(paste(output_dir,"enhancer_boxplots.pdf"),height = 10, width = 10)

for(i in 1:(ncol(data_set_all)-2)){
  
  feature = colnames(data_set_all)[i]
  column_all = cbind.data.frame(feature=data_set_all[,i],target=data_set_all$target, halftime=data_set_all$halftime, enhancer_set = factor("all"))
  column_strongest = cbind.data.frame(feature=data_set_strong[,i],target=data_set_strong$target, halftime=data_set_strong$halftime, enhancer_set = factor("strongest"))
  column_close = cbind.data.frame(feature=data_set_close[,i],target=data_set_close$target, halftime=data_set_close$halftime, enhancer_set = factor("closest"))
  
  if(is.factor(column_all$feature)){
    plot_binary_feature(feature,column_all,column_strongest,column_close,n)
  }else{
    plot_continuous_feature(feature,column_all,column_strongest,column_close,n)
    
  }
  if(feature %in% hic_inteactions){
    column_all$feature[column_all$feature > 0] = 1
    column_all$feature = as.factor(column_all$feature)
    
    column_strongest$feature[column_strongest$feature > 0] = 1
    column_strongest$feature = as.factor(column_strongest$feature)
    
    column_close$feature[column_close$feature > 0] = 1
    column_close$feature = as.factor(column_close$feature)
    
    plot_binary_feature(feature,column_all,column_strongest,column_close,n)
  }
}

dev.off()
