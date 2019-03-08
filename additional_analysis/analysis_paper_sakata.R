###################################################################################
#libraries
###################################################################################

library(Cairo)
library(ggplot2)

###################################################################################
#directories and files
###################################################################################

CAGdelta5_file = "/project/lncrna/Xist/data/annotation_files/CAGdelta5/GSE93031_allelic_expression.txt"
model_matrix_file = "/project/lncrna/Xist/data/modelling/feature_matrix/promoter_matrix_reannotated_normRAdjusted_pro_seq_genes.RData"
clustering_matrix = "/project/lncrna/Xist/data/modelling/model/silencing_dynamics_model/results_thr_0.5_0.9_1.3/best_clustering_data_set_k3.RData"
output_dir = "/project/lncrna/Xist/plots/additional_analysis/" 

###################################################################################
#load CAGdelta5 data
###################################################################################

table = read.table(file = CAGdelta5_file,sep='\t',header=T)
table_CAGdelta5 = table[,c(1,2,7)]
table_CAGdelta5_chrX = table_CAGdelta5[table_CAGdelta5$chrom=="chrX",]
table_CAGdelta5_chrX = na.omit(table_CAGdelta5_chrX)

#threhold for silenced vs not silenced as defined in Sakata et. al.
thr_silencing = 10

genes_CAGdetla5_silenced = table_CAGdelta5_chrX$gene[table_CAGdelta5_chrX$avePat.CAGdelta5 <= thr_silencing]
genes_CAGdetla5_silenced = as.character(genes_CAGdetla5_silenced)

genes_CAGdetla5_not_silenced = table_CAGdelta5_chrX$gene[table_CAGdelta5_chrX$avePat.CAGdelta5 > thr_silencing]
genes_CAGdetla5_not_silenced = as.character(genes_CAGdetla5_not_silenced)

###################################################################################
#load feature matrix
###################################################################################

load(model_matrix_file)
table_halftimes = data.frame(gene = rownames(data_set), halftime = halftime)

table_halftimes_silenced = table_halftimes[table_halftimes$gene %in% genes_CAGdetla5_silenced,]
table_halftimes_not_silenced = table_halftimes[table_halftimes$gene %in% genes_CAGdetla5_not_silenced,]

###################################################################################
#generate plots
###################################################################################

####plot distribution and kumulative distribution of halftime for silenced vs not silenced CAGdelta5 genes
plot_file = paste(output_dir,"analysis_cag_delta_5.pdf",sep="")
CairoPDF(file = plot_file, width = 10, height = 10)
par(mfrow=c(1,1),mar=c(5,5,5,5),oma=c(2,2,2,2))

plot(density(table_halftimes_silenced$halftime),main='Halftime distribution of CAgdelta5 genes',col='darkorange')
lines(density(table_halftimes_not_silenced$halftime),col='darkred')
legend('topright',legend = c('CAGdelta5 silenced genes','CAGdelta5 not silenced genes'), col = c('darkorange','darkred'),pch='_')

n = nrow(table_halftimes_silenced)
m = nrow(table_halftimes_not_silenced)
plot(sort(table_halftimes_silenced$halftime), (1:n)/n, type = 's', ylim = c(0, 1), xlab = 'halftime',col='darkorange', ylab = '', main = 'Empirical Cumluative Distribution\nCAGdelta5')
lines(sort(table_halftimes_not_silenced$halftime), (1:m)/m, type = 's',col='darkred')
legend('bottomright',legend = c('CAGdelta5 silenced genes','CAGdelta5 not silenced genes'), col = c('darkorange','darkred'),pch='_')


###barplot of into which clusters (silenced vs not silenced model) the CAGdelta5 silenced genes fall
load(clustering_matrix)
data_set_plot$cluster = as.numeric(data_set_plot$cluster)

###################################################################################
#calculate fration of repeat A dependent and independent genes per cluster
###################################################################################

#select genes that are silenced in CAGdelta5 mutant
matrix_CAGdelta5_silenced = data_set_plot[rownames(data_set_plot) %in% genes_CAGdetla5_silenced,]
cluster_silenced = matrix_CAGdelta5_silenced$cluster

matrix_CAGdelta5_not_silenced = data_set_plot[rownames(data_set_plot) %in% genes_CAGdetla5_not_silenced,]
cluster_not_silenced = matrix_CAGdelta5_not_silenced$cluster

total_n_of_genes_in_each_cluster = table(data_set_plot$cluster)

par(mfrow=c(1,1),mar=c(5,5,5,5),oma=c(2,2,2,2))

table = rbind(table(factor(cluster_silenced,levels = 1:max(data_set_plot$cluster)))/total_n_of_genes_in_each_cluster,table(factor(cluster_not_silenced,levels = 1:max(data_set_plot$cluster)))/total_n_of_genes_in_each_cluster)
table = rbind(table,1 - colSums(table))
rownames(table) = c("repA indenpendent genes", "repA dependent genes & escapees", "not covered")

barplot(table*100,xlab="cluster",col = c('black','lightgrey','white'), ylab="fraction of genes (%)",main = "CAGdelta5 genes distribution among cluster", ylim=c(0,100),
        legend = rownames(table))


###plot boxplots of silenced vs not silenced CAGdelta5 genes for each feature in our model
matrix_CAGdelta5_silenced = data_set[rownames(data_set) %in% genes_CAGdetla5_silenced,]
matrix_not_CAGdelta5_silenced = data_set[rownames(data_set) %in% genes_CAGdetla5_not_silenced,]

n = nrow(matrix_CAGdelta5_silenced)
m = nrow(matrix_not_CAGdelta5_silenced)

par(mfrow=c(1,1),mar=c(15,10,10,5),oma=c(5,5,5,5))

for(i in 1:ncol(matrix_CAGdelta5_silenced)){
  feature = colnames(matrix_CAGdelta5_silenced)[i]
  if(!is.factor(matrix_CAGdelta5_silenced[,i])){
    p_value = wilcox.test(matrix_CAGdelta5_silenced[,i],matrix_not_CAGdelta5_silenced[,i])$p.value
  }else{p_value = NA}
  print(feature)
  data_plot = data.frame(genes = c(rep("CAGdelta5_sil",n),rep("CAGdelta5_not_sil",m)), feature = c(matrix_CAGdelta5_silenced[,i],matrix_not_CAGdelta5_silenced[,i]))
  
  gg_box = ggplot(data_plot, aes(x=genes,y=feature)) + geom_boxplot(notch=FALSE,fill = "lightgrey", colour = "black",alpha = 0.7,outlier.shape = 16) + 
    ggtitle(paste(feature,"P<",round(p_value,4))) + theme_bw() + 
    theme(axis.text.x=element_text(angle = 90,hjust = 1,size=20),axis.text.y=element_text(size=20),axis.title = element_text(face="bold", size=20),
          plot.title = element_text(hjust = 0.5,size=35,face='bold'),plot.margin = unit(c(3,3,3,3), "cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
    scale_x_discrete(name = "genes") + scale_y_continuous(name = "signal")
  print(gg_box)
}

dev.off()

####Fisher test
cluster1 = 1
cluster2 = 3
#repeat A dependent
dependent = table(cluster_not_silenced)
#repeat independent + not covered
independent = total_n_of_genes_in_each_cluster - dependent
#fisher test
matrix = matrix(c(dependent[cluster1],dependent[cluster2],independent[cluster1],independent[cluster2]), nr=2)
fisher.test(matrix)


