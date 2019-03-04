####directories and files 
file_mutant_fc = "/project/lncrna/Xist/data/annotation_files/xist_mutants/FC_Dox_noDox_table_all_mutants_norm.txt"
file_cluster_matrix = "/project/lncrna/Xist/data/modelling/model/silencing_dynamics_model/results_thr_0.5_0.9_1.3/best_clustering_data_set_k3.RData"
output_dir = "/project/lncrna/Xist/plots/additional_analysis/" 

###load mutant data
table_mutants = read.table(file=file_mutant_fc,header=T,sep="\t")
mutants = c("WT","A","BC")

###load cluster data
load(file_cluster_matrix)
table_cluster = data_set_plot
table_cluster$gene = row.names(table_cluster)

table = merge(table_mutants,table_cluster,by="gene")

###plot foldchange vs cluster
CairoPDF(file = "/project/lncrna/Xist/plots/additional_analysis/analysis_paper_mutants.pdf", width = 15, height = 15)
par(mfrow=c(1,1),mar=c(15,10,10,5),oma=c(5,5,5,5))

for(i in 1:3){
  mutant = mutants[i]
  data_plot = data.frame(cluster=table$cluster,mutant=table[,colnames(table)==mutant])
  gg_box = ggplot(data_plot, aes(x=cluster,y=mutant)) + geom_boxplot(notch=FALSE,fill = "lightgrey", colour = "black",alpha = 0.7,outlier.shape = 16) + ggtitle(paste("mutant:",mutant)) + theme_bw() + 
    theme(axis.text.x=element_text(hjust = 0.5,size=20),axis.text.y=element_text(size=20),axis.title = element_text(face="bold", size=20),
          plot.title = element_text(hjust = 0.5,size=35,face='bold'),plot.margin = unit(c(3,3,3,3), "cm"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
    scale_x_discrete(name = "cluster") + scale_y_continuous(name = "log2 foldchange",limits = c(-4, 1))
  print(gg_box)
  anova = aov(mutant ~ cluster, data = data_plot)
  print(TukeyHSD(anova))
}

###get change in silencing between WT and mutant
###mutanta A
thr_mutant_A = 1
table$change_mutant_A = table$A - table$WT
table$reduced_silencing_mutant_A = 0
table$reduced_silencing_mutant_A[table$change_mutant_A > thr_mutant_A] = 1

plot(density(table$change_mutant_A),main="difference Mutant A - WT")
abline(v=thr_mutant_A)

mutant_A_silenced = as.numeric(table$cluster[table$reduced_silencing_mutant_A == 0])
mutant_A_impaired_silenced = as.numeric(table$cluster[table$reduced_silencing_mutant_A == 1])

total_n_of_genes_in_each_cluster = table(table_cluster$cluster)
number_of_clusters = as.numeric(names(rev(total_n_of_genes_in_each_cluster))[1])

par(mfrow=c(1,1),mar=c(5,5,5,5),oma=c(2,2,2,2))

table_plot = rbind(table(factor(mutant_A_silenced,levels = 1:number_of_clusters))/total_n_of_genes_in_each_cluster,table(factor(mutant_A_impaired_silenced,levels = 1:number_of_clusters))/total_n_of_genes_in_each_cluster)
table_plot = rbind(table_plot,1 - colSums(table_plot))
rownames(table_plot) = c("repeat indenpendent genes", "repeat dependent genes", "not covered")

barplot(table_plot*100,xlab="cluster",col = c('black','lightgrey','white'), ylab="fraction of genes (%)",main = "mutant A gene distribution among cluster", ylim=c(0,100),
        legend = rownames(table_plot))

####Fisher test
cluster1 = 3
cluster2 = 1
#repeat A dependent
dependent = table(mutant_A_silenced)
#repeat independent + not covered
independent = total_n_of_genes_in_each_cluster - dependent
#fisher test
matrix = matrix(c(dependent[cluster1],dependent[cluster2],independent[cluster1],independent[cluster2]), nr=2)
fisher.test(matrix)



###get change in silencing between WT and mutant
###mutanta B
thr_mutant_B = 0.5
table$change_mutant_BC = table$BC - table$WT
table$reduced_silencing_mutant_BC = 0
table$reduced_silencing_mutant_BC[table$change_mutant_BC > thr_mutant_B] = 1

plot(density(table$change_mutant_BC),main="difference Mutant BC - WT")
abline(v=thr_mutant_B)

mutant_BC_silenced = as.numeric(table$cluster[table$reduced_silencing_mutant_BC == 0])
mutant_BC_impaired_silenced = as.numeric(table$cluster[table$reduced_silencing_mutant_BC == 1])

total_n_of_genes_in_each_cluster = table(table_cluster$cluster)
number_of_clusters = as.numeric(names(rev(total_n_of_genes_in_each_cluster))[1])

par(mfrow=c(1,1),mar=c(5,5,5,5),oma=c(2,2,2,2))

table_plot = rbind(table(factor(mutant_BC_silenced,levels = 1:number_of_clusters))/total_n_of_genes_in_each_cluster,table(factor(mutant_BC_impaired_silenced,levels = 1:number_of_clusters))/total_n_of_genes_in_each_cluster)
table_plot = rbind(table_plot,1 - colSums(table_plot))
rownames(table_plot) = c("repeat indenpendent genes", "repeat dependent genes", "not covered")

barplot(table_plot*100,xlab="cluster",col = c('black','lightgrey','white'), ylab="fraction of genes (%)",main = "mutant BC gene distribution among cluster", ylim=c(0,100),
        legend = rownames(table_plot))


####Fisher test
cluster1 = 1
cluster2 = 3
#repeat BC dependent
dependent = table(mutant_BC_silenced)
#repeat independent + not covered
independent = total_n_of_genes_in_each_cluster - dependent
#fisher test
matrix = matrix(c(dependent[cluster1],dependent[cluster2],independent[cluster1],independent[cluster2]), nr=2)
fisher.test(matrix)

dev.off()

