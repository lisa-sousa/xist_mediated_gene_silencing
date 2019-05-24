library(here)
library(ggplot2)
library(Cairo)

###################################################################################
#directories and files
###################################################################################

file_mutant_fc = here("data/annotation_files/xist_mutants","FC_Dox_noDox_table_all_mutants_norm.txt")
file_cluster_matrix = here("data/modelling/model/silencing_dynamics_model/results_thr_0.5_0.9_1.3","best_clustering_data_set_k3.RData")

###################################################################################
#load mutant data
###################################################################################

table_mutants = read.table(file=file_mutant_fc,header=T,sep="\t")
mutants = c("WT","A","BC")

###################################################################################
#load cluster data
###################################################################################

load(file_cluster_matrix)
table_cluster = data_set_plot
table_cluster$gene = row.names(table_cluster)

table = merge(table_mutants,table_cluster,by="gene")

###################################################################################
#plot foldchange vs cluster
###################################################################################

CairoPDF(file = here("plots/additional_analysis","analysis_paper_bousard.pdf"), width = 4, height = 4)

for(i in 1:3){
  mutant = mutants[i]
  data_plot = data.frame(cluster=table$cluster,mutant=table[,colnames(table)==mutant])
  gg_box = ggplot(data_plot, aes(x=cluster,y=mutant)) + 
    geom_boxplot(colour = "#4d4d4d",alpha = 0.7,outlier.size=0.1,lwd=0.4) + 
    ggtitle(paste("mutant:",mutant)) + 
    scale_x_discrete(name = "cluster") + 
    scale_y_continuous(name = "log2 foldchange",limits = c(-4, 1)) +
    theme_minimal(base_family = "Source Sans Pro") + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=7), axis.text.y = element_text(size=8), 
          axis.title=element_text(size=8, margin = margin(t=0)),plot.title = element_text(size=7)) 
  print(gg_box)
  anova = aov(mutant ~ cluster, data = data_plot)
  print(TukeyHSD(anova))
}

dev.off()

###################################################################################
#analysis mutant repeat A
###################################################################################

###get change in silencing between WT and mutant
###mutanta A
thr_mutant_A = 1
table$change_mutant_A = table$A - table$WT
table$reduced_silencing_mutant_A = 0
table$reduced_silencing_mutant_A[table$change_mutant_A > thr_mutant_A] = 1

mutant_A_silenced = as.numeric(table$cluster[table$reduced_silencing_mutant_A == 0])
mutant_A_impaired_silenced = as.numeric(table$cluster[table$reduced_silencing_mutant_A == 1])

total_n_of_genes_in_each_cluster = table(table_cluster$cluster)
number_of_clusters = as.numeric(names(rev(total_n_of_genes_in_each_cluster))[1])

table_plot = rbind(table(factor(mutant_A_impaired_silenced,levels = 1:number_of_clusters))/total_n_of_genes_in_each_cluster,
                   table(factor(mutant_A_silenced,levels = 1:number_of_clusters))/total_n_of_genes_in_each_cluster)
table_plot = rbind(1 - colSums(table_plot),table_plot)
rownames(table_plot) = c("not covered", "repeat A dependent genes", "repeat A indenpendent genes")

cairo_pdf(here('plots/additional_analysis','paper_figures_cluster_repeatA_bousard.pdf'),width = 2,height = 3, onefile = TRUE)
ggplot(aes(x=Var2,y=value,fill=Var1),data=melt(table_plot*100)) +
  geom_bar(stat="identity",color="black") +
  scale_x_continuous(breaks=c(1,2,3), label=c("1", "2", "3"), name='cluster') +
  scale_y_continuous(breaks=c(0,25,50,75,100), name='fraction of genes (%)') +
  scale_fill_grey(name="group",start=0.1,end=0.9) +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8), axis.text.y = element_text(size=8), 
        axis.title=element_text(size=8),legend.text = element_text(size=8), legend.title = element_text(size=8), legend.position = "bottom") +
  guides(fill=guide_legend(nrow=3), col=guide_legend(nrow=3))
dev.off()

####Fisher test
cluster1 = 1
cluster2 = 2
#repeat A dependent
dependent = table(mutant_A_silenced)
#repeat independent + not covered
independent = total_n_of_genes_in_each_cluster - dependent
#fisher test
matrix = matrix(c(dependent[cluster1],dependent[cluster2],independent[cluster1],independent[cluster2]), nr=2)
fisher.test(matrix)

###################################################################################
#analysis mutant repeat BC
###################################################################################

###get change in silencing between WT and mutant
###mutanta B
thr_mutant_B = 0.5
table$change_mutant_BC = table$BC - table$WT
table$reduced_silencing_mutant_BC = 0
table$reduced_silencing_mutant_BC[table$change_mutant_BC > thr_mutant_B] = 1

mutant_BC_silenced = as.numeric(table$cluster[table$reduced_silencing_mutant_BC == 0])
mutant_BC_impaired_silenced = as.numeric(table$cluster[table$reduced_silencing_mutant_BC == 1])

total_n_of_genes_in_each_cluster = table(table_cluster$cluster)
number_of_clusters = as.numeric(names(rev(total_n_of_genes_in_each_cluster))[1])

table_plot = rbind(table(factor(mutant_BC_impaired_silenced,levels = 1:number_of_clusters))/total_n_of_genes_in_each_cluster,
                   table(factor(mutant_BC_silenced,levels = 1:number_of_clusters))/total_n_of_genes_in_each_cluster)
table_plot = rbind(1 - colSums(table_plot),table_plot)
rownames(table_plot) = c("not covered", "repeat BC dependent genes", "repeat BC indenpendent genes")


cairo_pdf(here('plots/additional_analysis','paper_figures_cluster_repeatBC_bousard.pdf'),width = 2,height = 3, onefile = TRUE)
ggplot(aes(x=Var2,y=value,fill=Var1),data=melt(table_plot*100)) +
  geom_bar(stat="identity",color="black") +
  scale_x_continuous(breaks=c(1,2,3), label=c("1", "2", "3"), name='cluster') +
  scale_y_continuous(breaks=c(0,25,50,75,100), name='fraction of genes (%)') +
  scale_fill_grey(name="group",start=0.1,end=0.9) +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8), axis.text.y = element_text(size=8), 
        axis.title=element_text(size=8),legend.text = element_text(size=8), legend.title = element_text(size=8), legend.position = "bottom") +
  guides(fill=guide_legend(nrow=3), col=guide_legend(nrow=3))
dev.off()

####Fisher test
cluster1 = 2
cluster2 = 1
#repeat BC dependent
dependent = table(mutant_BC_silenced)
#repeat independent + not covered
independent = total_n_of_genes_in_each_cluster - dependent
#fisher test
matrix = matrix(c(dependent[cluster1],dependent[cluster2],independent[cluster1],independent[cluster2]), nr=2)
fisher.test(matrix)

