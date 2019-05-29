library(here)
library(ggplot2)

###################################################################################
#load data
###################################################################################

predictions = read.table(here('data/modelling/model/xci_escape_model/results_thr_0.9_-_1.6/','random_forest_predictions_new_genes.txt'),header = T)
mrna_seq = read.table(here('data/silencing_halftimes/fitted_data','halftimes_mrna_seq_undiff_mm10.bed'))
colnames(mrna_seq) = c("chr","start","end","Genes","halftime","strand")

table = merge(predictions,mrna_seq,by="Genes")
table = table[!(table$class > 0 & table$class < 1),] #remove instable predictions
table$class = as.factor(table$class)

###################################################################################
#plotting function
###################################################################################

###scatterplot of normalized AER vs gene predictions
scatterplot_dense_colors <- function(x1, x2, xlab, ylab, main){
  
  df = data.frame(x1,x2)
  
  ## Use densCols() output to get density at each point
  x = densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens = as.factor(col2rgb(x)[1,] + 1L)
  
  ## Map densities to colors
  cols = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  cols = colorRampPalette(c("grey","black"))(256)
  df$col = cols[df$dens]
  
  cor = cor.test(x1,x2)
  
  ## Plot it, reordering rows so that densest points are plotted on top
  ggplot = ggplot(data=df[order(df$dens),]) +
    geom_point(aes(x1,x2,color=dens),size=0.5) +
    scale_color_grey(start=0.7,end=0) + 
    labs(title=main, subtitle = paste('r=',round(cor$estimate,2),"\np=",signif(cor$p.value,2),sep="")) +
    geom_smooth(aes(x1,x2),method = "lm", se = FALSE,color="#ff8080", size=0.5) +
    scale_x_continuous(name=xlab, limits = c(0,1)) +
    scale_y_continuous(name=ylab) +
    theme_minimal(base_family = "Source Sans Pro") + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8), axis.text.y = element_text(size=8), 
          axis.title=element_text(size=8), legend.position = "none",plot.title = element_text(size=8), plot.subtitle = element_text(size=7)) 
  print(ggplot)
}

###################################################################################
#generate plots
###################################################################################

cairo_pdf(here('plots/additional_analysis','paper_figures_mRNA_seq_vs_predictions.pdf'),height = 3,width = 2)

wilcox = signif(wilcox.test(table$halftime[table$class == 0],table$halftime[table$class==1])$p.value,2)

ggplot(data=table,aes(x=class, y=halftime)) + 
  geom_jitter(position=position_jitter(0.2), size=0.6) + 
  stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean,geom="crossbar", color='grey', width=0.7, size=0.2)+
  geom_text(aes(x=1.5,y=4.5), label=paste('p=',wilcox,sep=""), size=2) +
  geom_segment(aes(x=1,xend=2,y=4,yend=4)) +
  theme_minimal() + theme(text=element_text(size=8), panel.grid.minor = element_blank())+ scale_y_continuous(limits=c(0,5),breaks=c(0,2,4,6), name='Halftime [days]') + 
  scale_x_discrete (name='Predicted Class')

dev.off()

cairo_pdf(here('plots/additional_analysis','plots_mrna_seq_vs_predictions.pdf'),height = 2,width = 2)

vote = 1-table$vote
cor = cor.test(vote, table$halftime)
main = paste("Correlation:",round(cor$estimate,2),"\np-value:",signif(cor$p.value,2))
scatterplot_dense_colors(vote,table$halftime,"vote (probability class 1)","halftime",main)

dev.off()
