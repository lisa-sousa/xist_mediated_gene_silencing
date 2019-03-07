##################################
#read data
##################################

predictions = read.table("/project/lncrna/Xist/data/modelling/model/xci_escape_model/results_thr_0.9_-_1.6/random_forest_predictions_new_genes.txt",header = T)
mrna_seq = read.table("/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_mrna_seq_undiff_mm10.bed", header = F)
colnames(mrna_seq) = c("chr","start","end","Genes","halftime","strand")

table = merge(predictions,mrna_seq,by="Genes")
table = table[!(table$class > 0 & table$class < 1),] #remove instable predictions
table$class = as.factor(table$class)

##################################
#plotting function
##################################

scatterplot_dense_colors <- function(x1, x2, xlab, ylab, main){
  
  df = data.frame(x1,x2)
  
  ## Use densCols() output to get density at each point
  x = densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens = col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  cols = colorRampPalette(c("grey","black"))(256)
  df$col = cols[df$dens]
  
  ## Plot it, reordering rows so that densest points are plotted on top
  plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, xlim = c(0,1),ylim = c(0,3.5), ylab = ylab, xlab = xlab,main = main)
  abline(lm(x2~x1))
}

##################################
#generate plots
##################################

pdf("/project/lncrna/Xist/plots/additional_analysis/mRNA_seq_vs_predictions.pdf",height = 2,width = 2)

wilcox = signif(wilcox.test(table$halftime[table$class == 0],table$halftime[table$class==1])$p.value,2)

ggplot(data=table,aes(x=class, y=halftime)) + 
  geom_jitter(position=position_jitter(0.2), size=0.6) + 
  stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean,geom="crossbar", color='grey', width=0.7, size=0.2)+
  geom_text(aes(x=1.5,y=4.5), label=paste('p=',wilcox,sep=""), size=2) +
  geom_segment(aes(x=1,xend=2,y=4,yend=4)) +
  theme_minimal() + theme(text=element_text(size=8), panel.grid.minor = element_blank())+ scale_y_continuous(limits=c(0,5),breaks=c(0,2,4,6), name='Halftime [days]') + 
  scale_x_discrete (name='Predicted Class')
#boxplot(table$halftime ~ table$class, main=paste("Wilcox test:",wilcox))

vote = 1-table$vote
cor = cor.test(vote, table$halftime)
main = paste("Correlation:",round(cor$estimate,2),", p-value:",signif(cor$p.value,2))
scatterplot_dense_colors(vote,table$halftime,"vote (probability class 1)","halftime",main)

dev.off()


