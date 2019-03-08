###################################################################################
#load data
###################################################################################

pro_seq_file = "/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_pro_seq_mm10.bed"
pro_seq_table = read.table(pro_seq_file,sep='\t',header=F)
colnames(pro_seq_table) = c("Chromosomes","Start","End","Genes","halftime","Strand")

mrna_seq_undiff_file = "/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_mrna_seq_undiff_mm10.txt"
mrna_seq_undiff_table = read.table(mrna_seq_undiff_file,sep='\t',header=T)

mrna_seq_diff_file = "/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_mrna_seq_diff_mm10.txt"
mrna_seq_diff_table = read.table(mrna_seq_diff_file,sep='\t',header=T)

pro_mrna_undiff_table = merge(pro_seq_table,mrna_seq_undiff_table,by='Genes') 
pro_mrna_diff_table = merge(pro_seq_table,mrna_seq_diff_table,by='Genes') 
mrna_undiff_diff_table = merge(mrna_seq_undiff_table,mrna_seq_diff_table,by='Genes') 

###################################################################################
#plotting function
###################################################################################

scatterplot_dense_colors <- function(x1, x2, xlab, ylab){

  df = data.frame(x1,x2)
  
  ## Use densCols() output to get density at each point
  x = densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens = col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  cols = colorRampPalette(c("grey","black"))(256)
  df$col = cols[df$dens]
  
  ## Plot it, reordering rows so that densest points are plotted on top
  plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, xlim = c(0,3.5),ylim = c(0,3.5), ylab = ylab, xlab = xlab)
  abline(lm(x2~x1))
}

###################################################################################
#generate plots
###################################################################################

pdf(file="/project/lncrna/Xist/plots/additional_analysis/mrna_seq_vs_pro_seq.pdf")

scatterplot_dense_colors(pro_mrna_undiff_table$halftime.x,pro_mrna_undiff_table$halftime.y,'halftime pro-seq','halftime mRNA-Seq undiff')
scatterplot_dense_colors(pro_mrna_diff_table$halftime.x,pro_mrna_diff_table$halftime.y,'halftime pro-seq','halftime mRNA-Seq diff')
scatterplot_dense_colors(mrna_undiff_diff_table$halftime.x,mrna_undiff_diff_table$halftime.y,'halftime mRNA-Seq undiff','halftime mRNA-Seq diff')

dev.off()





