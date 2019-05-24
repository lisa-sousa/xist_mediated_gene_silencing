###################################################################################
#libraries
###################################################################################

library(ggplot2)
library(RColorBrewer)
library(here)

###################################################################################
#directories
###################################################################################

output_dir = "plots/additional_analysis"

input_dir_clones = "data/modelling/model/clones"
predictions_clones = c("predictions_clones_86.txt","predictions_clones_87.txt","predictions_clones_109.txt","predictions_clones_190.txt","predictions_clones_228.txt","predictions_clones_273.txt")

input_dir = "data/annotation_files/xist_transgenes"
clone = c("86","87","109","190","228","273")
chr = c("chrX","chrX","chrX","chrX","chr12","chr12")

###################################################################################
#parameters
###################################################################################

pseudocount = 0.001
thr_foldchange = c()
B=1000

###################################################################################
#functions to calculate normalized AER, do permutation test and plot clone distribution
###################################################################################

###read in clone data -> AER
get_clone_table <- function(file, chr){
  table_clone = read.table(file=file,header = T)
  table_clone = table_clone[grep(chr,table_clone$coord),]
  table_clone = table_clone[,c(1,6)]
  colnames(table_clone) = c("Genes","ratio")
  return(table_clone)
}

###calculate normalized AER
get_foldchange <- function(input_dir,chr,clone,pseudocount){
  table_clone_no_Dox_rep1 = get_clone_table(here(paste(input_dir,"/",chr,"_clones",sep=""),paste(clone,"_noDox_1.txt",sep="")),chr)
  table_clone_no_Dox_rep2 = get_clone_table(here(paste(input_dir,"/",chr,"_clones",sep=""),paste(clone,"_noDox_2.txt",sep="")),chr)
  table_clone_no_Dox = merge(table_clone_no_Dox_rep1,table_clone_no_Dox_rep2,by="Genes")
  table_clone_no_Dox$mean_ratio = rowMeans(table_clone_no_Dox[,2:3])
  
  table_clone_five_days_rep1 = get_clone_table(here(paste(input_dir,"/",chr,"_clones",sep=""),paste(clone,"_5Days_1.txt",sep="")),chr)
  table_clone_five_days_rep2 = get_clone_table(here(paste(input_dir,"/",chr,"_clones",sep=""),paste(clone,"_5Days_2.txt",sep="")),chr)
  table_clone_five_days = merge(table_clone_five_days_rep1,table_clone_five_days_rep2,by="Genes")
  table_clone_five_days$mean_ratio = rowMeans(table_clone_five_days[,2:3])
  
  table_clone = merge(table_clone_no_Dox,table_clone_five_days,by="Genes")
  table_clone = data.frame(Genes = table_clone$Genes, no_Dox = table_clone$mean_ratio.x+pseudocount, five_days = table_clone$mean_ratio.y+pseudocount)
  if(chr == "chrX"){table_clone = table_clone[table_clone$no_Dox > 0.2 & table_clone$no_Dox < 0.8,]}
  if(chr == "chr12"){table_clone = table_clone[table_clone$no_Dox > 0.46 & table_clone$no_Dox < 0.86,]}
  
  #table_clone$foldchange = table_clone$five_days/table_clone$no_Dox
  table_clone$foldchange = (table_clone$five_days/(1-table_clone$five_days))*((1-table_clone$no_Dox)/table_clone$no_Dox)
  return(table_clone)
}

####Permutation test
empirical_p_value <- function(B,table_clone_predictions,n,x){
  fraction = rep(0,B)
  for(i in 1:B){
    boot_idx = sample(nrow(table_clone_predictions),n,replace = F)
    boot_table = table_clone_predictions[boot_idx,]
    fraction[i] = nrow(boot_table[boot_table$class == 0,])/nrow(boot_table)
  }
  p_value = sum(fraction > x)/B
  return(list(p_value,fraction))
}

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
#perform permutation test and plot results
###################################################################################

cairo_pdf(here(output_dir,"analysis_paper_loda.pdf"),width = 4,height = 4,onefile = T)

clone_performance = NULL
clone_fc = NULL
bootstrap_performance = as.data.frame(matrix(0,nrow = B,ncol = length(clone)))

for(i in 1:length(clone)){
  print(paste("clone:",clone[i]))
  table_clone = get_foldchange(input_dir,chr[i],clone[i],pseudocount)
  print(paste("genes in clone:",nrow(table_clone)))
  
  predictions = read.table(file = here(input_dir_clones,predictions_clones[i]),header = T)
  table_clone_predictions = merge(table_clone,predictions,by="Genes")
  table_clone_predictions = table_clone_predictions[table_clone_predictions$class == 0 | table_clone_predictions$class == 1,]
  clone_fc = rbind(clone_fc,data.frame(foldchange=table_clone_predictions$foldchange,clone=clone[i],chr=chr[i]))
  print(paste("genes with predictions:",nrow(table_clone_predictions)))
  
  thr_foldchange = 0.9
  table_clone_predictions_fc = table_clone_predictions[table_clone_predictions$foldchange < thr_foldchange,]
  print(paste("genes with fc <",thr_foldchange,"are: ",nrow(table_clone_predictions_fc)))
  
  x = nrow(table_clone_predictions_fc[table_clone_predictions_fc$class == 0,])/nrow(table_clone_predictions_fc)
  p_value = empirical_p_value(B,table_clone_predictions,nrow(table_clone_predictions_fc),x)
  print(paste("accuracy:",round(x,3),"p-value:",p_value[[1]]))
  
  clone_performance = rbind(clone_performance,table(table_clone_predictions_fc$class)/nrow(table_clone_predictions_fc))
  bootstrap_performance[,i] = p_value[[2]]
  
  print(paste("silenced:",nrow(table_clone_predictions_fc[table_clone_predictions_fc$class == 0,])))
  print(paste("not silenced:",nrow(table_clone_predictions_fc[table_clone_predictions_fc$class == 1,])))
  
  cortest = cor.test(table_clone_predictions$vote,table_clone_predictions$foldchange)
  scatterplot_dense_colors(table_clone_predictions$vote,table_clone_predictions$foldchange,"vote(class0)","foldchange",paste("clone",clone[i]))

  table_clone_predictions$class = as.factor(table_clone_predictions$class)

  wilcox =wilcox.test(table_clone_predictions$foldchange[table_clone_predictions$class==0],table_clone_predictions$foldchange[table_clone_predictions$class==1])$p.value
  gg_box = ggplot(table_clone_predictions, aes(x=class,y=foldchange)) + 
    geom_boxplot(colour = "#4d4d4d",alpha = 0.7,outlier.size=0.1,lwd=0.4) + 
    labs(title=paste("clone",clone[i]), subtitle = paste("Wilcox-Test:",signif(wilcox,3))) +
    theme_minimal(base_family = "Source Sans Pro") + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size=8, hjust=1), axis.text.y = element_text(size=8), 
          axis.title=element_text(size=8),plot.title = element_text(size=8),plot.subtitle = element_text(size=7)) +
    scale_x_discrete(name = "silencing class") + scale_y_continuous(name = "folchange",limits = c(0,2))
  print(gg_box)

}

###plot histogram and boxplots
colnames(bootstrap_performance) = clone
performance = data.frame(ind = factor(clone), y = clone_performance[,1]*100)
mean = data.frame(ind = factor(clone), y = colMeans(bootstrap_performance)*100)

ggplot(stack(bootstrap_performance*100), aes(x = ind, y = values)) +  
  geom_boxplot(notch=FALSE,fill = "lightgrey", colour = "grey",alpha = 0.7,outlier.shape = "",width=0.5) + 
  geom_jitter(position=position_jitter(h = 1,w=0.25),size=0.01,color="grey") +
  stat_summary(fun.y="mean", geom="point", color = 'black',pch="_",size=18) +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=8), 
        axis.title=element_text(size=8), plot.title = element_text(size=8)) + 
  scale_x_discrete(name = "clone") + scale_y_continuous(name = "genes predicted as silenced (%)", limits = c(20,90)) +
  geom_point(data = data.frame(x = factor(clone), y = clone_performance[,1]*100),aes(x=x, y=y, size=15),color = 'red',pch="_",size=18)
dev.off()

#plot histogram of clone predictions
cairo_pdf(here(output_dir,"paper_figures_loda_histogram.pdf"),width = 2.7,height = 2.7,onefile = T)
ggplot(stack(bootstrap_performance*100), aes(x = values)) + 
  geom_histogram(binwidth=0.9, colour="black", fill="white",lwd=0.3) + 
  facet_grid(ind ~ .) +
  geom_vline(aes(xintercept = y),performance, size=0.4, colour="red") +
  geom_vline(aes(xintercept = y),mean, size=0.4, colour="black",linetype="dashed") +
  scale_y_continuous(breaks=c(0,300),name = "# permutations") + scale_x_continuous(name="genes predicted as silenced (%)") +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=8), 
        axis.title=element_text(size=8), plot.title = element_text(size=8))
dev.off()

#print AER distribution for each clone
cairo_pdf(here(output_dir,"paper_figures_loda_AER_density.pdf"),width = 5,height = 3,onefile = T)
ggplot(clone_fc, aes(x=foldchange,colour=clone)) +
  geom_density(aes(linetype=clone),data=clone_fc,alpha=.4, lwd=0.5) + 
  scale_x_continuous(limits=c(0,2), name='normalized allelic expression ratio (AER)') +
  scale_y_continuous(name='distribution of normalized AER') +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=8), 
        axis.title=element_text(size=8), plot.title = element_text(size=8), legend.title = element_text(size=8), 
        legend.text = element_text(size=8)) 
dev.off()




