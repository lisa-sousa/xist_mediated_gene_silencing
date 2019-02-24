
get_clone_table <- function(file, chr){
  table_clone = read.table(file=file,header = T)
  table_clone = table_clone[grep(chr,table_clone$coord),]
  table_clone = table_clone[,c(1,6)]
  colnames(table_clone) = c("Genes","ratio")
  return(table_clone)
}

get_foldchange <- function(input_dir,chr,clone,pseudocount){
  table_clone_no_Dox_rep1 = get_clone_table(paste(input_dir,chr,"_clones/",clone,"_noDox_1.txt",sep=""),chr)
  table_clone_no_Dox_rep2 = get_clone_table(paste(input_dir,chr,"_clones/",clone,"_noDox_2.txt",sep=""),chr)
  table_clone_no_Dox = merge(table_clone_no_Dox_rep1,table_clone_no_Dox_rep2,by="Genes")
  table_clone_no_Dox$mean_ratio = rowMeans(table_clone_no_Dox[,2:3])
  
  table_clone_five_days_rep1 = get_clone_table(paste(input_dir,chr,"_clones/",clone,"_5Days_1.txt",sep=""),chr)
  table_clone_five_days_rep2 = get_clone_table(paste(input_dir,chr,"_clones/",clone,"_5Days_2.txt",sep=""),chr)
  table_clone_five_days = merge(table_clone_five_days_rep1,table_clone_five_days_rep2,by="Genes")
  table_clone_five_days$mean_ratio = rowMeans(table_clone_five_days[,2:3])
  
  table_clone = merge(table_clone_no_Dox,table_clone_five_days,by="Genes")
  table_clone = data.frame(Genes = table_clone$Genes, no_Dox = table_clone$mean_ratio.x+pseudocount, five_days = table_clone$mean_ratio.y+pseudocount)
  if(chr == "chrX"){table_clone = table_clone[table_clone$no_Dox > 0.2 & table_clone$no_Dox < 0.8,]}
  if(chr == "chr12"){table_clone = table_clone[table_clone$no_Dox > 0.46 & table_clone$no_Dox < 0.86,]}
  
  table_clone$foldchange = table_clone$five_days/table_clone$no_Dox
  return(table_clone)
}

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
  plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, xlim = c(0,1), ylim = c(0,2), ylab = ylab, xlab = xlab,main = main)
  abline(lm(x2~x1))
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
  return(list(p_value,mean(fraction)))
}



output_dir = "/project/lncrna/Xist/plots/additional_analysis/"
input_dir_clones = "/project/lncrna/Xist/data/modelling/model/clones/"
predictions_clones = c("predictions_clones_86.txt","predictions_clones_87.txt","predictions_clones_109.txt","predictions_clones_190.txt","predictions_clones_228.txt","predictions_clones_273.txt")

input_dir = "/project/lncrna/Xist/data/annotation_files/xist_transgenes/"
clone = c("86","87","109","190","228","273")
chr = c("chrX","chrX","chrX","chrX","chr12","chr12")
pseudocount = 0.001
thr_foldchange = c()
B=1000

clone_performance = NULL
clone_fc = list()

pdf(paste(output_dir,"analysis_paper_loda.pdf",sep=""))

for(i in 1:length(clone)){
  table_clone = get_foldchange(input_dir,chr[i],clone[i],pseudocount)
  # if(chr[i]=="chr12"){
  #   x = table_clone$foldchange
  #   scaled_x = (x-0.5)/(1-0.5)
  #   table_clone$foldchange = scaled_x
  # }
  print(paste("genes in clone:",nrow(table_clone)))
  
  predictions = read.table(file = paste(input_dir_clones,predictions_clones[i],sep=""),header = T)
  table_clone_predictions = merge(table_clone,predictions,by="Genes")
  table_clone_predictions = table_clone_predictions[table_clone_predictions$class == 0 | table_clone_predictions$class == 1,]
  clone_fc[[i]] = table_clone_predictions$foldchange
  print(paste("genes with predictions:",nrow(table_clone_predictions)))
  
  thr_foldchange = as.numeric(quantile(table_clone_predictions$foldchange[table_clone_predictions$foldchange<1],0.85))
  table_clone_predictions_fc = table_clone_predictions[table_clone_predictions$foldchange < thr_foldchange,]
  print(paste("genes with fc <",thr_foldchange,"are: ",nrow(table_clone_predictions_fc)))
  
  x = nrow(table_clone_predictions_fc[table_clone_predictions_fc$class == 0,])/nrow(table_clone_predictions_fc)
  p_value = empirical_p_value(B,table_clone_predictions,nrow(table_clone_predictions_fc),x)
  print(paste("accuracy:",round(x,3),"p-value:",p_value[[1]]))
  
  clone_performance = rbind(clone_performance,table(table_clone_predictions_fc$class)/nrow(table_clone_predictions_fc))
  #clone_performance = rbind(clone_performance,c(p_value[[2]],1-p_value[[2]]))
  
  print(paste("silenced:",nrow(table_clone_predictions_fc[table_clone_predictions_fc$class == 0,])))
  print(paste("not silenced:",nrow(table_clone_predictions_fc[table_clone_predictions_fc$class == 1,])))
  
  scatterplot_dense_colors(table_clone_predictions$vote,table_clone_predictions$foldchange,"vote (class 0)","foldchange",paste("clone",clone[i]))
  print(cor.test(table_clone_predictions$vote,table_clone_predictions$foldchange))
  
  table_clone_predictions$class = as.factor(table_clone_predictions$class)
  
  wilcox = wilcox.test(table_clone_predictions$foldchange[table_clone_predictions$class==0],table_clone_predictions$foldchange[table_clone_predictions$class==1])$p.value
  gg_box = ggplot(table_clone_predictions, aes(x=class,y=foldchange)) + geom_boxplot(notch=FALSE,fill = "lightgrey", colour = "black",alpha = 0.7,outlier.shape = 16) + ggtitle(paste("clone",clone[i],"Wilcox-Test:",signif(wilcox,3))) + theme_bw() + 
    theme(axis.text.x=element_text(hjust = 0.5,size=10),axis.text.y=element_text(size=10),axis.title = element_text(face="bold", size=15),
          plot.title = element_text(hjust = 0.5,size=15,face='bold'),plot.margin = unit(c(2,2,2,2), "cm"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
    scale_x_discrete(name = "silencing class") + scale_y_continuous(name = "folchange")
  print(gg_box)
  
}


row.names(clone_performance) = clone
par(mfrow=c(1,1),mar=c(2,2,2,2),oma=c(5,5,5,5))
barplot(t(clone_performance),main="performance of clones",legend = c("true","false"),ylab="fraction of genes (%)",xlab="clone",las=2)

col = brewer.pal(6,"Set2")
max = c()
for(i in 1:length(clone)){
  density = density(clone_fc[[i]])
  max[i] = max(density$y)
}

plot(density(clone_fc[[1]]),col=col[1],xlab = "foldchange",main="foldchange distribution per clone",ylim=c(0,max(max)))
for(i in 2:length(clone)){
  lines(density(clone_fc[[i]]),col=col[i])
}
legend("topleft",legend = clone,col = col,lty=1)

dev.off()


