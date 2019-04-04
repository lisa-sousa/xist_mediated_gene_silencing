###################################################################################
#libraries
###################################################################################

library(rtracklayer)
library(GenomicRanges)
library(gplots)
library(Cairo)
library(ggplot2)
library(gridExtra)
library(cowplot)

###################################################################################
#input files
###################################################################################

input_dir_data = '/project/lncrna/Xist/data/silencing_halftimes/raw_data/'
input_dir_gene_annotation = '/project/lncrna/Xist/data/annotation_files/gene_annotation/'

file_pro_seq_raw = paste(input_dir_data,'PROseq.txt',sep='')
file_mrna_seq_undiff_raw = paste(input_dir_data,'UndifferenciatedCells.txt',sep='')
file_mrna_seq_diff_raw = paste(input_dir_data,'DifferenciatedCells.txt',sep='')
file_gencode_gene_annotation = paste(input_dir_gene_annotation,'gencode.vM9.annotation.chrX.genes.bed',sep='')

###################################################################################
#output files
###################################################################################

output_dir_plot = '/project/lncrna/Xist/plots/silencing_halftimes/'
output_dir_data = '/project/lncrna/Xist/data/silencing_halftimes/fitted_data/'

file_pro_seq_all = paste(output_dir_data,'halftimes_pro_seq_mm10_RSS_initial_ratio.txt',sep='')
file_pro_seq_halftimes = paste(output_dir_data,'halftimes_pro_seq_mm10.bed',sep='')
file_mrna_seq_undiff_all = paste(output_dir_data,'halftimes_mrna_seq_undiff_mm10_RSS_initial_ratio.txt',sep='')
file_mrna_seq_undiff_halftimes = paste(output_dir_data,'halftimes_mrna_seq_undiff_mm10.bed',sep='')
file_mrna_seq_diff_all = paste(output_dir_data,'halftimes_mrna_seq_diff_mm10_RSS_initial_ratio.txt',sep='')
file_mrna_seq_diff_halftimes = paste(output_dir_data,'halftimes_mrna_seq_diff_mm10.bed',sep='')

###################################################################################
#parameters files
###################################################################################

control=nls.control(warnOnly=TRUE) # prevent stopping the fit
pseudo_count = 0.001

#get gencode vM9 gene annotation
gene_annotation = read.table(file=file_gencode_gene_annotation)
colnames(gene_annotation) = c('Chromosomes','Start','End','Genes','score','Strand')
gene_annotation = gene_annotation[!duplicated(gene_annotation$Genes),]

###################################################################################
#functions
###################################################################################

#eponential decay function
fitting_function <- function(t, start) {
  return(exp(-start*t))
}

#normalize data and fit half-times with exponential decay function
fit_data <- function(data,pseudo_count,time,fitted_data,control){
  for (i in 1:nrow(data)) {
    
    # Transform the data such that you can fit relative reads from BL6
    data_trans = as.numeric(data[i,])
    data_trans = ((1-data_trans[1])/data_trans[1]) * (data_trans/(1-data_trans+pseudo_count))
    data_trans = data.frame(time = time, expression = data_trans)
    
    fitnl=nls(expression ~ fitting_function(time,k), data=data_trans, start=list(k=1),control=control, algorithm='port',lower=c(k=1/5))
    k = abs(coef(fitnl)[1])
    fitted_data$halftime[i] = log(2)/k
    fitted_data$RSS[i] = sqrt(sum(resid(fitnl)^2))
    
    plot_fit(data_trans,k,as.character(fitted_data$Genes[i]))
    
  }
  return(fitted_data)
}

#plot fitted data
plot_fit <- function(data_trans,k,gene){
  halftime_label = bquote(t[1/2] == .(as.double(round(log(2)/k,2)))~d)
  ggplot = ggplot(data_trans, aes(x=time, y=expression)) + 
    geom_point() + 
    stat_function(fun=fitting_function, args = list(k)) + 
    geom_text(aes(x=0.1,y=0.1),label=deparse(halftime_label),parse=T,size=2.5,hjust=0) +
    theme_minimal(base_family = "Source Sans Pro") + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=6), axis.title=element_text(size=7), 
          plot.title=element_text(size=8,hjust = 0.5,face="italic")) + 
    scale_x_continuous(limits=c(0,max(data_trans$time)),breaks=c(0,max(data_trans$time)/2,max(data_trans$time)), name='time [days]')  + 
    scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.5,1), name='norm. B6 expression') + 
    ggtitle(gene)
  print(ggplot)
}

#plot raw data sorted by genomic position
plot_ratios <- function(data){
  
  data_sorted = data[order(data$Start),]
  matrix = as.matrix(data_sorted[,grep('Ratio',colnames(data_sorted))])
  rownames(matrix) = data_sorted$Genes
  colnames(matrix)=gsub("Dox30min[A-z,0-9]*","Dox 0.5h",gsub("NoDox([A-z]{1})[A-z,0-9]*","No Dox \\1",gsub("h([A-z]+)","h \\1",gsub("Dox([0-9]*)hr([A-z]?)_[A-z,0-9]*","Dox \\1h\\2",colnames(matrix),perl = T))))
  matrix = matrix[,ncol(matrix):1]
  
  require(reshape2)
  matrix = melt(t(matrix))
  
  # Create the heatmap
  ggplot = ggplot(matrix, aes(Var2, Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5, limit = c(0,1), 
                         space = "Lab",breaks = c(0,0.5,1),name="fraction of reads\nexpressed from \nB6 allele") + # Change gradient color
    theme_minimal(base_family = "Source Sans Pro",base_size=8) +
    theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), legend.position="bottom", legend.margin = margin(t = 0, unit='cm'), 
          legend.text = element_text(size=7), legend.title = element_text(size=8), axis.title=element_text(size=8),plot.background=element_blank(),panel.border=element_blank())+
    coord_fixed(ratio=5) + xlab("genes (sorted by genomic location)") + ylab("time point") + 
    scale_y_discrete(position = "right",expand = c(0, 0)) + scale_x_discrete(breaks = c("Tsix"),labels = expression(italic(Xist)~'locus'),expand = c(0, 0)) +
    geom_vline(aes(xintercept = which(data_sorted$Genes=="Tsix")),size=0.7,alpha=0.6)
  return(ggplot)
}

###################################################################################
####fit PRO-Seq data
###################################################################################

pro_seq_data = read.table(file_pro_seq_raw,na.strings=c('NA','nd'),header=T,sep='\t')
time = c(0,0.5,1,2,4,8,12,24)/24

#remove rows with missing values and filter for X-chromosomal genes
pro_seq_data = pro_seq_data[complete.cases(pro_seq_data),]
pro_seq_data_chrX = pro_seq_data[pro_seq_data$Chromosomes=='chrX',]

#filter for genes that have at least 10 reads for every timestep
total_counts_data = pro_seq_data_chrX[,grep('total',colnames(pro_seq_data_chrX))]
expressed_genes = which(rowSums(total_counts_data < 10) == 0)

#filter for genes that have a more or less balanced initial ratio 
initial_ratio = rowMeans(pro_seq_data_chrX[,colnames(pro_seq_data_chrX) %in% c('NoDoxA_BL6_Ratio','NoDoxB_BL6_Ratio')])
balanced_initial_ratio_genes = which((initial_ratio > 0.2) & (initial_ratio < 0.8))

#use only expressed genes with balanced initial ratio
pro_seq_genes = intersect(expressed_genes,balanced_initial_ratio_genes)
pro_seq_data_chrX = pro_seq_data_chrX[pro_seq_genes,]

pro_seq_gene_ratios = pro_seq_data_chrX[,grep('Ratio',colnames(pro_seq_data_chrX))]
pro_seq_gene_ratios = data.frame(NoDox_BL6_Ratio = rowMeans(pro_seq_gene_ratios[,1:2]),pro_seq_gene_ratios[,3:8],Dox24hr_BL6_Ratio = rowMeans(pro_seq_gene_ratios[,9:10]))

#fit the data
cairo_pdf(paste(output_dir_plot,'pro_seq_fit.pdf',sep=''),width = 2,height = 2, onefile = TRUE)
fitted_data_pro_seq = data.frame(pro_seq_data_chrX[,c(2,4,5,1)],halftime = rep(-1,length(pro_seq_genes)),Strand = pro_seq_data_chrX[,3],
                                 RPKM = rowMeans(pro_seq_data_chrX[,c(10,15)]),initial_ratio = initial_ratio[pro_seq_genes],RSS = rep(-1,length(pro_seq_genes)))
fitted_data_pro_seq = fit_data(pro_seq_gene_ratios,pseudo_count,time,fitted_data_pro_seq,control)
dev.off()

#filter for not well fitted genes --> residuals
fitted_data_pro_seq_filtered = fitted_data_pro_seq[fitted_data_pro_seq$RSS < 1.5,]

gene_regions_pro_seq = GRanges(seqnames = fitted_data_pro_seq_filtered$Chromosomes, ranges = IRanges(fitted_data_pro_seq_filtered$Start,fitted_data_pro_seq_filtered$End), 
                               strand = fitted_data_pro_seq_filtered$Strand)
mcols(gene_regions_pro_seq) = data.frame(name = fitted_data_pro_seq_filtered$Genes, score = fitted_data_pro_seq_filtered$halftime)

#write data to output
write.table(fitted_data_pro_seq_filtered,file=file_pro_seq_all,sep='\t',col.names=T,row.names=F,quote=F)
export.bed(gene_regions_pro_seq, con=file_pro_seq_halftimes, format='bed')	

cairo_pdf(paste(output_dir_plot,'pro_seq_ratios.pdf',sep=''),width = 6,height = 2)
print(plot_ratios(pro_seq_data_chrX))
dev.off()

###################################################################################
####fit undifferentiated mRNA-Seq data
###################################################################################

mrna_seq_undiff_data = read.table(file_mrna_seq_undiff_raw,na.strings=c('NA','nd'),header=T,sep='\t')
time = c(0,2,4,8,12,24)/24

#filter out second replicate for each time point --> did not work very well
mrna_seq_undiff_data = mrna_seq_undiff_data[,c(1,2,3,5,4,8,12,16,20,24)]
colnames(mrna_seq_undiff_data) = c('Genes','chr','strand','RPKM','NoDox_Ratio','Dox2hr_Ratio','Dox4hr_Ratio','Dox8hr_Ratio','Dox12hr_Ratio','Dox24hr_Ratio')

#remove rows with missing values and filter for X-chromosomal genes
mrna_seq_undiff_data = mrna_seq_undiff_data[complete.cases(mrna_seq_undiff_data),]
mrna_seq_undiff_data_chrX = mrna_seq_undiff_data[mrna_seq_undiff_data$chr=='chrX',]

#filter for genes that have a more or less balanced initial ratio 
initial_ratio = mrna_seq_undiff_data_chrX$NoDox_Ratio
mrna_seq_undiff_genes = which((initial_ratio > 0.2) & (initial_ratio < 0.8))

#use only genes with balanced initial ratio
mrna_seq_undiff_data_chrX = mrna_seq_undiff_data_chrX[mrna_seq_undiff_genes,]
mrna_seq_undiff_data_chrX = merge(gene_annotation,mrna_seq_undiff_data_chrX,by='Genes')
mrna_seq_undiff_gene_ratios = mrna_seq_undiff_data_chrX[,grep('Ratio',colnames(mrna_seq_undiff_data_chrX))]

#fit the data
cairo_pdf(paste(output_dir_plot,'mrna_seq_undiff_fit.pdf',sep=''),width = 2,height = 2, onefile = TRUE)
fitted_data_mrna_seq_undiff = data.frame(mrna_seq_undiff_data_chrX[,c(2,3,4,1)],halftime = rep(-1,length(mrna_seq_undiff_genes)),Strand = mrna_seq_undiff_data_chrX[,6],
                                  RPKM = rowMeans(mrna_seq_undiff_data_chrX[,c(10,15)]),initial_ratio = initial_ratio[mrna_seq_undiff_genes],RSS = rep(-1,length(mrna_seq_undiff_genes)))
fitted_data_mrna_seq_undiff = fit_data(mrna_seq_undiff_gene_ratios,pseudo_count,time,fitted_data_mrna_seq_undiff,control)
dev.off()

#filter for not well fitted genes --> residuals
fitted_data_mrna_seq_undiff_filtered = fitted_data_mrna_seq_undiff[fitted_data_mrna_seq_undiff$RSS < 1.5,]

gene_regions_mrna_seq_undiff = GRanges(seqnames = fitted_data_mrna_seq_undiff_filtered$Chromosomes, ranges = IRanges(fitted_data_mrna_seq_undiff_filtered$Start,fitted_data_mrna_seq_undiff_filtered$End), 
                               strand = fitted_data_mrna_seq_undiff_filtered$Strand)
mcols(gene_regions_mrna_seq_undiff) = data.frame(name = fitted_data_mrna_seq_undiff_filtered$Genes, score = fitted_data_mrna_seq_undiff_filtered$halftime)

#write data to output
write.table(fitted_data_mrna_seq_undiff_filtered,file=file_mrna_seq_undiff_all,sep='\t',col.names=T,row.names=F,quote=F)
export.bed(gene_regions_mrna_seq_undiff, con=file_mrna_seq_undiff_halftimes, format='bed')	

#plot ratios
mrna_seq_undiff_data = read.table(file_mrna_seq_undiff_raw,na.strings=c('NA','nd'),header=T,sep='\t')
mrna_seq_undiff_data = mrna_seq_undiff_data[,-grep('_RPKM',colnames(mrna_seq_undiff_data))]
colnames(mrna_seq_undiff_data)[4:15] = c('NoDoxA_Ratio','NoDoxB_Ratio','Dox2hrA_Ratio','Dox2hrB_Ratio','Dox4hrA_Ratio','Dox4hrB_Ratio','Dox8hrA_Ratio','Dox8hrB_Ratio',
                                         'Dox12hrA_Ratio','Dox12hrB_Ratio','Dox24hrA_Ratio','Dox24hrB_Ratio')
#remove rows with missing values and filter for X-chromosomal genes
mrna_seq_undiff_data = mrna_seq_undiff_data[complete.cases(mrna_seq_undiff_data),]
mrna_seq_undiff_data_chrX = mrna_seq_undiff_data[mrna_seq_undiff_data$Chromosomes=='chrX',]
#filter for genes that have a more or less balanced initial ratio 
initial_ratio = rowMeans(mrna_seq_undiff_data[,colnames(mrna_seq_undiff_data) %in% c('NoDoxA_Ratio','NoDoxB_Ratio')])
mrna_seq_undiff_genes = which((initial_ratio > 0.2) & (initial_ratio < 0.8))
#use only genes with balanced initial ratio
mrna_seq_undiff_data_chrX = mrna_seq_undiff_data_chrX[mrna_seq_undiff_genes,]
mrna_seq_undiff_data_chrX = merge(gene_annotation,mrna_seq_undiff_data_chrX,by='Genes')

cairo_pdf(paste(output_dir_plot,'mrna_seq_undiff_ratios.pdf',sep=''),width = 6,height = 2)
print(plot_ratios(mrna_seq_undiff_data_chrX))
dev.off()

###################################################################################
####fit differentiated mRNA-Seq data
###################################################################################

mrna_seq_diff_data = read.table(file_mrna_seq_diff_raw,na.strings=c('NA','nd'),header=T,sep='\t')
time = c(0,8,16,24,48)/24

#only keep expression ratios
mrna_seq_diff_data = mrna_seq_diff_data[,c(1,2,3,5,7,4,6,8,10,12,14,16,18,20,22)]
colnames(mrna_seq_diff_data) = c('Genes','chr','strand','NoDoxA_RPKM','NoDoxB_RPKM','NoDoxA_Ratio','NoDoxB_Ratio','Dox8hrA_Ratio','Dox8hrB_Ratio','Dox16hrA_Ratio','Dox16hrB_Ratio',
                                       'Dox24hrA_Ratio','Dox24hrB_Ratio','Dox48hrA_Ratio','Dox48hrB_Ratio')

#remove rows with missing values and filter for X-chromosomal genes
mrna_seq_diff_data = mrna_seq_diff_data[complete.cases(mrna_seq_diff_data),]
mrna_seq_diff_data_chrX = mrna_seq_diff_data[mrna_seq_diff_data$chr=='chrX',]

#filter for genes that have a more or less balanced initial ratio 
initial_ratio = rowMeans(mrna_seq_diff_data_chrX[,colnames(mrna_seq_diff_data_chrX) %in% c('NoDoxA_Ratio','NoDoxB_Ratio')])
mrna_seq_diff_genes = which((initial_ratio > 0.2) & (initial_ratio < 0.8))

#use only genes with balanced initial ratio
mrna_seq_diff_data_chrX = mrna_seq_diff_data_chrX[mrna_seq_diff_genes,]
mrna_seq_diff_data_chrX = merge(gene_annotation,mrna_seq_diff_data_chrX,by='Genes')
mrna_seq_diff_gene_ratios = mrna_seq_diff_data_chrX[,grep('Ratio',colnames(mrna_seq_diff_data_chrX))]
mrna_seq_diff_gene_ratios = data.frame(NoDox_Ratio = rowMeans(mrna_seq_diff_gene_ratios[,1:2]),Dox8hr_Ratio = rowMeans(mrna_seq_diff_gene_ratios[,3:4]),
                                       Dox16hr_Ratio = rowMeans(mrna_seq_diff_gene_ratios[,5:6]),Dox24hr_Ratio = rowMeans(mrna_seq_diff_gene_ratios[,7:8]),
                                       Dox48hr_Ratio = rowMeans(mrna_seq_diff_gene_ratios[,9:10]))

#fit the data
cairo_pdf(paste(output_dir_plot,'mrna_seq_diff_fit.pdf',sep=''),width = 2,height = 2, onefile = TRUE)
fitted_data_mrna_seq_diff = data.frame(mrna_seq_diff_data_chrX[,c(2,3,4,1)],halftime = rep(-1,length(mrna_seq_diff_genes)),Strand = mrna_seq_diff_data_chrX[,6],
                                  RPKM = rowMeans(mrna_seq_diff_data_chrX[,c(9,10)]),initial_ratio = initial_ratio[mrna_seq_diff_genes],RSS = rep(-1,length(mrna_seq_diff_genes)))
fitted_data_mrna_seq_diff = fit_data(mrna_seq_diff_gene_ratios,pseudo_count,time,fitted_data_mrna_seq_diff,control)
dev.off()

#filter for not well fitted genes --> residuals
fitted_data_mrna_seq_diff_filtered = fitted_data_mrna_seq_diff[fitted_data_mrna_seq_diff$RSS < 1.5,]

gene_regions_mrna_seq_diff = GRanges(seqnames = fitted_data_mrna_seq_diff_filtered$Chromosomes, ranges = IRanges(fitted_data_mrna_seq_diff_filtered$Start,fitted_data_mrna_seq_diff_filtered$End), 
                                       strand = fitted_data_mrna_seq_diff_filtered$Strand)
mcols(gene_regions_mrna_seq_diff) = data.frame(name = fitted_data_mrna_seq_diff_filtered$Genes, score = fitted_data_mrna_seq_diff_filtered$halftime)

#write data to output
write.table(fitted_data_mrna_seq_diff_filtered,file=file_mrna_seq_diff_all,sep='\t',col.names=T,row.names=F,quote=F)
export.bed(gene_regions_mrna_seq_diff, con=file_mrna_seq_diff_halftimes, format='bed')	


cairo_pdf(paste(output_dir_plot,'mrna_seq_diff_ratios.pdf',sep=''),width = 6,height = 2)
print(plot_ratios(mrna_seq_diff_data_chrX))
dev.off()

###################################################################################
####plots for paper
###################################################################################

ggplot = list()

#### plot example for fitting half-times
time = c(0,0.5,1,2,4,8,12,24)/24
expression = c(1,0.9,0.85,0.8,0.7,0.6,0.5,0.4)
data_trans =  data.frame(time = time, expression = expression)
fitnl=nls(expression ~ fitting_function(time,k), data=data_trans, start=list(k=1),control=control, algorithm='port',lower=c(k=1/5))
k = abs(coef(fitnl)[1])

halftime_label = bquote(t[1/2] == ln(2)/k)
function_label  = expression(paste("N(t)=e"^"-k*t"))
ggplot[[1]] = ggplot(data_trans, aes(x=time, y=expression)) + 
  geom_point(aes(colour="data"),size=1) + 
  stat_function(aes(colour="fitted curve"),fun=fitting_function, args = list(k)) + 
  geom_text(aes(x=1,y=0.6),label=deparse(halftime_label),parse=T,size=2.5,hjust=1,vjust=0) +
  geom_text(aes(x=1,y=0.9),label=function_label,size=2.5,hjust=1,vjust=1) +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=8), axis.title=element_text(size=8), 
        legend.text = element_text(size=8),legend.position = "top") + 
  scale_x_continuous(limits=c(0,max(data_trans$time)),breaks=c(0,max(data_trans$time)/2,max(data_trans$time)), name='time [days]')  + 
  scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.5,1), name='norm. B6 expression') +
  geom_segment(aes(x = 0, y = 0.5, xend = 0.52, yend = 0.5), linetype="dashed",size = 0.3) +
  geom_segment(aes(x = 0.52, y = 0.5, xend = 0.52, yend = 0), linetype="dashed",size = 0.3) +
  scale_colour_manual(name='', values=c('data'='#71c837', 'fitted curve'='black', guide='legend')) +
  guides(fill = guide_legend(override.aes = list(linetype = 1, shape=''),nrow=3), colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16,NA)),nrow=3))

#####plot example genes for fitting
example_genes = c("Otud5", "Piga", "Stard8")
data_trans_all = pro_seq_gene_ratios[pro_seq_data_chrX$Genes %in% example_genes,]
data_trans_all = ((1-data_trans_all[,1])/data_trans_all[,1]) * (data_trans_all/(1-data_trans_all+pseudo_count))
k = c()

for (i in 1:length(example_genes)) {
  data_trans = data.frame(time = time, expression = as.numeric(data_trans_all[i,]))
  fitnl=nls(expression ~ fitting_function(time,k), data=data_trans, start=list(k=1),control=control, algorithm='port',lower=c(k=1/5))
  k[i] = abs(coef(fitnl)[1])
}

data_trans_all$gene = example_genes
plot_data = melt(data_trans_all,id.vars = "gene", measure.vars = colnames(data_trans_all)[1:(ncol(data_trans_all)-1)])
plot_data$variable = rep(time,each=3)

halftime_label = c(bquote(italic(.(example_genes[1]))~":"~t[1/2] == .(as.double(round(log(2)/k[1],2)))~d),
                   bquote(italic(.(example_genes[2]))~":"~t[1/2] == .(as.double(round(log(2)/k[2],2)))~d),
                   bquote(italic(.(example_genes[3]))~":"~t[1/2] == .(as.double(round(log(2)/k[3],2)))~d))

ggplot[[2]] = ggplot(data=plot_data, aes(x=variable, y=value, shape=gene, colour=gene)) +
  geom_point() +
  scale_colour_manual(name = "",values=c("#16502d", "#00aa88", "#55ddff"), label=halftime_label) +
  scale_shape_manual(name = "",values=c(16,17,18), label=halftime_label) +
  stat_function(fun=fitting_function, args = list(k[1]),color="#16502d",linetype="solid") +
  stat_function(fun=fitting_function, args = list(k[2]),color="#00aa88",linetype="dashed") +
  stat_function(fun=fitting_function, args = list(k[3]),color="#55ddff",linetype="dotted") +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=8), axis.title=element_text(size=8), 
        legend.text = element_text(size=8),legend.position = "top") +
  scale_x_continuous(limits=c(0,max(time)),breaks=c(0,max(time)/2,max(time)), name='time [days]') +
  scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.5,1), name='norm. B6 expression') +
  guides(fill=guide_legend(nrow=3), col=guide_legend(nrow=3))

#####plot ratio chrX vs autosomes
pro_seq_data$Chromosomes = factor(pro_seq_data$Chromosomes, levels = c("autosomes",levels(pro_seq_data$Chromosomes)))
pro_seq_data$Chromosomes[pro_seq_data$Chromosomes != "chrX"] = "autosomes"

ratios = pro_seq_data[,grep('Ratio',colnames(pro_seq_data))]
plot_data = data.frame(chromosome = pro_seq_data$Chromosomes, NoDox_BL6_Ratio = rowMeans(ratios[,1:2]), ratios[,3:8], Dox24hr_BL6_Ratio = rowMeans(ratios[,9:10]))


colnames(plot_data)[2:9]=gsub("Dox30min[A-z,0-9]*","0.5h",gsub("NoDox([A-z]{1})[A-z,0-9]*","0h",gsub("h([A-z]+)","h \\1",gsub("Dox([0-9]*)hr([A-z]?)_[A-z,0-9]*","\\1h\\2",colnames(plot_data)[2:9],perl = T))))
require(reshape2)
plot_data = melt(plot_data)

ggplot[[3]] = ggplot(plot_data,aes(y=value, x=variable, fill=chromosome)) + 
  stat_boxplot(geom ='errorbar',lwd=0.3) +
  geom_boxplot(outlier.size=-1,lwd=0.4) +
  scale_fill_grey(start=0.4, end=0.7, label=c("autosomes","X-chromosome")) +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=8), axis.title=element_text(size=8, margin = margin(t = 0)), axis.text.x = element_text(size=7,angle = 45,margin = margin(t = 0, b = 0)), 
        legend.text = element_text(size=8),legend.position = "top", legend.title = element_blank()) +
  scale_x_discrete(name='timepoints') +
  scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.5,1), name='fraction B6 reads') +
  guides(fill=guide_legend(nrow=3), col=guide_legend(nrow=3))



cairo_pdf(paste(output_dir_plot,'paper_example_fitting.pdf',sep=''),width = 6,height = 3, onefile = TRUE)
grid.arrange(grobs=ggplot, ncol=3)
dev.off()

#####plot all ratios
pro = plot_ratios(pro_seq_data_chrX)
mrna_und = plot_ratios(mrna_seq_undiff_data_chrX)
mrna_diff = plot_ratios(mrna_seq_diff_data_chrX)

pro = pro + theme(legend.position="none")
mrna_und = mrna_und + theme(legend.position="none")
mrna_diff = mrna_diff + theme(legend.position="bottom")
legend = get_legend(mrna_diff)
mrna_diff = mrna_diff + theme(legend.position="none")

cairo_pdf(paste(output_dir_plot,'paper_ratios.pdf',sep=''),width = 7,height = 5, onefile = TRUE)
grid.arrange(pro,mrna_und,mrna_diff,legend,ncol=1,heights=c(1.3,1.56,1.3,0.5))
dev.off()










