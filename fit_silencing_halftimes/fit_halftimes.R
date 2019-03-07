###################################################################################
#libraries
###################################################################################

library(rtracklayer)
library(GenomicRanges)
library(gplots)
library(Cairo)

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
    
    plot(data_trans$time,data_trans$expression,xlim=c(-0.2, 1.2),ylim=c(0,1.2),xlab='time [days]', ylab='BL6',main=fitted_data$Genes[i])
    lines(seq.int(0,1,0.01), fitting_function(seq.int(0,1,0.01),k),lwd=2)
    legend('topright',legend = paste('square root RSS=',round(sqrt(sum(resid(fitnl)^2)),2),'\nhalftime=',round(log(2)/k,2)))
  }
  return(fitted_data)
}

#plot raw data sorted by genomic position
plot_ratios <- function(data){

  data_sorted = data[order(data$Start),]
  matrix = as.matrix(data_sorted[,grep('Ratio',colnames(data_sorted))])
  rownames(matrix) = data_sorted$Genes
  
  #define colors ranges 
  idx_low = round(min(matrix),2)*100
  idx_high = round(max(matrix),2)*100
  color = bluered(100)[idx_low:idx_high]
  
  par(oma=c(8,1,1,1))
  heatmap.2(x=matrix, Rowv=NULL,Colv=NULL, col = color, scale='none', margins=c(3,0), trace='none', symkey=FALSE, symbreaks=FALSE, dendrogram='none', density.info='histogram', 
            denscol='black', keysize=1, key.par=list(mar=c(3.5,0,3,0)), lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(1, 5), lwid=c(3, 10, 3), cexCol = 1,cexRow = 0.25)
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
CairoPDF(paste(output_dir_plot,'pro_seq_fit.pdf',sep=''),width = 7,height = 7)
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

CairoPDF(paste(output_dir_plot,'pro_seq_ratios.pdf',sep=''),width = 10,height = 20)
plot_ratios(pro_seq_data_chrX)
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
CairoPDF(paste(output_dir_plot,'mrna_seq_undiff_fit.pdf',sep=''),width = 7,height = 7)
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

CairoPDF(paste(output_dir_plot,'mrna_seq_undiff_ratios.pdf',sep=''),width = 10,height = 20)
plot_ratios(mrna_seq_undiff_data_chrX)
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
CairoPDF(paste(output_dir_plot,'mrna_seq_diff_fit.pdf',sep=''),width = 7,height = 7)
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


CairoPDF(paste(output_dir_plot,'mrna_seq_diff_ratios.pdf',sep=''),width = 10,height = 20)
plot_ratios(mrna_seq_diff_data_chrX)
dev.off()
