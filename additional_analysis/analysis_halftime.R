library(Cairo)

######
#pro-seq
######
halftimes_file = '/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_pro_seq_mm9_reannotated_with_rr.bed'
halftimes_table = read.table(file = halftimes_file, sep = '\t', header = F)
colnames(halftimes_table) = c('chr','start','end','gene_name','halftime','strand')
halftimes_table$TSS[halftimes_table$strand == '+'] = halftimes_table$start[halftimes_table$strand == '+']
halftimes_table$TSS[halftimes_table$strand == '-'] = halftimes_table$end[halftimes_table$strand == '-']
halftimes_table = halftimes_table[order(halftimes_table$TSS),]

file_mm9_chrom_sizes = '/project/lncrna/Xist/data/annotation_files/mouse_genome/mm9.chrom.sizes'
mm9_chrom_sizes = read.table(file = file_mm9_chrom_sizes, header=F, sep='\t')
colnames(mm9_chrom_sizes) = c('chr','length')

initiation_sites_file = '/project/lncrna/Xist/data/annotation_files/xist_initiation_sites/1hr-Xist-Initiation-Sites_mm9.bed'
initiation_sites = read.table(initiation_sites_file)

CairoPDF(file = "/project/lncrna/Xist/plots/additional_analysis/analysis_halftime.pdf", width = 10, height = 10)

#plot distribution of halftimes
plot(density(halftimes_table$halftime),ylab = 'Density',xlab = 'Halftime [days]',main='PRO-Seq')
hist(halftimes_table$halftime,breaks = seq(0,3.6,by = 0.1),ylab = 'Frequency',xlab = 'Halftime [days]',main='PRO-Seq')

#plot halftime vs genomic position
plot(halftimes_table$TSS, halftimes_table$halftime, main='PRO-Seq', pch = '*', col = 'darkgray', xlab = 'Genomic Position (TSS of genes)', 
     ylab = 'half-time',xlim = c(0,mm9_chrom_sizes$length[mm9_chrom_sizes$chr == 'chrX']))
lines(halftimes_table$TSS,smooth(halftimes_table$halftime),col='black')
abline(v=100678572,col='darkorange') #TSS of Xist on mm9



######
#mrna-seq undiff
######
mrna_seq_undiff_file = "/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_mrna_seq_undiff_mm9_reannotated_with_rr.bed"
mrna_seq_undiff_table = read.table(file = mrna_seq_undiff_file, sep = '\t', header = F)
colnames(mrna_seq_undiff_table) = c('chr','start','end','gene_name','halftime','strand')
mrna_seq_undiff_table$TSS[mrna_seq_undiff_table$strand == '+'] = mrna_seq_undiff_table$start[mrna_seq_undiff_table$strand == '+']
mrna_seq_undiff_table$TSS[mrna_seq_undiff_table$strand == '-'] = mrna_seq_undiff_table$end[mrna_seq_undiff_table$strand == '-']
mrna_seq_undiff_table = mrna_seq_undiff_table[order(mrna_seq_undiff_table$TSS),]

#plot distribution of halftimes
plot(density(mrna_seq_undiff_table$halftime),ylab = 'Density',xlab = 'Halftime [days]',col='black',main='mRNA-Seq undifferentiated')
hist(mrna_seq_undiff_table$halftime,breaks = seq(0,3.6,by = 0.1),ylab = 'Frequency',xlab = 'Halftime [days]',main='mRNA-Seq undifferentiated')

#plot halftime vs genomic position
plot(mrna_seq_undiff_table$TSS, mrna_seq_undiff_table$halftime, main='mRNA-Seq undifferentiated', pch = '*', col = 'darkgray', xlab = 'Genomic Position (TSS of genes)', 
     ylab = 'half-time',xlim = c(0,mm9_chrom_sizes$length[mm9_chrom_sizes$chr == 'chrX']))
lines(mrna_seq_undiff_table$TSS,smooth(mrna_seq_undiff_table$halftime),col='black')
abline(v=100678572,col='darkorange') #TSS of Xist on mm9


######
#mrna-seq diff
######
mrna_seq_diff_file = "/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_mrna_seq_diff_mm9_reannotated_with_rr.bed"
mrna_seq_diff_table = read.table(file = mrna_seq_diff_file, sep = '\t', header = F)
colnames(mrna_seq_diff_table) = c('chr','start','end','gene_name','halftime','strand')
mrna_seq_diff_table$TSS[mrna_seq_diff_table$strand == '+'] = mrna_seq_diff_table$start[mrna_seq_diff_table$strand == '+']
mrna_seq_diff_table$TSS[mrna_seq_diff_table$strand == '-'] = mrna_seq_diff_table$end[mrna_seq_diff_table$strand == '-']
mrna_seq_diff_table = mrna_seq_diff_table[order(mrna_seq_diff_table$TSS),]

#plot distribution of halftimes
plot(density(mrna_seq_diff_table$halftime),ylab = 'Density',xlab = 'Halftime [days]',col='black',main='mRNA-Seq differentiated')
hist(mrna_seq_diff_table$halftime,breaks = seq(0,3.6,by = 0.1),ylab = 'Frequency',xlab = 'Halftime [days]',main='mRNA-Seq differentiated')

#plot halftime vs genomic position
plot(mrna_seq_diff_table$TSS, mrna_seq_diff_table$halftime, main='mRNA-Seq differentiated', pch = '*', col = 'darkgray', xlab = 'Genomic Position (TSS of genes)', 
     ylab = 'half-time',xlim = c(0,mm9_chrom_sizes$length[mm9_chrom_sizes$chr == 'chrX']))
lines(mrna_seq_diff_table$TSS,smooth(mrna_seq_diff_table$halftime),col='black')
abline(v=100678572,col='darkorange') #TSS of Xist on mm9

dev.off()
