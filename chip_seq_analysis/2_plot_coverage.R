library(Cairo)

samtools = '/home/lisasous/tools/samtools-1.9/samtools'
file_metadata = '/project/lncrna/Xist/data/chip_seq/metadata/metadata_merged_replicates.txt'
experiment_dir = '/project/ngs_marsico/Xist/bam/experiment/'
control_dir = '/project/ngs_marsico/Xist/bam/control/'
output_dir = '/project/lncrna/Xist/plots/chip_seq_analysis/library_coverage/'

metadata = read.table(file_metadata,sep='\t',header=T,comment.char='')
colnames(metadata)[1] = 'feature'
metadata$feature = as.character(metadata$feature)
metadata$accession_number = as.character(metadata$accession_number)
metadata$experiment_file = as.character(metadata$experiment_file)
metadata$control_file = as.character(metadata$control_file)

coverage_df = data.frame(title = unlist(lapply(strsplit(metadata$experiment_file,"[.]"),'[[',1)), coverage_experiment = 0, coverage_control = 0)

for (i in 1:nrow(metadata)) {
  feature = metadata$feature[i]
  print(feature)
  exp = paste(experiment_dir,metadata$experiment_file[i],sep='')
  ctr = paste(control_dir,metadata$control_file[i],sep='')

  
  cmd = paste(samtools, " view -F 0x904 -c ", exp, sep='')
  coverage_exp = as.numeric(system(cmd, intern = T))
  
  
  cmd = paste(samtools, " view -F 0x904 -c ", ctr, sep='')
  coverage_ctr = as.numeric(system(cmd, intern = T))
  
  coverage_df$coverage_experiment[i] = coverage_exp
  coverage_df$coverage_control[i] = coverage_ctr
}

data = coverage_df[,c(2,3)]
rownames(data) = coverage_df$title

file = paste(output_dir,'barplot_coverage.pdf',sep='')
CairoPDF(file = file, width = 20, height = 50)
par(mar = c(6,16,6,6),oma=c(4,4,4,4))
barplot(t(as.matrix(data)),main='library coverage',cex.main=2,xlab='coverage',beside=T,las=2,col=terrain.colors(2),horiz=T)
legend('topright', c("experiment","control"), cex=2, fill=terrain.colors(2))
abline(v=3000000,col='red')
dev.off()
