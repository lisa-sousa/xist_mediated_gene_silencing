###################################################################################
#input directories
###################################################################################

library(here)
deepTools = '/home/lisasous/tools/deepTools2.0/bin/' #set path to deepTools 2.0 /bin/ folder
file_metadata = here('data/chip_seq/metadata','metadata_filtered_by_coverage.txt')
experiment_dir = '/project/ngs_marsico/Xist/bam/experiment/' #set path to directory with chip-seq bam files folder
control_dir = '/project/ngs_marsico/Xist/bam/control/' #set path to directory with control bam files folder
output_dir_data = here('data/chip_seq/bigwigs/raw/')
output_dir_plot = here('plots/chip_seq_analysis/heatmap/data_raw/')
gene_region = here('data/silencing_halftimes/fitted_data','halftimes_pro_seq_mm9_new_gene_annotation.bed')

###################################################################################
#parameters
###################################################################################

a = 2000 #bp downstream of TSS
b = 2000 #bp upstream of TSS
binsize.bigwig = 25
binsize.heatmap = 100

###################################################################################
#load metadata
###################################################################################

metadata = read.table(file_metadata,sep='\t',header=T,comment.char='')
colnames(metadata)[1] = 'feature'
metadata$feature = as.character(metadata$feature)
metadata$accession_number = as.character(metadata$accession_number)
metadata$experiment_file = as.character(metadata$experiment_file)
metadata$control_file = as.character(metadata$control_file)
titles = unlist(lapply(strsplit(metadata$experiment_file,"[.]"),'[[',1))

#create temporary directory
cmd = paste("mkdir ",output_dir_data,"tmp",sep="")
print(cmd)
system(cmd)

###################################################################################
#create deepTools plots
###################################################################################

for (i in 1:nrow(metadata)) {
  print(metadata$feature[i])
  title = titles[i]
  exp = paste(experiment_dir,metadata$experiment_file[i],sep='')
  exp_bw = paste(output_dir_data,unlist(strsplit(metadata$experiment_file[i],split='\\.'))[1],'.bw',sep='')
  ctr = paste(control_dir,metadata$control_file[i],sep='')
  ctr_bw = paste(output_dir_data,unlist(strsplit(metadata$control_file[i],split='\\.'))[1],'.bw',sep='')

  if(!file.exists(exp_bw)){
    cmd = paste(deepTools,'bamCoverage -b ',exp,' -o ',exp_bw,' -of bigwig --binSize ', binsize.bigwig, ' --numberOfProcessors max/2',sep='')
    #cmd = paste(deepTools,'bamCoverage -b ',exp,' -o ',exp_bw,' -of bigwig --binSize ', binsize.bigwig, ' --normalizeUsing RPKM --numberOfProcessors max/2',sep='')
    print(cmd)
    system(cmd)
  }
  if(!file.exists(ctr_bw)){
    cmd = paste(deepTools,'bamCoverage -b ',ctr,' -o ',ctr_bw,' -of bigwig --binSize ', binsize.bigwig, ' --numberOfProcessors max/2',sep='')
    #cmd = paste(deepTools,'bamCoverage -b ',ctr,' -o ',ctr_bw,' -of bigwig --binSize ', binsize.bigwig, ' --normalizeUsing RPKM --numberOfProcessors max/2',sep='')
    print(cmd)
    system(cmd)
  }  
  
  
  #plot a deepTools heatmap for given regions
  mat_file = paste(output_dir_data,'tmp/',title,'.mat.gz',sep='')
  cmd = paste(deepTools,'computeMatrix reference-point -S ',ctr_bw,' ',exp_bw,' -R ',gene_region,' --referencePoint TSS -a ',a,' -b ',b,
              ' --binSize ', binsize.heatmap, ' --numberOfProcessors max/2 --outFileName ',mat_file,sep='')
  print(cmd)
  system(cmd)
  
  
  plotFile = paste(output_dir_plot,title,'.pdf',sep='')
  cmd = paste(deepTools,'plotHeatmap -m ',mat_file,' -out ',plotFile,
                ' --averageTypeSummaryPlot mean --sortRegions ascend --colorList white,darksalmon,brown --colorNumber 20',sep='')
  print(cmd)
  system(cmd)
}

cmd = paste('pdfunite ',output_dir_plot,'*.pdf ',output_dir_plot,'heatmaps.pdf',sep='')
system(cmd)

cmd = paste("rm -rf ",output_dir_data,"tmp",sep="")
print(cmd)
system(cmd)
