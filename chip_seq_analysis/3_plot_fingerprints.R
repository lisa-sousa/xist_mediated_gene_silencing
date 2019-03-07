###################################################################################
#input directories
###################################################################################

samtools = '/home/lisasous/tools/samtools-1.9/samtools'
plotFingerprint = '/home/lisasous/tools/deepTools2.0/bin/plotFingerprint'
file_metadata = '/project/lncrna/Xist/data/chip_seq/metadata/metadata_filtered_by_coverage.txt'
experiment_dir = '/project/ngs_marsico/Xist/bam/experiment/'
control_dir = '/project/ngs_marsico/Xist/bam/control/'
output_dir = '/project/lncrna/Xist/plots/chip_seq_analysis/fingerprints/'

###################################################################################
#load metadata
###################################################################################

metadata = read.table(file_metadata,sep='\t',header=T,comment.char='')
colnames(metadata)[1] = 'feature'
metadata$feature = as.character(metadata$feature)
metadata$accession_number = as.character(metadata$accession_number)
metadata$experiment_file = as.character(metadata$experiment_file)
metadata$control_file = as.character(metadata$control_file)
plotTitles = unlist(lapply(strsplit(metadata$experiment_file,"[.]"),'[[',1))

###################################################################################
#plot fingerprints
###################################################################################

for (i in 1:nrow(metadata)) {
  feature = metadata$feature[i]
  print(feature)
  acNr = metadata$accession_number[i]
  exp = paste(experiment_dir,metadata$experiment_file[i],sep='')
  exp_index = paste(exp,'.bai',sep='')
  ctr = paste(control_dir,metadata$control_file[i],sep='')
  ctr_index = paste(ctr,'.bai',sep='')
  if(!file.exists(exp_index)){
    cmd = paste(samtools,'index',exp)
    system(cmd)
  }
  if(!file.exists(ctr_index)){
    cmd = paste(samtools,'index',ctr)
    system(cmd)
  }
  
  plotTitle = plotTitles[i]
  plotFile = paste(output_dir,plotTitle,'.pdf',sep='')
  if(!file.exists(plotFile)){
	  cmd = paste(plotFingerprint,'-b',exp,ctr,'--numberOfProcessors max/2','--plotTitle',plotTitle,'--labels',feature,'input','--plotFile',plotFile)
    #cmd = paste(plotFingerprint,'-b',exp,ctr,'-r chrX','--numberOfProcessors max/2','--plotTitle',plotTitle,'--labels',feature,'input','--plotFile',plotFile)
    system(cmd)
  }

}

cmd = paste('pdfunite ',output_dir,'*.pdf ',output_dir,'fingerprints.pdf',sep='')
system(cmd)
