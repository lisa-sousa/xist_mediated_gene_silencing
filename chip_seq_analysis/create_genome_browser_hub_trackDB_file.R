file_metadata = '/project/lncrna/Xist/data_lisa/chip_seq/metadata/metadata_filtered_by_deepTools.txt'
output_dir = '/home/lisasous/'

metadata = read.table(file_metadata,sep='\t',header=T,comment.char='')
colnames(metadata)[1] = 'feature'
metadata$feature = as.character(metadata$feature)
metadata$accession_number = as.character(metadata$accession_number)
metadata$experiment_file = as.character(metadata$experiment_file)
metadata$control_file = as.character(metadata$control_file)
titles = unlist(lapply(strsplit(metadata$experiment_file,"[.]"),'[[',1))
metadata$bigwig = paste(lapply(strsplit(metadata$experiment_file,split='\\.'), `[[`, 1),'.bw',sep='')

trackDB_file = file(paste(output_dir,'trackDb.txt',sep=''), open = "a")

for(i in 1:nrow(metadata)){
  feature = paste(metadata$feature[i],"_",metadata$accession_number[i],sep="")
  bigwig_file = metadata$bigwig[i]
  cat(paste("track chip_seq_",feature,sep=""), file = trackDB_file, sep = "\n")
  cat(paste("bigDataUrl",bigwig_file,sep=" "), file = trackDB_file, sep = "\n")
  cat(paste("shortLabel",feature,sep=" "), file = trackDB_file, sep = "\n")
  cat(paste("longLabel raw chip-seq data of",feature,sep=" "), file = trackDB_file, sep = "\n")
  cat("type bigWig", file = trackDB_file, sep = "\n")
  cat("color 0,0,400", file = trackDB_file, sep = "\n")
  cat("visibility full", file = trackDB_file, sep = "\n")
  cat("autoScale off", file = trackDB_file, sep = "\n")
  cat("graphType bar", file = trackDB_file, sep = "\n")
  cat("smoothingWindow 10", file = trackDB_file, sep = "\n")
  cat("windowingFunction mean", file = trackDB_file, sep = "\n")
  cat("viewLimits 0:10", file = trackDB_file, sep = "\n")
  cat("", file = trackDB_file, sep = "\n")
}

close(trackDB_file)

