###################################################################################
#libraries
###################################################################################

library(normr)
library(GenomicRanges)
library(rtracklayer)
library(bamsignals)
library(here)

###################################################################################
#input directories
###################################################################################

file_genes = here('data/silencing_halftimes/fitted_data','halftimes_pro_seq_mm9_reannotated_with_rr.bed')
#file_genes = here('data/annotation_files/gene_annotation','gencode.vM9.annotation.chrX.genes.reannotated.with.rr.mm9.bed') #all genes on chrX
#file_genes = here('data/annotation_files/enhancers','gene_enhancers.bed') #enhancer regions
file_metadata = here('data/chip_seq/metadata','metadata_normalization_regions.txt')
file_mm9_chrom_sizes = here('data/annotation_files/mouse_genome','mm9.chrom.sizes')
ChIP_dir = '/project/ngs_marsico/Xist/bam/experiment/' #set path to directory with chip-seq bam files folder
ctrl_dir = '/project/ngs_marsico/Xist/bam/control/' #set path to directory with control bam files folder
output_dir = here('data/chip_seq/normalized_counts/')

###################################################################################
#load metadata
###################################################################################

metadata = read.table(file_metadata,sep='\t',header=T,comment.char='')
colnames(metadata)[1] = 'feature'
metadata$feature = as.character(metadata$feature)
metadata$accession_number = as.character(metadata$accession_number)
metadata$experiment_file = as.character(metadata$experiment_file)
metadata$control_file = as.character(metadata$control_file)
#for enhancer analysis
#metadata$start = NA
#metadata$end = NA

#chromosome
mm9_chrom_sizes = read.table(file = file_mm9_chrom_sizes, header=F, sep='\t')
colnames(mm9_chrom_sizes) = c('chr','length')

###################################################################################
#calculate normalized enrichment in given window
###################################################################################

for (i in 1:nrow(metadata)) {
  feature = metadata$feature[i]
  print(feature)
  acNr = metadata$accession_number[i]
  a = as.numeric(as.character(metadata$start[i]))
  b = as.numeric(as.character(metadata$end[i]))
  binsize_normR = as.numeric(as.character(metadata$binsize_normR[i]))
  title = paste(feature,acNr,a,b,'normalized_normr',sep='_')
  
  bed_file_out = paste(output_dir,title,'.bed',sep='')
  
  if(!file.exists(bed_file_out)){
  
    #gene regions
    table_genes = read.table(file = file_genes, sep = '\t', header = F)
    colnames(table_genes) = c('chr','start','end','gene_name','halftime','strand')
    if(is.na(a)){
      if(is.na(b)){
        gene_regions = GRanges(seqnames = table_genes$chr, ranges = IRanges(table_genes$start,table_genes$end), strand = table_genes$strand)
      }else{
        table_genes$gene_length = table_genes$end - table_genes$start
        table_genes$new_start = table_genes$start
        table_genes$new_end = table_genes$end
        table_genes$new_end[table_genes$strand == '+' & table_genes$gene_length > b] = table_genes$start[table_genes$strand == '+' & table_genes$gene_length > b] + b
        table_genes$new_start[table_genes$strand == '-' & table_genes$gene_length > b] = table_genes$end[table_genes$strand == '-' & table_genes$gene_length > b] - b
        gene_regions = GRanges(seqnames = table_genes$chr, ranges = IRanges(table_genes$new_start,table_genes$new_end), strand = table_genes$strand)
      }
    }else{
      table_genes$new_start[table_genes$strand == '+'] = table_genes$start[table_genes$strand == '+'] + a
      table_genes$new_end[table_genes$strand == '+'] = table_genes$start[table_genes$strand == '+'] + b
      table_genes$new_start[table_genes$strand == '-'] = table_genes$end[table_genes$strand == '-'] - b
      table_genes$new_end[table_genes$strand == '-'] = table_genes$end[table_genes$strand == '-'] - a
      gene_regions = GRanges(seqnames = table_genes$chr, ranges = IRanges(table_genes$new_start,table_genes$new_end), strand = table_genes$strand)
    }
    
    #count reads in gene regions
    bamfileChIP = paste(ChIP_dir,metadata$experiment_file[i],sep='')
    bamfileCtrl = paste(ctrl_dir,metadata$control_file[i],sep='')
    chipGeneCts = bamCount(bamfileChIP, gene_regions, mapqual = 10)
    ctrlGeneCts = bamCount(bamfileCtrl, gene_regions, mapqual = 10)
    
    #calculate normr normalization factors
    fit <- enrichR(bamfileChIP, bamfileCtrl, genome=mm9_chrom_sizes, countConfigSingleEnd(binsize = binsize_normR))
    ctsChIP <- getCounts(fit)$treatment
    ctsCtrl <- getCounts(fit)$control
    post <- getPosteriors(fit)[,1]
    
    #adjust binsize of fit to size of enrichment window
    binsize.fit = width(getRanges(fit))[2]
    if(is.na(a)){
      binsize.feature = table_genes$end - table_genes$start
      if(!is.na(b)){
        binsize.feature[binsize.feature > b] = b
      }
    }else{
      binsize.feature = abs(a-b) 
    }
    
    pseudoChIP <- sum(post * ctsChIP)/sum(post)  / (binsize.fit / binsize.feature)
    pseudoCtrl <- sum(post * ctsCtrl)/sum(post)  / (binsize.fit / binsize.feature)
    regul <- log(pseudoCtrl/pseudoChIP)
    stdrz <- log(fit@theta[2]/(1-fit@theta[2]) * (1-fit@theta[1])/fit@theta[1])
    
    #normalize read counts with normalization factors
    fc <- log((chipGeneCts + pseudoChIP)/(ctrlGeneCts + pseudoCtrl))
    normEnrichment <- (fc + regul)/stdrz
    normEnrichment[normEnrichment < 0] = 0
    normEnrichment[normEnrichment > 1] = 1
    
    #export normalized counts to bed
    normalized_counts = gene_regions
    mcols(normalized_counts) = data.frame(name = as.character(table_genes$gene_name), score = normEnrichment)
    export.bed(normalized_counts,bed_file_out,format='bed')
  }
}





