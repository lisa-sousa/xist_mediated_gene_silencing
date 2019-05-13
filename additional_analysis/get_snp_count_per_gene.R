####################################################
#libraries and tools
####################################################

library(rtracklayer)
library(GenomicRanges)

bedTools = '/home/lisasous/tools/bedtools2/bin/'

####################################################
#directories and files
####################################################

output_dir = "/project/lncrna/Xist/data/annotation_files/SNPs/"

#input files
file_mm10_chrom_sizes = '/project/lncrna/Xist/data/annotation_files/mouse_genome/mm10.chrom.sizes'
file_SNPs_raw = "/project/lncrna/Xist/data/annotation_files/SNPs/all_TX1072_cellLine_B6_CAST.txt"
file_gencode = "/project/lncrna/Xist/data/annotation_files/gene_annotation/gencode.vM9.annotation.gtf"

#processing files
file_SNPs_bed = paste(output_dir,rev(unlist(strsplit(unlist(strsplit(file_SNPs_raw,'[.]'))[1],"/")))[1],'.bed',sep='')
file_exons_bed = paste(output_dir,rev(unlist(strsplit(unlist(strsplit(file_gencode,'.gtf'))[1],"/")))[1],'.exons.bed',sep='')
file_overlap_exons_SNPs = paste(output_dir,rev(unlist(strsplit(unlist(strsplit(file_SNPs_raw,'[.]'))[1],"/")))[1],'.overlap',sep='')
file_genes_SNP_count = paste(output_dir,rev(unlist(strsplit(unlist(strsplit(file_gencode,'.gtf'))[1],"/")))[1],'.SNP.count.bed',sep='')

####################################################
#functions
####################################################

get_attributes <- function(gencode_table, info_gene_name, info_gene_type){
  gencode_table$V9 = as.character(gencode_table$V9)
  gencode_table$gene_name = ''
  
  for(i in 1:nrow(gencode_table)){
    info = unlist(strsplit(gencode_table$V9[i],split = "; "))
    gencode_table$gene_name[i] = unlist(strsplit(info[info_gene_name],split = " "))[2]
  }
  gencode_table$V9 = NULL
  return(gencode_table)
}

####################################################
#get number of SNPs per gene
####################################################

#create mm10 sequence info
mm10_chrom_sizes = read.table(file = file_mm10_chrom_sizes, header=F, sep='\t')
mm10_seq_info = Seqinfo(seqnames = as.character(mm10_chrom_sizes[,1]), seqlengths = mm10_chrom_sizes[,2], isCircular = rep(FALSE,nrow(mm10_chrom_sizes)), genome = 'mm10')

#convert SNP caller file type to bed file
table_SNPs = read.table(file_SNPs_raw,header=F,sep='\t')
bed_format = sort(GRanges(seqnames = table_SNPs$V2, ranges = IRanges(start = table_SNPs$V3, end = table_SNPs$V3),seqinfo = mm10_seq_info))
export.bed(bed_format,file_SNPs_bed,format='bed')

#create exon list
annotation_gencode = read.table(file=file_gencode, sep = "\t")
annotation_gencode = annotation_gencode[annotation_gencode$V1 == "chrX",]

exon_annotation_gencode = annotation_gencode[annotation_gencode$V3 == "exon",]
exon_annotation_gencode = get_attributes(exon_annotation_gencode, 5, 3)
colnames(exon_annotation_gencode) = c("seqname","source","feature","start","end","score","strand","frame","attribute")
granges_exons = makeGRangesFromDataFrame(exon_annotation_gencode)
mcols(granges_exons) = data.frame(name = exon_annotation_gencode$attribute)
export.bed(granges_exons,file_exons_bed,format='bed')

#overlap exons with SNPs
cmd = paste(bedTools, "intersectBed -a ", file_exons_bed, " -b ", file_SNPs_bed, " -c > ", file_overlap_exons_SNPs, sep='')
system(cmd)

table_overlap = read.table(file = file_overlap_exons_SNPs, sep = '\t', header = F, comment.char = '#')
table_overlap = table_overlap[table_overlap$V7 != 0,]

#create SNPs list
gene_annotation_gencode = annotation_gencode[annotation_gencode$V3 == "gene",]
gene_annotation_gencode = get_attributes(gene_annotation_gencode, 4, 2)
gene_annotation_gencode = gene_annotation_gencode[,c(1,4,5,9,7)]
gene_annotation_gencode$SNP_count = 0

for(i in 1:nrow(gene_annotation_gencode)){
  gene_annotation_gencode$SNP_count[i] = sum(table_overlap$V7[table_overlap$V4 == gene_annotation_gencode$gene_name[i]])
}
gene_annotation_gencode = gene_annotation_gencode[order(gene_annotation_gencode$gene_name),]

write.table(gene_annotation_gencode[,c(1,2,3,4,6,5)], file = file_genes_SNP_count, sep = '\t', col.names = F, row.names = F, quote = F)
