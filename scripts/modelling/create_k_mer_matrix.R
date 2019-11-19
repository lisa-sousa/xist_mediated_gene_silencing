###################################################################################
#libraries and tools
###################################################################################

library(here)
library(GenomicRanges)
library(rtracklayer)
library(ape)
library(kmer)
bedTools = '/home/lisasous/tools/bedtools2/bin/' #set path to bedtools 2 /bin/ folder

###################################################################################
#directories
###################################################################################

feature_matrix_name = 'promoter_pro_seq_genes_kmer1kb'
output_dir = 'data/modelling/feature_matrix'
tmp_dir = here(output_dir,'tmp')
system(paste('mkdir',tmp_dir))
x = 1000 #region around promoter

###################################################################################
#load gene annotation
###################################################################################

##########silencing halftimes
halftimes_file = here('data/silencing_halftimes/fitted_data','halftimes_pro_seq_mm9_reannotated_with_rr.bed')
#halftimes_file = here('data/annotation_files/gene_annotation','gencode.vM9.annotation.chrX.genes.reannotated.with.rr.mm9.bed') #all genes
halftimes_table = read.table(file = halftimes_file, sep = '\t', header = F)
colnames(halftimes_table) = c('chr','start','end','gene_name','halftime','strand')

halftimes_table$TSS[halftimes_table$strand == '+'] = halftimes_table$start[halftimes_table$strand == '+']
halftimes_table$TSS[halftimes_table$strand == '-'] = halftimes_table$end[halftimes_table$strand == '-']


##########gene regions from start to end and x bp around the promoter
gene_regions = GRanges(seqnames = halftimes_table$chr, ranges = IRanges(halftimes_table$start,halftimes_table$end), strand = halftimes_table$strand)
mcols(gene_regions) = data.frame(name = halftimes_table$gene_name, score = halftimes_table$halftime)
#mcols(gene_regions) = data.frame(name = halftimes_table$gene_name, score = 0) #all genes

tmp_gene_regions_file = paste(tmp_dir,'/tmp_gene_regions_file.bed',sep='')
export.bed(gene_regions,tmp_gene_regions_file,format='bed')

gene_regions_x = promoters(gene_regions, upstream = x/2, downstream = x/2)
tmp_gene_regions_x_file = paste(tmp_dir,'/tmp_gene_regions_x_file.bed',sep='')
export.bed(gene_regions_x,tmp_gene_regions_x_file,format='bed')

###################################################################################
#create matrix
###################################################################################

##########get kmers
mouse_genome_file = here('data/annotation_files/mouse_genome/mm9.fa')
tmp_sequence_file = paste(tmp_dir,'/tmp_sequences.txt',sep='')
cmd = paste(bedTools,'fastaFromBed -fi ',mouse_genome_file,' -bed ',tmp_gene_regions_x_file,' -fo ',tmp_sequence_file, ' -s -name -tab',sep='')
system(cmd)

sequences = read.table(file=tmp_sequence_file,sep="\t")
sequence_matrix = t(sapply(strsplit(as.character(sequences[,2]),''),as.character))
rownames(sequence_matrix) = sequences[,1]
sequence_matrix = as.DNAbin(sequence_matrix)

kmers_k3 = kcount(sequence_matrix, k = 3,residues = 'DNA')
kmers_k5 = kcount(sequence_matrix, k = 5,residues = 'DNA')

##########LINE density around each gene
LINEs_file = here('data/annotation_files/LINEs','LINEs_mm9.bed')
LINEs_table = read.table(LINEs_file,sep='\t')
colnames(LINEs_table) = c("chr","start","end","name","score","strand","dispStart","dispEnd","col")
LINEs_table = LINEs_table[LINEs_table$chr == 'chrX',]

tmp_LINEs_file = paste(tmp_dir,'/tmp_LINEs_file.bed',sep='')
export.bed(LINEs_table,tmp_LINEs_file,format='bed')

tmp_LINE_density_file = paste(tmp_dir,'/tmp_LINE_density.txt',sep='')
cmd = paste(bedTools,'windowBed -a ',tmp_gene_regions_file,' -b ',tmp_LINEs_file,' -w ',(x/2),' -c > ',tmp_LINE_density_file,sep='') #-w 	Base pairs added upstream and downstream of each entry in A when searching for overlaps in B. Default is 1000 bp.
system(cmd)
LINE_density_table = read.table(tmp_LINE_density_file,header=F,sep='\t')[,c(4,7)]
colnames(LINE_density_table) = c('gene_name','LINE_density')

#gene body
tmp_LINE_density_genebody_file = paste(tmp_dir,'/tmp_LINE_density_genebody.txt',sep='')
cmd = paste(bedTools,'intersectBed -a ',tmp_gene_regions_file,' -b ',tmp_LINEs_file,' -c > ',tmp_LINE_density_genebody_file,sep='') 
system(cmd)
LINE_density_genebody_table = read.table(tmp_LINE_density_genebody_file,header=F,sep='\t')[,c(4,7)]
colnames(LINE_density_genebody_table) = c('gene_name','LINE_density_genebody')


##########ALU density around each gene
ALUs_file = here('data/annotation_files/ALUs','ALUs_mm9.bed')

tmp_ALUs_density_file = paste(tmp_dir,'/tmp_ALUsE_density.txt',sep='')
cmd = paste(bedTools,'windowBed -a ',tmp_gene_regions_file,' -b ',ALUs_file,' -w ',(x/2),' -c > ',tmp_ALUs_density_file,sep='') #-w 	Base pairs added upstream and downstream of each entry in A when searching for overlaps in B. Default is 1000 bp.
system(cmd)
ALUs_density_table = read.table(tmp_ALUs_density_file,header=F,sep='\t')[,c(4,7)]
colnames(ALUs_density_table) = c('gene_name','ALUs_density')

#gene body
tmp_ALUs_density_genebody_file = paste(tmp_dir,'/tmp_ALUs_density_genebody.txt',sep='')
cmd = paste(bedTools,'intersectBed -a ',tmp_gene_regions_file,' -b ',ALUs_file,' -c > ',tmp_ALUs_density_genebody_file,sep='') 
system(cmd)
ALUs_density_genebody_table = read.table(tmp_ALUs_density_genebody_file,header=F,sep='\t')[,c(4,7)]
colnames(ALUs_density_genebody_table) = c('gene_name','ALUs_density_genebody')

##########save data matrix
data_set = cbind(kmers_k3,kmers_k5,LINE_density_table[,2],LINE_density_genebody_table[,2],ALUs_density_table[,2],ALUs_density_genebody_table[,2])
halftime = halftimes_table$halftime
save(data_set,halftime,file = here(output_dir,paste(feature_matrix_name,'.RData',sep='')))


##########remove tmp files
system(paste('rm -rf',tmp_dir))



















