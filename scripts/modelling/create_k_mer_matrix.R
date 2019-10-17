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

feature_matrix_name = 'promoter_matrix_reannotated_3_mer_pro_seq_genes_1000'
output_dir = 'data/modelling/feature_matrix'
tmp_dir = here(output_dir,'tmp')
system(paste('mkdir',tmp_dir))

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


##########gene regions from start to end and 1000 bp around the promoter
gene_regions = GRanges(seqnames = halftimes_table$chr, ranges = IRanges(halftimes_table$start,halftimes_table$end), strand = halftimes_table$strand)
mcols(gene_regions) = data.frame(name = halftimes_table$gene_name, score = halftimes_table$halftime)
#mcols(gene_regions) = data.frame(name = halftimes_table$gene_name, score = 0) #all genes

tmp_gene_regions_file = paste(tmp_dir,'/tmp_gene_regions_file.bed',sep='')
export.bed(gene_regions,tmp_gene_regions_file,format='bed')

gene_regions_1000 = promoters(gene_regions, upstream = 500, downstream = 500)
tmp_gene_regions_1000_file = paste(tmp_dir,'/tmp_gene_regions_1000_file.bed',sep='')
export.bed(gene_regions_1000,tmp_gene_regions_1000_file,format='bed')

###################################################################################
#get sequences
###################################################################################

mouse_genome_file = here('data/annotation_files/mouse_genome/mm9.fa')
tmp_sequence_file = paste(tmp_dir,'/tmp_sequences.txt',sep='')
cmd = paste(bedTools,'fastaFromBed -fi ',mouse_genome_file,' -bed ',tmp_gene_regions_1000_file,' -fo ',tmp_sequence_file, ' -s -name -tab',sep='')
system(cmd)

sequences = read.table(file=tmp_sequence_file,sep="\t")
sequence_matrix = t(sapply(strsplit(as.character(sequences[,2]),''),as.character))
rownames(sequence_matrix) = sequences[,1]
sequence_matrix = as.DNAbin(sequence_matrix)

kmers = kcount(sequence_matrix, k = 3,residues = 'DNA')

data_set = kmers
halftime = halftimes_table$halftime

save(data_set,halftime,file = here(output_dir,paste(feature_matrix_name,'.RData',sep='')))

##########remove tmp files
system(paste('rm -rf',tmp_dir))



















