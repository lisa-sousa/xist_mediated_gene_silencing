
library(energy)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)

##############################
#READ AND MERGE DATA
##############################

region_name = 'promoter'
bedTools = '/home/lisasous/tools/bedtools2/bin/'
output_dir = '/project/lncrna/Xist/data/modelling/feature_matrix/'
tmp_dir = paste(output_dir,'tmp',sep='')
cmd = paste('mkdir',tmp_dir)
system(cmd)


##########genes mm9 from gencode
genes_mm9_gencode_file = '/project/lncrna/Xist/data/annotation_files/gene_annotation/gencode.vM9.annotation.chrX.genes.mm9.bed'
genes_mm9_gencode_table = read.table(genes_mm9_gencode_file,sep = '\t',header=F)
colnames(genes_mm9_gencode_table) = c('chr','start','end','gene_name','score','strand')


##########silencing halftimes
halftimes_file = '/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_pro_seq_mm9_reannotated_with_rr.bed'
#halftimes_file = '/project/lncrna/Xist/data_lisa/annotation_files/gene_annotation/gencode.vM9.annotation.chrX.genes.reannotated.with.rr.mm9.bed' #all genes
halftimes_table = read.table(file = halftimes_file, sep = '\t', header = F)
colnames(halftimes_table) = c('chr','start','end','gene_name','halftime','strand')

halftimes_table$TSS[halftimes_table$strand == '+'] = halftimes_table$start[halftimes_table$strand == '+']
halftimes_table$TSS[halftimes_table$strand == '-'] = halftimes_table$end[halftimes_table$strand == '-']


##########gene regions
gene_regions = GRanges(seqnames = halftimes_table$chr, ranges = IRanges(halftimes_table$start,halftimes_table$end), strand = halftimes_table$strand)
mcols(gene_regions) = data.frame(name = halftimes_table$gene_name, score = halftimes_table$halftime)
#mcols(gene_regions) = data.frame(name = halftimes_table$gene_name, score = 0) #all genes

tmp_gene_regions_file = paste(tmp_dir,'/tmp_gene_regions_file.bed',sep='')
export.bed(gene_regions,tmp_gene_regions_file,format='bed')

gene_regions_1000 = promoters(gene_regions, upstream = 500, downstream = 500)
tmp_gene_regions_1000_file = paste(tmp_dir,'/tmp_gene_regions_1000_file.bed',sep='')
export.bed(gene_regions_1000,tmp_gene_regions_1000_file,format='bed')


##########ChIP-Seq data
#metadata
input_dir = '/project/lncrna/Xist/data/chip_seq/normalized_counts/reannotated_normRAdjusted_pro_seq_genes/'
#input_dir = '/project/lncrna/Xist/data/chip_seq/normalized_counts/reannotated_normRAdjusted_all_chrX_genes/' #all genes
file_metadata = '/project/lncrna/Xist/data/chip_seq/metadata/metadata_normalization_regions.txt'
metadata = read.table(file_metadata,sep='\t',header=T,comment.char='')
colnames(metadata)[1] = 'feature'
metadata$feature = as.character(metadata$feature)
metadata$accession_number = as.character(metadata$accession_number)
metadata$experiment_file = as.character(metadata$experiment_file)
metadata$control_file = as.character(metadata$control_file)

chip_seq_table = data.frame(gene_name = halftimes_table$gene_name)

for (i in 1:nrow(metadata)) {
  feature = metadata$feature[i]
  print(feature)
  acNr = metadata$accession_number[i]
  a = metadata$start[i]
  b = metadata$end[i]
  title = paste(feature,acNr,a,b,sep='_')
  bed_file = paste(input_dir,title,'_normalized_normr.bed',sep='')
  
  chip_seq = read.table(bed_file,sep='\t')
  chip_seq = data.frame(gene_name = chip_seq$V4, score = chip_seq$V5)
  colnames(chip_seq) = c('gene_name',title)
  
  chip_seq_table = merge(chip_seq_table,chip_seq,by='gene_name')
}


##########distance to closest topological domain (TAD) (mapped to mm9)
TADs_file = '/project/lncrna/Xist/data/annotation_files/TADs/HindIII_combined/topological_domains_chrX.domains'
TADs = read.table(TADs_file,sep = '\t',header=F)
TADs = c(TADs$V2,TADs$V3)
all_distances <- abs(outer(halftimes_table$TSS,TADs,FUN="-"))
distance_TAD_table = data.frame(gene_name = halftimes_table$gene_name, distance_to_TAD_border = apply(all_distances,1,min))


##########gene density in 1mb around each gene
tmp_gene_density_file = paste(tmp_dir,'/tmp_gene_density.txt',sep='')
cmd = paste(bedTools,'windowBed -a ',tmp_gene_regions_file,' -b ',genes_mm9_gencode_file,' -w 500000 -c > ',tmp_gene_density_file,sep='') #-w 	Base pairs added upstream and downstream of each entry in A when searching for overlaps in B. Default is 1000 bp.
system(cmd)
gene_density_gencode_table = read.table(tmp_gene_density_file,header=F,sep='\t')[,c(4,7)]
colnames(gene_density_gencode_table) = c('gene_name','gene_density')


##########overlap with Xist initiation site (mapped to mm9)
initiation_sites_file = '/project/lncrna/Xist/data/annotation_files/xist_initiation_sites/1hr-Xist-Initiation-Sites_mm9.bed'
tmp_overlap_initiation_file = paste(tmp_dir,'/tmp_initiation_sites.txt',sep='')
cmd = paste(bedTools,'intersectBed -a ',tmp_gene_regions_file,' -b ',initiation_sites_file,' -c > ',tmp_overlap_initiation_file,sep='')
system(cmd)
overlap_initiation_table = read.table(tmp_overlap_initiation_file, sep='\t', header=F)
overlap_initiation_table = data.frame(gene_name = overlap_initiation_table$V4, overlap_with_Xist_entry_sites = overlap_initiation_table$V7)
overlap_initiation_table$overlap_with_Xist_entry_sites[overlap_initiation_table$overlap_with_Xist_entry_sites > 0] = 1


##########distance to the Xist locus
TSS_Xist = genes_mm9_gencode_table$end[genes_mm9_gencode_table$gene_name == 'Xist'] #Xist is on the minus strand
distance_xist_table = data.frame(gene_name = halftimes_table$gene_name, distance_to_Xist = abs(halftimes_table$TSS - TSS_Xist))


##########distance (closest) and overlap to Lamina Associated Domain (LADs)
LADs_file = '/project/lncrna/Xist/data/annotation_files/LADs/lads.csv'
LADs_table = read.table(LADs_file,sep='\t',header=T)
LADs_table = LADs_table[LADs_table$seqname == 'chrX',]
tmp_LADs_file = paste(tmp_dir,'/tmp_LADs_file.bed',sep='')
export.bed(LADs_table,tmp_LADs_file,format='bed')

tmp_overlap_LADs_file = paste(tmp_dir,'/tmp_overlap_LADs_file.bed',sep='')
cmd = paste(bedTools,'intersectBed -a ',tmp_gene_regions_1000_file,' -b ',tmp_LADs_file,' -c > ',tmp_overlap_LADs_file,sep='')
system(cmd)
overlap_LADs_table = read.table(tmp_overlap_LADs_file, sep='\t', header=F)
overlap_LADs_table = data.frame(gene_name = overlap_LADs_table$V4, overlap_with_LADs = overlap_LADs_table$V7)

LADs = c(LADs_table$start,LADs_table$end)
all_distances <- abs(outer(halftimes_table$TSS,LADs,FUN="-"))
distance_LAD_table = data.frame(gene_name = halftimes_table$gene_name, distance_to_LADs = apply(all_distances,1,min))


##########Hi-C interaction between gene promoter and other regulatory regions (mapped to mm9)
hic_promoter_other_file = '/project/lncrna/Xist/data/annotation_files/Hi-C_data/ESC_promoter_other_significant_interactions_HIC.txt'
hic_promoter_other_table = read.table(hic_promoter_other_file,sep='\t',header=T)

hic_promoter_promoter_file = '/project/lncrna/Xist/data/annotation_files/Hi-C_data/ESC_promoter_promoter_significant_interactions_HIC.txt'
hic_promoter_promoter_table = read.table(hic_promoter_promoter_file,sep='\t',header=T)

hic_table = data.frame(gene_name = halftimes_table$gene_name, number_interactions_HiC_all = 0, mean_interaction_strength_HiC_all = 0, number_interactions_HiC_promoter = 0, 
                       mean_interaction_strength_HiC_promoter = 0, mean_interaction_strength_HiC_xist = 0)
for(i in 1:nrow(hic_table)){
  interaction_other = hic_promoter_other_table[grepl(hic_table$gene_name[i],hic_promoter_other_table$Symbol),]
  interaction_promoter = hic_promoter_promoter_table[(grepl(hic_table$gene_name[i],hic_promoter_promoter_table$Symbol.1)) | (grepl(hic_table$gene_name[i],hic_promoter_promoter_table$Symbol)),]
  hic_table$number_interactions_HiC_all[i] = nrow(interaction_other)
  hic_table$mean_interaction_strength_HiC_all[i] = sum(interaction_other$log.observed.expected.)/nrow(interaction_other)
  hic_table$number_interactions_HiC_promoter[i] = nrow(interaction_promoter)
  hic_table$mean_interaction_strength_HiC_promoter[i] = sum(interaction_promoter$log.observed.expected.)/nrow(interaction_promoter)
  hic_table$mean_interaction_strength_HiC_xist[i] = sum(interaction_other$log.observed.expected.[(interaction_other$chr == 'chrX') & (interaction_other$start == '100647076')],
                                          interaction_promoter$log.observed.expected.[(grepl(as.factor('Xist'),interaction_promoter$Symbol.1)) | (grepl(as.factor('Xist'),interaction_promoter$Symbol))])
}
hic_table[is.na(hic_table)] = 0

##########HICap interaction between gene promoter and other regulatory regions (mapped to mm9)
hicap_promoter_enhancer_file = '/project/lncrna/Xist/data/annotation_files/HICap_data/Promoter_Enhancer_Interactions.csv'
hicap_promoter_enhancer_table = read.table(hicap_promoter_enhancer_file,sep='\t',header=T)

hicap_promoter_promoter_file = '/project/lncrna/Xist/data/annotation_files/HICap_data/Promoter_Promoter_Interactions.csv'
hicap_promoter_promoter_table = read.table(hicap_promoter_promoter_file,sep='\t',header=T)


hicap_table = data.frame(gene_name = halftimes_table$gene_name, number_interactions_HiCap_promoter = 0, number_interactions_HiCap_enhancer = 0, number_interactions_HiCap_all = 0)
for(i in 1:nrow(hicap_table)){
  interactions_promoter = hicap_promoter_promoter_table[as.character(hicap_promoter_promoter_table$gene_name) == as.character(hicap_table$gene_name[i]),]
  interactions_enhancer = hicap_promoter_enhancer_table[as.character(hicap_promoter_enhancer_table$gene_name) == as.character(hicap_table$gene_name[i]),]
  hicap_table$number_interactions_HiCap_promoter[i] = sum(interactions_promoter$read_pairs_in_HiCap_replicate1,interactions_promoter$read_pairs_in_HiCap_replicate2)/(2*nrow(interactions_promoter)) #mean connectivity for two replicates
  hicap_table$number_interactions_HiCap_enhancer[i] = sum(interactions_enhancer$read_pairs_in_HiCap_replicate1,interactions_enhancer$read_pairs_in_HiCap_replicate2)/(2*nrow(interactions_enhancer)) #mean connectivity for two replicates
  hicap_table[is.na(hicap_table)] = 0
  hicap_table$number_interactions_HiCap_all[i] = hicap_table$number_interactions_HiCap_promoter[i] + hicap_table$number_interactions_HiCap_enhancer[i]
}


##########DNA methylation 
dna_methylation_file = '/project/lncrna/Xist/data/annotation_files/Bisulfite_Seq_data/GSE30202_BisSeq_ES_CpGmeth.tsv'
dna_methylation_table = read.table(dna_methylation_file,sep='\t',header=T)
dna_methylation_table = dna_methylation_table[dna_methylation_table$chr == 'chrX',]
dna_methylation_table = data.frame(chr = dna_methylation_table$chr, start = dna_methylation_table$position, end = dna_methylation_table$position, name = 'meth_base',
                                   score = dna_methylation_table$nMeth/dna_methylation_table$nTot, strand = '.')

tmp_dna_methylation_file = paste(tmp_dir,'/tmp_dna_methylation_file.bed',sep='')
export.bed(dna_methylation_table,tmp_dna_methylation_file,format='bed')

tmp_overlap_dna_methylation_file = paste(tmp_dir,'/tmp_overlap_dna_methylation_file.bed',sep='')
cmd = paste(bedTools,'intersectBed -a ',tmp_gene_regions_1000_file,' -b ',tmp_dna_methylation_file,' -wo > ',tmp_overlap_dna_methylation_file,sep='')
system(cmd)

overlap_dna_methylation_table = read.table(tmp_overlap_dna_methylation_file,sep='\t')
overlap_dna_methylation_table$V4 = as.character(overlap_dna_methylation_table$V4)

dna_methylation_table = data.frame(gene_name = as.character(halftimes_table$gene_name), DNA_methylation_BS_Seq = 0)
for(i in 1:nrow(dna_methylation_table)){
  dna_methylation = overlap_dna_methylation_table[overlap_dna_methylation_table$V4 == dna_methylation_table$gene_name[i],]
  dna_methylation_table$DNA_methylation_BS_Seq[i] = mean(dna_methylation$V11)
}
dna_methylation_table$DNA_methylation_BS_Seq[is.na(dna_methylation_table$DNA_methylation_BS_Seq)] = 0

##########CpG islands
CpG_island_file = '/project/lncrna/Xist/data/annotation_files/CpG_islands/CGI_mm9_ucsc'
CpG_island_table = read.table(CpG_island_file,sep='\t')
CpG_island_table = data.frame(chr = CpG_island_table$V2, start = CpG_island_table$V3, end = CpG_island_table$V4, name = 'CpG_island', score = 0, strand = '+')

tmp_CpG_island_file = paste(tmp_dir,'/tmp_CpG_island_file.bed',sep='')
export.bed(CpG_island_table,tmp_CpG_island_file,format='bed')

tmp_overlap_CpG_island_file = paste(tmp_dir,'/tmp_overlap_CpG_island_file.bed',sep='')
cmd = paste(bedTools,'intersectBed -a ',tmp_gene_regions_1000_file,' -b ',tmp_CpG_island_file,' -c > ',tmp_overlap_CpG_island_file,sep='')
system(cmd)
overlap_CpG_island_table = read.table(tmp_overlap_CpG_island_file, sep='\t', header=F)
CpG_island_table = data.frame(gene_name = overlap_CpG_island_table$V4, overlap_with_CpG_islands = overlap_CpG_island_table$V7)


##########CpG content
mm9_fasta = '/project/lncrna/Xist/data/annotation_files/mouse_genome/mm9.fa'

tmp_gene_regions_1000_fasta = paste(tmp_dir,'/tmp_gene_regions_1000_file.fa',sep='')
cmd = paste(bedTools,'fastaFromBed -fi ',mm9_fasta,' -bed ',tmp_gene_regions_1000_file,' -fo ',tmp_gene_regions_1000_fasta,sep='')
system(cmd)

sequences = readDNAStringSet(tmp_gene_regions_1000_fasta)
dinucleotide_frequency = dinucleotideFrequency(sequences,step = 1)
frequencies = data.frame(CpG = dinucleotide_frequency[,colnames(dinucleotide_frequency) == 'CG'], G = letterFrequency(sequences,letters = 'G'), C = letterFrequency(sequences,letters = 'C'))
L = width(sequences)
CpG_content_table = data.frame(gene_name = halftimes_table$gene_name, CpG_content = (frequencies$CpG/L) / (((frequencies$G+frequencies$C)/(2*L))^2)) #normalization of CpG content taken from: http://nar.oxfordjournals.org/content/37/19/6305.full.pdf+html


##########distance (closest) to LINE and LINE density in 750kb around each gene
LINEs_file = '/project/lncrna/Xist/data/annotation_files/LINEs/LINEs_mm9.bed'
LINEs_table = read.table(LINEs_file,sep='\t')
colnames(LINEs_table) = c("chr","start","end","name","score","strand","dispStart","dispEnd","col")
LINEs_table = LINEs_table[LINEs_table$chr == 'chrX',]

tmp_LINEs_file = paste(tmp_dir,'/tmp_LINEs_file.bed',sep='')
export.bed(LINEs_table,tmp_LINEs_file,format='bed')

tmp_LINE_density_file = paste(tmp_dir,'/tmp_gene_density.txt',sep='')
cmd = paste(bedTools,'windowBed -a ',tmp_gene_regions_file,' -b ',tmp_LINEs_file,' -w 375000 -c > ',tmp_LINE_density_file,sep='') #-w 	Base pairs added upstream and downstream of each entry in A when searching for overlaps in B. Default is 1000 bp.
system(cmd)
LINE_density_table = read.table(tmp_LINE_density_file,header=F,sep='\t')[,c(4,7)]
colnames(LINE_density_table) = c('gene_name','LINE_density')

LINEs = c(LINEs_table$start,LINEs_table$end)
all_distances = abs(outer(halftimes_table$TSS,LINEs,FUN="-"))
distance_LINE_table = data.frame(gene_name = halftimes_table$gene_name, distance_to_LINE = apply(all_distances,1,min))


# ##########PRO-seq RPKM at time point 0
# proseq_table = read.table(file = "/project/lncrna/Xist/data_lisa/silencing_halftimes/raw_data/PROseq.txt", sep='\t', header=T)
# proseq_rpkm_table = data.frame(gene_name=proseq_table$Genes, proseq_rpkm = rowMeans(cbind(proseq_table$NoDoxA_RPKM,proseq_table$NoDoxB_RPKM)))


# ##########Pausing Index
# pausing_index_table = read.table(file = "/project/lncrna/Xist/data_lisa/annotation_files/pausing_data/NoDox_all_highPI_pausedata", sep='\t', header=T)
# pausing_index_table = pausing_index_table[pausing_index_table$chr == 'chrX',]
# pausing_index_table = pausing_index_table[pausing_index_table$PI != -1,]
#
# #only keep one transcript per gene (most 5') because we do not have transcript information for silencing genes
# duplicates = pausing_index_table$name[duplicated(pausing_index_table$name)]
# for(i in 1:length(duplicates)){
#   duplicated_gene_name = duplicates[i]
#   transcripts = pausing_index_table[pausing_index_table$name == duplicated_gene_name,]
#   five_prime_transcript = transcripts[which.min(transcripts$pr_start),]
#   pausing_index_table = pausing_index_table[!(pausing_index_table$name == duplicated_gene_name),]
#   pausing_index_table = rbind(pausing_index_table,five_prime_transcript)
# }
#
# pausing_index_table = data.frame(gene_name = pausing_index_table$name, promoter_density = pausing_index_table$prDens,
#                                  gene_density_pro_seq = pausing_index_table$gbDens, pausing_index = pausing_index_table$PI)
#

# ##########bivalency
# bivalency_table = read.table(file = "/project/lncrna/Xist/data_lisa/annotation_files/bivalency/bivalncy_table.csv", sep='\t', header=T)
# bivalency_table$start = 0
#
# #only keep one transcript per gene (most 5') because we do not have transcript information for silencing genes
# duplicates = bivalency_table$Gene.s.[duplicated(bivalency_table$Gene.s.)]
# for(i in 1:length(duplicates)){
#   duplicated_gene_name = duplicates[i]
#   transcripts = bivalency_table[bivalency_table$Gene.s. == duplicated_gene_name,]
#   for(j in 1:nrow(transcripts)){
#     transcripts$start[j] = as.numeric(unlist(strsplit(unlist(strsplit(as.character(transcripts$TSS.Pos[j]),split = ':'))[2],split = '-'))[1])
#   }
#   five_prime_transcript = transcripts[which.min(transcripts$start),]
#   bivalency_table = bivalency_table[!(bivalency_table$Gene.s. == duplicated_gene_name),]
#   bivalency_table = rbind(bivalency_table,five_prime_transcript)
# }
#
# bivalency_table = data.frame(gene_name = bivalency_table$Gene.s., CpG_class = bivalency_table$Class, bivalent = bivalency_table$ESC)
# bivalency_table$CpG_class =as.character(bivalency_table$CpG_class)
# bivalency_table$CpG_class[bivalency_table$CpG_class != 'HCP'] = 0
# bivalency_table$CpG_class[bivalency_table$CpG_class == 'HCP'] = 1
# bivalency_table$CpG_class =as.numeric(bivalency_table$CpG_class)
# bivalency_table$bivalent = as.character(bivalency_table$bivalent)
# bivalency_table$bivalent[bivalency_table$bivalent != 'K4+K27'] = 0
# bivalency_table$bivalent[bivalency_table$bivalent == 'K4+K27'] = 1
# bivalency_table$bivalent = as.numeric(bivalency_table$bivalent)


##########RNA-Seq at time point 0 --> not mappable bause of missing transcript information for rates and missing start/stop information to select the most 5' one
# rnaseq_rpkm_table = read.table(file = "/project/lncrna/Xist/old/Data/all_refgene_rpkm_final.data", sep='\t', header=T)
# rnaseq_rpkm_table = na.omit(rnaseq_rpkm_table)
# rnaseq_rpkm_table = data.frame(gene_name=rnaseq_rpkm_table$gene_name, trans_name = rnaseq_rpkm_table$trans_name,
#                                rnaseq_rpkm=rowMeans(cbind(rnaseq_rpkm_table$A134T9_RPKM,rnaseq_rpkm_table$A134T10_RPKM)))


##########remove tmp files
cmd = paste('rm -rf',tmp_dir)
system(cmd)


##########merge all tables
merged_tables = merge(halftimes_table[,c(4,5)], chip_seq_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, distance_TAD_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, gene_density_gencode_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, overlap_initiation_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, distance_xist_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, overlap_LADs_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, distance_LAD_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, hic_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, hicap_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, dna_methylation_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, CpG_island_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, CpG_content_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, LINE_density_table, by = 'gene_name', all.x = T, all.y = F)
merged_tables = merge(merged_tables, distance_LINE_table, by = 'gene_name', all.x = T, all.y = F)
#merged_tables = merge(merged_tables, proseq_rpkm_table, by = 'gene_name', all.x = T, all.y = F)
#merged_tables = merge(merged_tables, pausing_index_table, by = 'gene_name', all.x = T, all.y = F)
#merged_tables = merge(merged_tables, bivalency_table, by = 'gene_name', all.x = T, all.y = F)

#for all empty entries (NA) write a 0 into data frame
#merged_tables[is.na(merged_tables)] = 0

#create matrix
matrix_table = merged_tables
row.names(matrix_table) = matrix_table$gene_name
matrix_table$gene_name = NULL
print(nrow(matrix_table))

#remove halftime from tablebecause it is used as target
halftime = matrix_table$halftime
matrix_table$halftime = NULL

#make certain features as factors for nomalization
matrix_table$overlap_with_Xist_entry_sites = as.factor(matrix_table$overlap_with_Xist_entry_sites)
matrix_table$overlap_with_CpG_islands = as.factor(matrix_table$overlap_with_CpG_islands)
matrix_table$overlap_with_LADs = as.factor(matrix_table$overlap_with_LADs)
#matrix_table$CpG_class = as.factor(matrix_table$CpG_class)
#matrix_table$bivalent = as.factor(matrix_table$bivalent)


#########Save data_set and target
data_set = matrix_table
halftime = halftime

save(data_set,halftime,file = paste(output_dir,region_name,'_matrix.RData',sep=''))


#standardizing
data_set = matrix_table
halftime = halftime
for(j in 1:(ncol(data_set))){
  if(!is.factor(data_set[,j])){
    data_set[,j] = (data_set[,j] - mean(data_set[,j])) / sd(data_set[,j])
  }
}
save(data_set,halftime,file = paste(output_dir,region_name,'_matrix_standardized.RData',sep=''))


#log transformation
data_set = matrix_table
halftime = halftime
for(j in 1:(ncol(data_set))){ 
  if(!is.factor(data_set[,j])){
    data_set[,j] = log(data_set[,j]+1)
  }
}

save(data_set,halftime,file = paste(output_dir,region_name,'_matrix_log_transformed.RData',sep=''))


  

#######boxplots
# hic_inteactions = c("mean_interaction_strength_HiC_all","mean_interaction_strength_HiC_promoter","mean_interaction_strength_HiC_xist")
# 
# pdf(paste("/project/lncrna/Xist/plots/additional_analysis/basic_analysis/boxplots_",region_name,".pdf", sep=""), 8.5,7)
# 
# thr_silencing_lower = 0.5
# thr_silencing_middle = 0.9
# thr_silencing_uppper = 1.3
# 
# for(i in 1:ncol(matrix_table)){
#   
#   print(colnames(matrix_table)[i])
#   column = matrix_table[,i]
#   feature = colnames(matrix_table)[i]
#   
#   if(is.factor(column)){
# 
#     column_target = cbind.data.frame(column,halftime)
#     column_target = column_target[complete.cases(column_target),]
# 
#     feature0 = column_target[column_target[,1]==0,] # the feature
#     feature1 = column_target[column_target[,1]==1,] # the feature
#     wilcox = wilcox.test(feature0[,2],feature1[,2])$p.value
#     
#     par(mar=c(15,7,7,7))
#     title_box_plot = paste(feature,"\nwilcox p-value:",round(wilcox,5))
#     boxplot(feature1[,2],feature0[,2], names=c("overlap\n(interaction)", "no overlap\n(no interaction)"), main=title_box_plot, ylab="Half time", lwd=2, cex.axis = 1.75, cex.main = 1.9, cex.lab = 1.9, las=2)
#     
#     
#   }else{
#     
#     #exclude outliers from the data 
#     if(sum(range(column) == range(c(0,1))) != 2){
#       qnt=quantile(column,probs=c(.05,0.95),na.rm = T)
#       column[column<(qnt[1])] <-NA
#       column[column>(qnt[2])] <-NA 
#     }
#     
#     #if removal of outliers destroys bimodal distribution(e.g. binary variables) then don't remove outliers
#     if(length(unique(column[complete.cases(column)])) == 1){
#       column = matrix_table[,i]
#     }
#     column_target = cbind.data.frame(column,halftime)
#     column_target = column_target[complete.cases(column_target),]
#     
#     #calculate overall correlations
#     sper_corr = cor(column_target[,1],column_target[,2], method="spearman")
#     sper_pValue = cor.test(column_target[,1],column_target[,2], method="spearman")$p.value
#     
#     pears_corr = cor(column_target[,1],column_target[,2])
#     pears_pValue = cor.test(column_target[,1],column_target[,2])$p.value
#     
#     dist_corr = dcor(column_target[,1],column_target[,2])
#     
#     #calculate correlation for genes in classes
#     if(thr_silencing_middle == "-"){
#       class1 = column_target[column_target[,2] < thr_silencing_lower,] #silenced/early silenced
#       class2 = column_target[column_target[,2] > thr_silencing_uppper,] #not silenced/late silenced
#       column_target_class = rbind(column_target[column_target[,2] < thr_silencing_lower,], column_target[column_target[,2] > thr_silencing_uppper,])
#     }else{
#       class1 = column_target[column_target[,2] < thr_silencing_lower,]
#       class2 = column_target[column_target[,2] > thr_silencing_middle & column_target[,2] < thr_silencing_uppper,]
#       column_target_class = rbind(column_target[column_target[,2] < thr_silencing_lower,], column_target[column_target[,2] > thr_silencing_middle & column_target[,2] < thr_silencing_uppper,])
#     }
#     
#     sper_corr_class = cor(column_target_class[,1],column_target_class[,2], method="spearman")
#     sper_pValue_class = cor.test(column_target_class[,1],column_target_class[,2], method="spearman")$p.value
#     
#     pears_corr_class = cor(column_target_class[,1],column_target_class[,2])
#     pears_pValue_class = cor.test(column_target_class[,1],column_target_class[,2])$p.value
#     
#     dist_corr_class = dcor(column_target_class[,1],column_target_class[,2])
#     
#     wilcox = wilcox.test(class1[,1],class2[,1])$p.value
#     
#     #print(paste(feature,pears_corr,sper_corr,dist_corr,sep=","))
#     
#     par(mar=c(10,7,10,7))
#     density_class1 = density(class1[,1])
#     density_class2 = density(class2[,1])
#     plot(density_class1, main = feature, ylim = c(0,max(c(density_class1$y,density_class2$y))), col="lightblue")
#     lines(density_class2,col="orange")
#     legend("topright", legend = c("class1: silenced/early silenced","class2: not silenced/late silenced"), col=c("lightblue","orange"),pch = "*")
#     
#     par(mar=c(15,7,7,7))
#     title_box_plot = paste(feature,"\npearson cor p-value (all genes with halftime):",round(pears_pValue,5),"\n pearson cor p-value (only class genes):",round(pears_pValue_class,5))
#     boxplot(names=c("class1 \n(silenced/early silenced)","class2 \n(not silenced/late silenced)"), class1[,1],class2[,1], main=title_box_plot, lwd=2, cex.axis = 1, cex.main = 1.9, cex.lab = 1.9, las=2)
#     
#   }
#   if(feature %in% hic_inteactions){
#     
#     column_target = cbind.data.frame(column,halftime)
#     column_target = column_target[complete.cases(column_target),]
#     column_target$column[column_target$column > 0] = 1
#     column_target$column = as.factor(column_target$column)
#     
#     feature0 = column_target[column_target[,1]==0,] # the feature
#     feature1 = column_target[column_target[,1]==1,] # the feature
#     wilcox = wilcox.test(feature0[,2],feature1[,2])$p.value
#     
#     par(mar=c(15,7,7,7))
#     title_box_plot = paste(feature,"\nwilcox p-value:",round(wilcox,5))
#     boxplot(feature1[,2],feature0[,2], names=c("overlap\n(interaction)", "no overlap\n(no interaction)"), main=title_box_plot, ylab="Half time", lwd=2, cex.axis = 1.75, cex.main = 1.9, cex.lab = 1.9, las=2)
#   }
# }
# dev.off()







