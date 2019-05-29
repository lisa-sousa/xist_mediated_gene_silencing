###################################################################################
#libraries
###################################################################################

library(biomaRt)
library(rtracklayer)
library(here)

###################################################################################
#process downloaded gencode file https://www.gencodegenes.org/mouse_releases/
###################################################################################

gene_annotation_gencode = read.table(file=here("data/annotation_files/gene_annotation","gencode.vM9.annotation.gtf"),comment.char = "#",sep="\t")
gene_annotation_gencode = gene_annotation_gencode[gene_annotation_gencode$V3 == "gene",]
gene_annotation_gencode = gene_annotation_gencode[gene_annotation_gencode$V1 == "chrX",]
gene_annotation_gencode = gene_annotation_gencode[,c(1,4,5,7,9)]
gene_annotation_gencode$V9 = as.character(gene_annotation_gencode$V9)
gene_annotation_gencode$gene_name = ''
gene_annotation_gencode$score = "."

for(i in 1:nrow(gene_annotation_gencode)){
  info = unlist(strsplit(gene_annotation_gencode$V9[i],split = "; "))
  status = unlist(strsplit(info[3],split = " "))[2]
  if(status != "KNOWN"){print(status)}
  gene_name = unlist(strsplit(info[4],split = " "))[2]
  gene_annotation_gencode$gene_name[i] = gene_name
}
gene_annotation_gencode$V9 = NULL


gene_annotation_gencode = gene_annotation_gencode[,c(1,2,3,5,6,4)]
colnames(gene_annotation_gencode) = c("chr","start","end","name","score","strand")
gene_annotation_gencode = gene_annotation_gencode[order(gene_annotation_gencode$name),]
write.table(gene_annotation_gencode, file=here("data/annotation_files/gene_annotation","gencode.vM9.annotation.chrX.genes.bed"), sep='\t', quote = F, col.names = F, row.names = F)

###################################################################################
#get genes from ensemble
###################################################################################

##mm9 genome
ensembl <- useEnsembl(biomart="ensembl", version = 67, dataset = "mmusculus_gene_ensembl")
listDatasets(ensembl)[listDatasets(ensembl)$"dataset" == "mmusculus_gene_ensembl",]

##mm10 genome
ensembl <- useEnsembl(biomart="ensembl", version = 87, dataset = "mmusculus_gene_ensembl")
listDatasets(ensembl)[listDatasets(ensembl)$"dataset" == "mmusculus_gene_ensembl",]

#list of possible filters
#listFilters(ensembl)

#list of possible attributes
#listAttributes(ensembl)

#
gene_annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand', 'external_gene_name','ucsc'), 
                         filters = 'chromosome_name', values = 'X', mart = ensembl)

gene_annotation$chromosome_name = paste('chr',gene_annotation$chromosome_name,sep='')
gene_annotation$strand[gene_annotation$strand == -1] = '-'
gene_annotation$strand[gene_annotation$strand == 1] = '+'

###################################################################################
#get genes from ucsc
###################################################################################

session <- browserSession()
genome(session) <- "mm9"

#list of all tracks
#trackNames(session)

query <- ucscTableQuery(session, "UCSC Genes", GRangesForUCSCGenome("mm9", "chrX"))

#list the table names
#tableNames(query)


#get genes
tableName(query) <- "knownGene"

#get a GRanges Object
ucsc_genes = track(query)

#get a table --> downloads all transcripts!!
ucsc_genes = getTable(query)


#get gene names
tableName(query) <- "kgXref"
ucsc_gene_names = getTable(query)


#merge genes and gene names
colnames(ucsc_gene_names)[1] = 'name'    #name is the transcript id --> geneSymbol is the gene name
merged_table = merge(ucsc_genes,ucsc_gene_names,by='name')
ucsc_genes = GRanges(seqnames = merged_table$chrom, ranges = IRanges(start = merged_table$txStart, end = merged_table$txEnd), strand = merged_table$strand)
mcols(ucsc_genes) = data.frame(name = merged_table$geneSymbol)




