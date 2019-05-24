###################################################
#tools
###################################################

library(here)
bedTools = '/home/lisasous/tools/bedtools2/bin/' #set path to bedtools 2 /bin/ folder

###################################################
#input/output files and directories
###################################################

dir_gene_annotation = here("/data/annotation_files/gene_annotation/")
file_gencode = paste(dir_gene_annotation,"gencode.vM9.annotation.gtf",sep="")

file_transcript_annotation = paste(dir_gene_annotation,"gencode.vM9.annotation.chrX.transcripts.bed",sep="")
file_transcript_annotation_200bp_aroundTSS = paste(dir_gene_annotation,"gencode.vM9.annotation.chrX.transcripts.200bp.aroundTSS.bed",sep="")
file_transcript_annotation_2000bp_aroundTSS = paste(dir_gene_annotation,"gencode.vM9.annotation.chrX.transcripts.2000bp.aroundTSS.bed",sep="")


dir_regulatory_regions = here("data/annotation_files/regulatory_regions/")
file_rr_replicate1 = paste(dir_regulatory_regions,"NoDoxA_0.8Named_regs250merge.bed",sep="")
file_rr_replicate2 = paste(dir_regulatory_regions,"NoDoxB_0.8Named_regs250merge.bed",sep="")

file_rr_merged = paste(dir_regulatory_regions,"NoDox_A_B_chrX_0.8Named_regs250merge.bed",sep="")
file_rr_sorted = paste(dir_regulatory_regions,"NoDox_A_B_chrX_0.8Named_regs250merge.sorted.bed",sep="")
file_rr_collapsed = paste(dir_regulatory_regions,"NoDox_merged_chrX_0.8Named_regs250merge.bed",sep="")


file_overlap_transcript_200bp_aroundTSS_rr = paste(dir_gene_annotation,"gencode.vM9.annotation.chrX.transcripts.200bp.aroundTSS.overlap.rr.bed",sep="")
file_overlap_transcript_2000bp_aroundTSS_rr = paste(dir_gene_annotation,"gencode.vM9.annotation.chrX.transcripts.2000bp.aroundTSS.overlap.rr.bed",sep="")

file_output = paste(dir_gene_annotation,"gencode.vM9.annotation.chrX.genes.reannotated.with.rr.bed",sep="")

####################################################
#get transcript annotation from gencode for mm10
####################################################

#process downloaded gencode file https://www.gencodegenes.org/mouse_releases/
gene_annotation_gencode = read.table(file=file_gencode, sep = "\t")
gene_annotation_gencode = gene_annotation_gencode[gene_annotation_gencode$V3 == "transcript",]
gene_annotation_gencode = gene_annotation_gencode[gene_annotation_gencode$V1 == "chrX",]
gene_annotation_gencode = gene_annotation_gencode[,c(1,4,5,7,9)]
gene_annotation_gencode$V9 = as.character(gene_annotation_gencode$V9)
gene_annotation_gencode$transcript_id = ''
gene_annotation_gencode$gene_name = ''
gene_annotation_gencode$gene_type = ''


for(i in 1:nrow(gene_annotation_gencode)){
  info = unlist(strsplit(gene_annotation_gencode$V9[i],split = "; "))
  transcript_id = unlist(strsplit(info[2],split = " "))[2]
  gene_type = unlist(strsplit(info[3],split = " "))[2]
  gene_name = unlist(strsplit(info[5],split = " "))[2]
  gene_annotation_gencode$transcript_id[i] = transcript_id
  gene_annotation_gencode$gene_type[i] = gene_type
  gene_annotation_gencode$gene_name[i] = gene_name
}
gene_annotation_gencode$V9 = NULL
#gene_annotation_gencode = gene_annotation_gencode[gene_annotation_gencode$gene_type == "protein_coding",]

gene_annotation_gencode = gene_annotation_gencode[,c(1,2,3,5,6,4,7)]
colnames(gene_annotation_gencode) = c("chr","start","end","transcript_id","gene_name","strand","gene_type")
write.table(gene_annotation_gencode, file=file_transcript_annotation, sep='\t', quote = F, col.names = F, row.names = F)

get_xbp_around_TSS = function(gene_annotation_xbp_around_TSS, x, output_file){
  gene_annotation_xbp_around_TSS$TSS[gene_annotation_xbp_around_TSS$strand == "-"] = gene_annotation_xbp_around_TSS$end[gene_annotation_xbp_around_TSS$strand == "-"]
  gene_annotation_xbp_around_TSS$TSS[gene_annotation_xbp_around_TSS$strand == "+"] = gene_annotation_xbp_around_TSS$start[gene_annotation_xbp_around_TSS$strand == "+"]
  gene_annotation_xbp_around_TSS$start = gene_annotation_xbp_around_TSS$TSS - x
  gene_annotation_xbp_around_TSS$end = gene_annotation_xbp_around_TSS$TSS + x
  gene_annotation_xbp_around_TSS$TSS = NULL
  write.table(gene_annotation_xbp_around_TSS, file = output_file, sep='\t' ,quote = F, col.names = F, row.names = F)
  
}

get_xbp_around_TSS(gene_annotation_gencode, 100, file_transcript_annotation_200bp_aroundTSS)
get_xbp_around_TSS(gene_annotation_gencode, 1000, file_transcript_annotation_2000bp_aroundTSS)

####################################################
#process regulatory regions
####################################################

#merge replicates
rr_replicate1 = read.table(file=file_rr_replicate1, sep = "\t")
rr_replicate2 = read.table(file=file_rr_replicate2, sep = "\t")

rr_merged = rbind(rr_replicate1, rr_replicate2)
rr_merged = rr_merged[rr_merged$V1 == "chrX", ]
write.table(rr_merged, file = file_rr_merged, sep='\t' ,quote = F, col.names = F, row.names = F)

#sort bed file and merge ovelapping rr and use the mean strength of merged rr
cmd = paste("sort -k1,1 -k2,2n",file_rr_merged,">",file_rr_sorted)
system(cmd)

cmd = paste(bedTools, "bedtools merge -i ",file_rr_sorted," -c 4,5,6 -o collapse,mean,distinct > ",file_rr_collapsed, sep="")
system(cmd)


####################################################
#overlap transcript annotation with regulatory regions
####################################################

cmd = paste(bedTools, "bedtools intersect -wao -a ",file_transcript_annotation_200bp_aroundTSS," -b ",file_rr_collapsed," > ",file_overlap_transcript_200bp_aroundTSS_rr,sep="")
system(cmd)

cmd = paste(bedTools, "bedtools intersect -wao -a ",file_transcript_annotation_2000bp_aroundTSS," -b ",file_rr_collapsed," > ",file_overlap_transcript_2000bp_aroundTSS_rr,sep="")
system(cmd)

####################################################
#get active promoters
####################################################

overlap_table = read.table(file_overlap_transcript_200bp_aroundTSS_rr)

original_transcript_annotation = read.table(file_transcript_annotation) #original transcript annotation (not 200bp surrounding TSS)
original_transcript_annotation$transcript_length = original_transcript_annotation$V3 - original_transcript_annotation$V2

transcript_annotation = merge(overlap_table,original_transcript_annotation,by="V4")
transcript_annotation = transcript_annotation[,c(2,16,17,1,5,6,7,14,9,10,12,21)]
colnames(transcript_annotation) = c("chr","start","end","transcript_id","gene_name","strand","gene_type","overlap_rr","start_rr","end_rr","score_rr","transcript_length")
transcript_annotation$gene_name = as.character(transcript_annotation$gene_name)
transcript_annotation$TSS[transcript_annotation$strand == "+"] = transcript_annotation$start[transcript_annotation$strand == "+"]
transcript_annotation$TSS[transcript_annotation$strand == "-"] = transcript_annotation$end[transcript_annotation$strand == "-"]
transcript_annotation = transcript_annotation[order(transcript_annotation$gene_name),]

#get list of genes (not transcripts)
gene_list = unique(transcript_annotation$gene_name)

#sort transcripts into three lists
one_entry_genes = NULL  #genes with one transcript overlapping a rr
multiple_entires_genes = NULL  #genes with multiple transcripts overlapping multiple rr
no_overlap_genes = NULL #gene without any transcript overlapping a rr

for(i in 1:length(gene_list)){
  all_transcripts_for_gene_i = transcript_annotation[transcript_annotation$gene_name == gene_list[i],]
  if(nrow(all_transcripts_for_gene_i) > 1){
    num_transcripts_overlap_rr = nrow(all_transcripts_for_gene_i[all_transcripts_for_gene_i$overlap_rr > 0,])
    if(num_transcripts_overlap_rr > 0){
      all_transcripts_for_gene_i = all_transcripts_for_gene_i[all_transcripts_for_gene_i$overlap_rr > 0,]
      if(nrow(all_transcripts_for_gene_i) > 1){
        multiple_entires_genes = rbind(multiple_entires_genes,all_transcripts_for_gene_i)
      }else{one_entry_genes = rbind(one_entry_genes,all_transcripts_for_gene_i)}
    }else{no_overlap_genes = rbind(no_overlap_genes,all_transcripts_for_gene_i)}
  }else{
    if(all_transcripts_for_gene_i$overlap_rr > 0){
      one_entry_genes = rbind(one_entry_genes,all_transcripts_for_gene_i)
    }else{no_overlap_genes = rbind(no_overlap_genes,all_transcripts_for_gene_i)}
  }
}

print(length(unique(one_entry_genes$gene_name)))
print(length(unique(multiple_entires_genes$gene_name)))
print(length(unique(no_overlap_genes$gene_name)))


#if multiple TSS have same start coordinates only keep the TSS with longest transcript
gene_list_multiple_entry_genes = as.character(unique(multiple_entires_genes$gene_name))
multiple_entires_genes_filtered_longest_transcript = NULL

for(i in 1:length(gene_list_multiple_entry_genes)){
  all_transcripts_for_gene_i = multiple_entires_genes[multiple_entires_genes$gene_name == gene_list_multiple_entry_genes[i],]
  unique_TSS = unique(all_transcripts_for_gene_i$TSS)
  longest_transcipts = NULL
  for(j in 1:length(unique_TSS)){
    transcripts_with_same_TSS = all_transcripts_for_gene_i[all_transcripts_for_gene_i$TSS == unique_TSS[j],]
    longest_transcipt = transcripts_with_same_TSS[which.max(transcripts_with_same_TSS$transcript_length),]
    longest_transcipts = rbind(longest_transcipts,longest_transcipt)
  }
  if(nrow(longest_transcipts) > 1){
    multiple_entires_genes_filtered_longest_transcript = rbind(multiple_entires_genes_filtered_longest_transcript,longest_transcipts)
  }else{one_entry_genes = rbind(one_entry_genes,longest_transcipts)}
}

print(length(unique(one_entry_genes$gene_name)))
print(length(unique(multiple_entires_genes_filtered_longest_transcript$gene_name)))

#if multiple TSS with different start coordinates take the one with strongest rr
gene_list_multiple_entry_genes = as.character(unique(multiple_entires_genes_filtered_longest_transcript$gene_name))
multiple_entires_genes_filtered_strongest_rr = NULL

for(i in 1:length(gene_list_multiple_entry_genes)){
  all_transcripts_for_gene_i = multiple_entires_genes_filtered_longest_transcript[multiple_entires_genes_filtered_longest_transcript$gene_name == gene_list_multiple_entry_genes[i],]
  strongest_rr = max(unique(all_transcripts_for_gene_i$score_rr))
  transcripts_at_strongest_rr = all_transcripts_for_gene_i[all_transcripts_for_gene_i$score_rr == strongest_rr,]
  if(nrow(transcripts_at_strongest_rr) > 1){
    multiple_entires_genes_filtered_strongest_rr = rbind(multiple_entires_genes_filtered_strongest_rr,transcripts_at_strongest_rr)
  }else{one_entry_genes = rbind(one_entry_genes,transcripts_at_strongest_rr)}
}

print(length(unique(one_entry_genes$gene_name)))
print(length(unique(multiple_entires_genes_filtered_strongest_rr$gene_name)))

#if multiple TSS start at the same rr: take TSS that is closest to the center of the rr
multiple_entires_genes_filtered_strongest_rr$center_rr = multiple_entires_genes_filtered_strongest_rr$start_rr + 
                                                         round((multiple_entires_genes_filtered_strongest_rr$end_rr - multiple_entires_genes_filtered_strongest_rr$start_rr)/2)
multiple_entires_genes_filtered_strongest_rr$distance_TSS_center_rr = abs(multiple_entires_genes_filtered_strongest_rr$TSS - multiple_entires_genes_filtered_strongest_rr$center_rr)
gene_list_multiple_entry_genes = as.character(unique(multiple_entires_genes_filtered_strongest_rr$gene_name))

for(i in 1:length(gene_list_multiple_entry_genes)){
  all_transcripts_for_gene_i = multiple_entires_genes_filtered_strongest_rr[multiple_entires_genes_filtered_strongest_rr$gene_name == gene_list_multiple_entry_genes[i],]
  closest_TSS = all_transcripts_for_gene_i[which.min(all_transcripts_for_gene_i$distance_TSS_center_rr),1:13]
  one_entry_genes = rbind(one_entry_genes,closest_TSS)
}

print(length(unique(one_entry_genes$gene_name)))

#for genes where NO TSS overlaps a rr: extend search space to 1000bp around the TSS
overlap_table_2000bp = read.table(file_overlap_transcript_2000bp_aroundTSS_rr)
transcript_annotation_2000bp = merge(overlap_table_2000bp,original_transcript_annotation,by="V4")
transcript_annotation_2000bp = transcript_annotation_2000bp[,c(2,16,17,1,5,6,7,14,9,10,12,21)]
colnames(transcript_annotation_2000bp) = c("chr","start","end","transcript_id","gene_name","strand","gene_type","overlap_rr","start_rr","end_rr","score_rr","transcript_length")
transcript_annotation_2000bp$gene_name = as.character(transcript_annotation_2000bp$gene_name)
transcript_annotation_2000bp$TSS[transcript_annotation_2000bp$strand == "+"] = transcript_annotation_2000bp$start[transcript_annotation_2000bp$strand == "+"]
transcript_annotation_2000bp$TSS[transcript_annotation_2000bp$strand == "-"] = transcript_annotation_2000bp$end[transcript_annotation_2000bp$strand == "-"]
transcript_annotation_2000bp = transcript_annotation_2000bp[order(transcript_annotation_2000bp$gene_name),]

transcript_annotation_no_overlap_genes = transcript_annotation_2000bp[transcript_annotation_2000bp$transcript_id %in% no_overlap_genes$transcript_id ,]
transcript_annotation_no_overlap_genes = transcript_annotation_no_overlap_genes[order(transcript_annotation_no_overlap_genes$gene_name),]

gene_list_no_overlap_genes = unique(transcript_annotation_no_overlap_genes$gene_name[transcript_annotation_no_overlap_genes$overlap_rr == 0])
print(length(gene_list_no_overlap_genes))

transcript_annotation = transcript_annotation_no_overlap_genes[transcript_annotation_no_overlap_genes$overlap_rr != 0,]
gene_list = unique(transcript_annotation$gene_name)
for(i in 1:length(gene_list)){
  all_transcripts_for_gene_i = transcript_annotation[transcript_annotation$gene_name == gene_list[i],]
  TSS_with_strongest_rr = all_transcripts_for_gene_i[which.max(all_transcripts_for_gene_i$score_rr),]
  TSS_with_strongest_rr = TSS_with_strongest_rr[which.max(TSS_with_strongest_rr$transcript_length),]
  new_TSS = round(mean(c(TSS_with_strongest_rr$start_rr,TSS_with_strongest_rr$end_rr)))
  TSS_with_strongest_rr$start[TSS_with_strongest_rr$strand == "+"] = new_TSS
  TSS_with_strongest_rr$end[TSS_with_strongest_rr$strand == "-"] = new_TSS
  if(TSS_with_strongest_rr$end > TSS_with_strongest_rr$start){
    one_entry_genes = rbind(one_entry_genes,TSS_with_strongest_rr)
  }
}

one_entry_genes = one_entry_genes[order(one_entry_genes$gene_name),]
print(length(unique(one_entry_genes$gene_name)))

output = one_entry_genes[,c(1,2,3,5,7,6)]
write.table(output,file=file_output, sep="\t", col.names = F, row.names = F, quote = F)





