
dir = "/project/lncrna/Xist/data/annotation_files/enhancers/"
file_HiCap_enhancers = "Promoter_Enhancer_Interactions.txt"
file_output = "gene_enhancers.bed"

hicap_enhancers = read.table(paste(dir,file_HiCap_enhancers,sep=""),header = T, sep = "\t")
hicap_enhancers_chrX = hicap_enhancers[hicap_enhancers$Promoter.chr == "chrX" & hicap_enhancers$Fragment.chromosome == "chrX",]
enhancer_bed = data.frame(chr = hicap_enhancers_chrX$Fragment.chromosome, start = hicap_enhancers_chrX$Fragment.start.coordinate, end = hicap_enhancers_chrX$Fragment.end.coordinate,
                          Genes = hicap_enhancers_chrX$Gene.Name, score = rowMeans(cbind(hicap_enhancers_chrX$Read.pairs.in.HiCap.replicate.1,hicap_enhancers_chrX$Read.pairs.in.HiCap.replicate.2)),
                          strand = hicap_enhancers_chrX$Promoter.Strand)

enhancer_length = abs(enhancer_bed$start-enhancer_bed$end)
enhancer_center = rowMeans(enhancer_bed[,2:3])
enhancer_bed$start[enhancer_length < 1000] = enhancer_center[enhancer_length < 1000] - 500
enhancer_bed$end[enhancer_length < 1000] = enhancer_center[enhancer_length < 1000] + 500

write.table(enhancer_bed,paste(dir,file_output,sep=""),col.names = F, row.names = F, sep="\t",quote = F)


#######
#-> map chip-seq data
#-> create data matrix
#######


source("/project/lncrna/Xist/scripts/modelling/model_functions.R")
output_dir = "/project/lncrna/Xist/plots/additional_analysis/"
file_halftimes = "/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_pro_seq_mm10_RSS_initial_ratio.txt"
file_feature_matrix = "/project/lncrna/Xist/data/modelling/feature_matrix/promoter_matrix_normRAdjusted_enhancer.RData"


thr_silencing_lower = 0.9
thr_silencing_middle = "-"
thr_silencing_upper = 1.6


halftimes = read.table(file_halftimes,header=T)

load(file_feature_matrix)
colnames(data_set)[1] = "Genes"
data_set$enhancer_strength = halftime
data_set = merge(halftimes,data_set,by="Genes")
data_set = data_set[data_set$halftime < thr_silencing_lower | data_set$halftime > thr_silencing_upper,]
data_set$target = 0
data_set$target[data_set$halftime > 2] = 1


data=list()
data[[1]] = data_set[,10:86]
data[[2]] = as.factor(data_set$target)
data[[3]] = data_set$halftime

plot_data_boxpots(output_dir, data)

genes = unique(as.character(data_set$Genes[duplicated(data_set$Genes)]))
for(i in 1:length(genes)){
  entry = data_set[data_set$Genes==genes[i],]
  max = max(entry$enhancer_strength)
  data_set = data_set[!(data_set$Genes==genes[i] & data_set$enhancer_strength != max),]
}
data_set = data_set[!duplicated(data_set),]  

  
data=list()
data[[1]] = data_set[,10:86]
data[[2]] = as.factor(data_set$target)
data[[3]] = data_set$halftime

plot_data_boxpots(output_dir, data)
  
  
