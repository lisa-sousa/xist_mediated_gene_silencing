
################
#functions
################

remove_na_columns <- function(x){
  x = x[sapply(x, function(x) !any(is.na(x)))] 
  return(x)
}

get_expression_pro_seq <- function(data){
  RPKM_index = grep("RPKM",colnames(data))
  ratio_index = grep("Ratio",colnames(data))
  b6_expression = remove_na_columns(data[,RPKM_index]*data[,ratio_index])
  b6_expression = data.frame(NoDox_RPKM = rowMeans(as.data.frame(b6_expression[,grep("NoDox",colnames(b6_expression))])),b6_expression[-c(grep("NoDox",colnames(b6_expression)),grep("24h",colnames(b6_expression)))],Dox24h_RPKM = rowMeans(as.data.frame(b6_expression[,grep("24h",colnames(b6_expression))])))
  cast_expression = remove_na_columns(data[,RPKM_index]*(1-data[,ratio_index]))
  cast_expression = data.frame(NoDox_RPKM = rowMeans(as.data.frame(cast_expression[,grep("NoDox",colnames(cast_expression))])),cast_expression[-c(grep("NoDox",colnames(cast_expression)),grep("24h",colnames(cast_expression)))],Dox24h_RPKM = rowMeans(as.data.frame(cast_expression[,grep("24h",colnames(cast_expression))])))
  return(list(b6_expression,cast_expression))
}

get_expression_mrna_seq_undiff <- function(data){
  RPKM_index = grep("RPKM",colnames(data))
  ratio_index = RPKM_index-1
  b6_expression = remove_na_columns(data[,RPKM_index]*data[,ratio_index])
  cast_expression = remove_na_columns(data[,RPKM_index]*(1-data[,ratio_index]))
  return(list(b6_expression,cast_expression))
}

get_expression_mrna_seq_diff <- function(data){
  RPKM_index = grep("RPKM",colnames(data))
  ratio_index = RPKM_index-1
  b6_expression = remove_na_columns(data[,RPKM_index]*data[,ratio_index])
  b6_expression = data.frame(nodox = rowMeans(b6_expression[,c(1,2)]),dox8 = rowMeans(b6_expression[,c(3,4)]),dox16 = rowMeans(b6_expression[,c(5,6)]),
                             dox24 = rowMeans(b6_expression[,c(7,8)]),dox48 = rowMeans(b6_expression[,c(9,10)]))
  
  cast_expression = remove_na_columns(data[,RPKM_index]*(1-data[,ratio_index]))
  cast_expression = data.frame(nodox = rowMeans(cast_expression[,c(1,2)]),dox8 = rowMeans(cast_expression[,c(3,4)]),dox16 = rowMeans(cast_expression[,c(5,6)]),
                               dox24 = rowMeans(cast_expression[,c(7,8)]),dox48 = rowMeans(cast_expression[,c(9,10)]))
  return(list(b6_expression,cast_expression))
}

plot_expression <- function(time,expression,gene){
  plot(time,expression[[1]],pch=20,col="darkgreen",main=paste(gene,"expression over time"),ylab=paste(gene,"expression [RPKM]"),xlab="time [days]",ylim=c(0,max(sapply(expression, max))))
  points(time,expression[[2]],pch=20,col="darkblue")
  lines(time,expression[[1]],col="darkgreen")
  lines(time,expression[[2]],col="darkblue")
}

################
#PRO-Seq
################

pdf("/project/lncrna/Xist/plots/additional_analysis/allelic_expression.pdf")

pro_seq_file = "/project/lncrna/Xist/data/silencing_halftimes/raw_data/PROseq.txt"
pro_seq_data = read.table(pro_seq_file,na.strings=c("NA","nd"),header=T,sep='\t')
time = c(0,0.5,1,2,4,8,12,24)/24

####Xist
xist_data = pro_seq_data[pro_seq_data$Genes == "Xist",]
expression = get_expression_pro_seq(xist_data)
plot_expression(time,expression,"Xist")


####Tsix
tsix_data = pro_seq_data[pro_seq_data$Genes == "Tsix",]
expression = get_expression_pro_seq(tsix_data)
plot_expression(time,expression,"Tsix")

dev.off()

################
#mRNA-Seq undiff
################

pdf("/project/lncrna/Xist/plots/additional_analysis/allelic_expression_undifferentiated.pdf")

mrna_seq_file = "/project/lncrna/Xist/data/silencing_halftimes/raw_data/UndifferenciatedCells.txt"
mrna_seq_file = read.table(mrna_seq_file,na.strings=c("NA","nd"),header=T,sep='\t')
time = c(0,2,4,8,12,24)/24
first_replicate = c(1,4,5,8,9,12,13,16,17,20,21,24,25)
mrna_seq_file = mrna_seq_file[,first_replicate]

####Xist
xist_data = mrna_seq_file[mrna_seq_file$Genes == "Xist",]
expression = get_expression_mrna_seq_undiff(xist_data)
plot_expression(time,expression,"Xist")


####Tsix
tsix_data = mrna_seq_file[mrna_seq_file$Genes == "Tsix",]
expression = get_expression_mrna_seq_undiff(tsix_data)
plot_expression(time,expression,"Tsix")

dev.off()


################
#mRNA-Seq diff
################

pdf("/project/lncrna/Xist/plots/additional_analysis/allelic_expression_differentiated.pdf")

mrna_seq_file = "/project/lncrna/Xist/data/silencing_halftimes/raw_data/DifferenciatedCells.txt"
mrna_seq_file = read.table(mrna_seq_file,na.strings=c("NA","nd"),header=T,sep='\t')
time = c(0,8,16,24,48)/24

####Xist
xist_data = mrna_seq_file[mrna_seq_file$Genes == "Xist",]
expression = get_expression_mrna_seq_diff(xist_data)
plot_expression(time,expression,"Xist")


####Tsix
tsix_data = mrna_seq_file[mrna_seq_file$Genes == "Tsix",]
expression = get_expression_mrna_seq_diff(tsix_data)
plot_expression(time,expression,"Tsix")

dev.off()


