
################
#functions
################

remove_na_columns <- function(x){
  x = x[sapply(x, function(x) !any(is.na(x)))] 
  return(x)
}

get_expression <- function(data){
  RPKM_index = grep("RPKM",colnames(data))
  ratio_index = grep("Ratio",colnames(data))
  b6_expression = remove_na_columns(data[,RPKM_index]*data[,ratio_index])
  b6_expression = data.frame(NoDox_RPKM = rowMeans(as.data.frame(b6_expression[,grep("NoDox",colnames(b6_expression))])),b6_expression[-c(grep("NoDox",colnames(b6_expression)),grep("24h",colnames(b6_expression)))],Dox24h_RPKM = rowMeans(as.data.frame(b6_expression[,grep("24h",colnames(b6_expression))])))
  cast_expression = remove_na_columns(data[,RPKM_index]*(1-data[,ratio_index]))
  cast_expression = data.frame(NoDox_RPKM = rowMeans(as.data.frame(cast_expression[,grep("NoDox",colnames(cast_expression))])),cast_expression[-c(grep("NoDox",colnames(cast_expression)),grep("24h",colnames(cast_expression)))],Dox24h_RPKM = rowMeans(as.data.frame(cast_expression[,grep("24h",colnames(cast_expression))])))
  return(list(b6_expression,cast_expression))
}

plot_expression <- function(time,expression,gene){
  plot(time,expression[[1]],pch=20,col="darkgreen",main=paste(gene,"expression over time"),ylab=paste(gene,"expression [RPKM]"),xlab="time [days]",ylim=c(0,max(sapply(expression, max))))
  points(time,expression[[2]],pch=20,col="darkblue")
  lines(time,expression[[1]],col="darkgreen")
  lines(time,expression[[2]],col="darkblue")
}

################
#output
################

pro_seq_file = "/project/lncrna/Xist/data/silencing_halftimes/raw_data/PROseq.txt"
pro_seq_data = read.table(pro_seq_file,na.strings=c("NA","nd"),header=T,sep='\t')
time = c(0,0.5,1,2,4,8,12,24)/24

pdf("/project/lncrna/Xist/plots/additional_analysis/allelic_expression.pdf")

####Xist
xist_data = pro_seq_data[pro_seq_data$Genes == "Xist",]
expression = get_expression(xist_data)
plot_expression(time,expression,"Xist")


####Tsix
tsix_data = pro_seq_data[pro_seq_data$Genes == "Tsix",]
expression = get_expression(tsix_data)
plot_expression(time,expression,"Tsix")

dev.off()
