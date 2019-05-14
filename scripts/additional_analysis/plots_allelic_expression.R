###################################################################################
#functions
###################################################################################

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
  library(ggplot2)
  library(gridExtra)
  library(cowplot)
  data_plot = data.frame(allele = rep(c("B6","Cast"),each=length(time)),time=c(time,time),expression = unlist(c(expression[[1]],expression[[2]])))
  ggplot = ggplot(data_plot, aes(x = time, y = expression, group = allele)) + 
    geom_line(aes(color=allele)) +
    geom_point(aes(color=allele)) +
    scale_color_manual(values=c("#71c837", "#003380")) +
    theme_minimal(base_family = "Source Sans Pro") + 
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),legend.title = element_blank(), legend.text = element_text(size=8),
          axis.text=element_text(size=8), axis.title=element_text(size=8)) + 
    scale_x_continuous(limits=c(0,max(time)), breaks=c(0,max(time)/2,max(time)), name='time [days]') +
    scale_y_continuous(name=substitute({italic(gene)} * " expression [RPKM]", list(gene=gene))) 
  return(ggplot)
}

###################################################################################
#plot allele-specific expression for PRO-seq data
###################################################################################


pro_seq_file = "/project/lncrna/Xist/data/silencing_halftimes/raw_data/PROseq.txt"
pro_seq_data = read.table(pro_seq_file,na.strings=c("NA","nd"),header=T,sep='\t')
time = c(0,0.5,1,2,4,8,12,24)/24

####Xist
xist_data = pro_seq_data[pro_seq_data$Genes == "Xist",]
expression = get_expression_pro_seq(xist_data)
xist_pro = plot_expression(time,expression,"Xist")
xist_pro = xist_pro + theme(legend.position="none")

####Tsix
tsix_data = pro_seq_data[pro_seq_data$Genes == "Tsix",]
expression = get_expression_pro_seq(tsix_data)
tsix_pro = plot_expression(time,expression,"Tsix")
legend = get_legend(tsix_pro)
tsix_pro = tsix_pro + theme(legend.position="none")

cairo_pdf("/project/lncrna/Xist/plots/additional_analysis/plots_allelic_expression_pro_seq.pdf",height = 2,width = 5)
grid.arrange(xist_pro, tsix_pro, legend, ncol=3, widths=c(2, 2, 0.5))
dev.off()

###################################################################################
#plot allele-specific expression for undifferentiated mRNA-seq data
###################################################################################

mrna_seq_file = "/project/lncrna/Xist/data/silencing_halftimes/raw_data/UndifferenciatedCells.txt"
mrna_seq_file = read.table(mrna_seq_file,na.strings=c("NA","nd"),header=T,sep='\t')
time = c(0,2,4,8,12,24)/24
first_replicate = c(1,4,5,8,9,12,13,16,17,20,21,24,25)
mrna_seq_file = mrna_seq_file[,first_replicate]

####Xist
xist_data = mrna_seq_file[mrna_seq_file$Genes == "Xist",]
expression = get_expression_mrna_seq_undiff(xist_data)
xist_mrna_und = plot_expression(time,expression,"Xist")
xist_mrna_und = xist_mrna_und + theme(legend.position="none")

####Tsix
tsix_data = mrna_seq_file[mrna_seq_file$Genes == "Tsix",]
expression = get_expression_mrna_seq_undiff(tsix_data)
tsix_mrna_und = plot_expression(time,expression,"Tsix")
legend = get_legend(tsix_mrna_und)
tsix_mrna_und = tsix_mrna_und + theme(legend.position="none")

cairo_pdf("/project/lncrna/Xist/plots/additional_analysis/plots_allelic_expression_mRNA_undifferentiated.pdf",height = 2,width = 5)
grid.arrange(xist_mrna_und, tsix_mrna_und, legend, ncol=3, widths=c(2, 2, 0.5))
dev.off()

###################################################################################
#plot allele-specific expression for differentiated mRNA-seq data
###################################################################################

mrna_seq_file = "/project/lncrna/Xist/data/silencing_halftimes/raw_data/DifferenciatedCells.txt"
mrna_seq_file = read.table(mrna_seq_file,na.strings=c("NA","nd"),header=T,sep='\t')
time = c(0,8,16,24,48)/24

####Xist
xist_data = mrna_seq_file[mrna_seq_file$Genes == "Xist",]
expression = get_expression_mrna_seq_diff(xist_data)
xist_mrna_diff = plot_expression(time,expression,"Xist")
xist_mrna_diff = xist_mrna_diff + theme(legend.position="none")

####Tsix
tsix_data = mrna_seq_file[mrna_seq_file$Genes == "Tsix",]
expression = get_expression_mrna_seq_diff(tsix_data)
tsix_mrna_diff = plot_expression(time,expression,"Tsix")
legend = get_legend(tsix_mrna_diff)
tsix_mrna_diff = tsix_mrna_diff + theme(legend.position="none")

cairo_pdf("/project/lncrna/Xist/plots/additional_analysis/plots_allelic_expression_mRNA_differentiated.pdf",height = 2,width = 5)
grid.arrange(xist_mrna_diff, tsix_mrna_diff, legend, ncol=3, widths=c(2, 2, 0.5))
dev.off()


###################################################################################
####plots for paper
###################################################################################

cairo_pdf("/project/lncrna/Xist/plots/additional_analysis/paper_figures_allelic_expression.pdf",height = 7,width = 4.5)
xist_pro = xist_pro + theme(legend.position = c(0.5,0.3))
grid.arrange(arrangeGrob(xist_pro, tsix_pro, top="undifferentiated PRO-seq",ncol=2), 
             arrangeGrob(xist_mrna_und, tsix_mrna_und, top="undifferentiated mRNA-seq",ncol=2), 
             arrangeGrob(xist_mrna_diff, tsix_mrna_diff, top = "differentiated mRNA-seq",ncol=2), ncol=1, widths=c(4))
dev.off()

cairo_pdf("/project/lncrna/Xist/plots/additional_analysis/paper_figures_allelic_expression_pro_seq.pdf",height = 3,width = 2)
xist_pro = xist_pro + theme(legend.position = "top") + guides(fill=guide_legend(nrow=3), col=guide_legend(nrow=3))
xist_pro
dev.off()
