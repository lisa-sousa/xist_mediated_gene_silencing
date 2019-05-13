###################################################################################
#libraries
###################################################################################

library(Cairo)
library(gplots)
library(corrplot)

###################################################################################
#data and directories
###################################################################################

output_directory = '/project/lncrna/Xist/plots/additional_analysis/'

load('/project/lncrna/Xist/data/modelling/feature_matrix/promoter_matrix_reannotated_normRAdjusted_pro_seq_genes.RData')

for(i in 1:ncol(data_set)){
  if(is.factor(data_set[,i])){
    data_set[,i] = as.numeric(data_set[,i])
  }
}

full_data_set = cbind.data.frame(data_set,halftime)
data_matrix = data.matrix(full_data_set)

colnames(data_matrix) = gsub("HiC ","Hi-C ",gsub("number interactions","number",gsub("mean interaction ","",
                      gsub("_"," ",gsub("8WG16","unphosphorylated",gsub('_GLIB_*-?[0-9]*_-?[0-9]*','',gsub('_ENCODE_[A-z]{4}_*-?[0-9]*_-?[0-9]*','',
                      gsub('_GSE[0-9]*_-?[0-9,A-z]*_-?[0-9,A-z]*','',colnames(data_matrix),perl=T))))))))
###################################################################################
#correlation plots
###################################################################################

cairo_pdf(file=paste(output_directory,"plots_correlation.pdf",sep="/"), width = 16, height = 16, onefile = T)

#correlation matrix with correlation values printed in each box
corResults = cor(data_matrix)
corResults = round(corResults,2)

n = ncol(corResults)
m = nrow(corResults)

a = seq(0,1,length=n)
x = c()
for(i in 1:n){x = c(x,rep(a[i],m))}

b = seq(0,1,length=m)
y = rep(b,n)

par(mfrow=c(1,1),mar=c(15,15,4,1),oma=c(2,2,0,0))
image(corResults,xaxt="n",yaxt="n",col=bluered(20),main="correlation between features")
axis(1,at=b,labels=rownames(corResults),las=2)
axis(2,at=a,labels=colnames(corResults),las=2)
text(y,x,round(unlist(as.list(corResults)),2),cex=0.4)


#correlation matrix with sorted entries (based on correlation)
par(mar=c(40,4,4,2)) 
heatmap.2(cor(data_matrix), cexRow=0.8, cexCol=0.8, col="bluered", trace="none", margins=c(22,20), dendrogram="none", key=TRUE, symkey=FALSE, density.info="none")

#correlation matrix with sorted entries (based on correlation) and points in size of correlation
par(mar=c(10,10,4,2),oma=c(2,2,0,0)) 
corrplot(cor(data_matrix),order='hclust', col=colorRampPalette(c('blue','white','red'))(200),diag=F,cl.cex=1.7,tl.cex=1,tl.col = "black")

dev.off()

