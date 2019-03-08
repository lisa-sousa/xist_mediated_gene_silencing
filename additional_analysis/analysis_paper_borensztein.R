###################################################################################
#libraries
###################################################################################

library(Cairo)
library(ggplot2)

###################################################################################
#load data
###################################################################################

load('/project/lncrna/Xist/data/modelling/feature_matrix/promoter_matrix_reannotated_normRAdjusted_pro_seq_genes.RData')
table_halftimes = data.frame(gene = rownames(data_set), halftime = halftime)

table_NSMB_paper = read.table(file = '/project/lncrna/Xist/data_lisa/annotation_files/imprinted_xci_silencing_rates/consistent_genes_categories.txt',sep='\t',header = F)
colnames(table_NSMB_paper) = c('gene','silencing_class')

table = merge(table_halftimes,table_NSMB_paper,by='gene')

###################################################################################
#plot boxplots
###################################################################################

CairoPDF(file = "/project/lncrna/Xist/plots/additional_analysis/analysis_paper_borensztein.pdf", width = 15, height = 15)
par(mfrow=c(1,1),mar=c(15,10,10,5),oma=c(5,5,5,5))

feature = 'half-time'
data_plot = table[,c(2,3)]
data_plot$silencing_class = factor(data_plot$silencing_class,levels = c('Early','Inter','Late','Esc','Bias'),ordered = TRUE)

gg_box = ggplot(data_plot, aes(x=silencing_class,y=halftime)) + geom_boxplot(notch=FALSE,fill = "lightgrey", colour = "black",alpha = 0.7,outlier.shape = 16) + ggtitle(feature) + theme_bw() + 
    theme(axis.text.x=element_text(hjust = 0.5,size=20),axis.text.y=element_text(size=20),axis.title = element_text(face="bold", size=20),
          plot.title = element_text(hjust = 0.5,size=35,face='bold'),plot.margin = unit(c(3,3,3,3), "cm"), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
    scale_x_discrete(name = "silencing class") + scale_y_continuous(name = "half-time")
print(gg_box)

dev.off()



