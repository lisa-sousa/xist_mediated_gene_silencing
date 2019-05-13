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

table_NSMB_paper = read.table(file = '/project/lncrna/Xist/data/annotation_files/imprinted_xci_silencing_rates/consistent_genes_categories.txt',sep='\t',header = F)
colnames(table_NSMB_paper) = c('gene','silencing_class')

table = merge(table_halftimes,table_NSMB_paper,by='gene')
table = table[table$silencing_class != "Bias",2:3]
table$silencing_class = factor(table$silencing_class,levels = c('Early','Inter','Late','Esc'),ordered = TRUE)

###################################################################################
#plot boxplots
###################################################################################

cairo_pdf("/project/lncrna/Xist/plots/additional_analysis/analysis_paper_borensztein.pdf",width = 2,height = 2.3, onefile = TRUE)
ggplot(table, aes(x=silencing_class,y=halftime)) + 
  geom_boxplot(colour = "#4d4d4d",alpha = 0.7,outlier.size=0.1,lwd=0.4) + 
  ggtitle("Pre-implantation embryos") + 
  scale_x_discrete(name = "silencing class",labels=c("early","interm.","late","escapee")) + 
  scale_y_continuous(breaks=c(0,1,2,3,3.5), label=c("0","1","2","3",">3.5"), name='half-time [days]') +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text.x = element_text(size=8, angle = 45, hjust=1, margin = margin(t=0,b=0)), axis.text.y = element_text(size=8), 
        axis.title=element_text(size=8, margin = margin(t=0)),plot.title = element_text(size=7)) 
dev.off()
