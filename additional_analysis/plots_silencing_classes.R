###################################################################################
#libraries
###################################################################################

library(Cairo)
library(ggplot2)

###################################################################################
#load data
###################################################################################

table_halftimes = read.table("/project/lncrna/Xist/data/silencing_halftimes/fitted_data/halftimes_pro_seq_mm9_reannotated_with_rr.bed")
halftimes = table_halftimes$V5

early = halftimes[halftimes < 0.5]
silenced = halftimes[halftimes < 0.9]
late = halftimes[halftimes > 0.9 & halftimes < 1.3]
not_silenced = halftimes[halftimes > 1.6]

data_plot = data.frame(halftime = early, silencing_class = rep("early",length(early)))
data_plot = rbind(data_plot, data.frame(halftime = late, silencing_class = rep("late",length(late))))
data_plot = rbind(data_plot, data.frame(halftime = silenced, silencing_class = rep("silenced",length(silenced))))
data_plot = rbind(data_plot, data.frame(halftime = not_silenced, silencing_class = rep("not_silenced",length(not_silenced))))

###################################################################################
#boxplot of different silncing classes
###################################################################################

CairoPDF(file = "/project/lncrna/Xist/plots/additional_analysis/boxplots_silencing_classes.pdf", width = 15, height = 15)
par(mfrow=c(1,1),mar=c(15,10,10,5),oma=c(5,5,5,5))

feature = 'half-time'

gg_box = ggplot(data_plot, aes(x=silencing_class,y=halftime)) + geom_boxplot(notch=FALSE,fill = "lightgrey", colour = "black",alpha = 0.7,outlier.shape = 16) + ggtitle(feature) + theme_bw() + 
  theme(axis.text.x=element_text(hjust = 0.5,size=20),axis.text.y=element_text(size=20),axis.title = element_text(face="bold", size=20),
        plot.title = element_text(hjust = 0.5,size=35,face='bold'),plot.margin = unit(c(3,3,3,3), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
  scale_x_discrete(name = "silencing class") + scale_y_continuous(name = "half-time")
print(gg_box)

dev.off()
