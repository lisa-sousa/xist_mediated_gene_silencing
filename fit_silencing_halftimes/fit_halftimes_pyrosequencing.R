library(ggplot2)
library(gridExtra)
library(tidyr)
library(plyr)
library(dplyr)

control=nls.control(warnOnly=TRUE) # prevent stopping the fit
not.silenced.genes = c('Maoa','Kcnd1','Pim2','Pqbp1','Wdr13')
silenced.genes = c('Taf9b', 'Rnf128', 'Ripply1', 'Foxo4', 'Sat1','Nxt2')
all.genes = c(silenced.genes,not.silenced.genes)

# read in data
data = read.table('/project/lncrna/Xist/data/silencing_halftimes/raw_data/pyrosequencing.txt',header=T,dec=',')

# get mean fraction B6 at t=0
t0 = data %>% filter(time==0) %>% group_by(gene) %>% summarize(t0=mean(b6)) %>% ungroup()

# normalize to t=0
data.norm = data %>% left_join(t0) %>% mutate(norm=b6/(100-b6)*(100-t0)/t0)
data.norm$class = "none"
data.norm$class[data.norm$gene %in% silenced.genes] = "silenced"
data.norm$class[data.norm$gene %in% not.silenced.genes] = "not silenced"

# define the function used to fit the data
fitting_function <- function(t, start) {
  return(exp(-start*t))
}

# Fit each gene with an exponetial and plot
halftime = NULL
expression_all = NULL
fit = NULL
label = c()

for (g in c(1:length(all.genes))) {
  expression = data.norm %>% filter(gene==all.genes[g]) %>% mutate(time_d=time/24) %>% select(time_d,norm,gene,class)
  fitnl=nls(norm ~ fitting_function(time_d,k), data=expression, start=list(k=1),control=control, algorithm='port',lower=c(k=1/5))
  k = abs(coef(fitnl)[1])
  halftime = rbind(halftime,data.frame(gene = all.genes[g], halftime = log(2)/k))
  
  expression_all = rbind(expression_all,expression)
  fit = rbind(fit,data.frame(time_d = seq.int(0,1,0.01),fitted_data = fitting_function(seq.int(0,1,0.01),k), gene = all.genes[g]))
  label[g] = deparse(bquote(t[1/2] == .(as.double(round(log(2)/k,2)))~d))
}

expression_all$gene = factor(expression_all$gene,levels = c(silenced.genes,not.silenced.genes),ordered = TRUE)
halftime_label = data.frame(gene = all.genes, x= 0.1,y= 0.1,label = label)

cairo_pdf(file='/project/lncrna/Xist/plots/silencing_halftimes/paper_figures_pyro_genes.pdf',height=3, width=6)
ggplot(expression_all, aes(x=time_d,y=norm)) +
  geom_point(size=0.6)+
  geom_line(data=fit,aes(y=fitted_data)) +
  facet_wrap(.~gene,labeller = label_wrap_gen(width = 20, multi_line = TRUE),nrow=2) +
  geom_text(data=halftime_label,aes(x,y,label=label),parse=T,size=2.5,hjust=0) +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),axis.text=element_text(size=6), axis.title=element_text(size=8), 
        strip.text=element_text(size=8,hjust = 0.5,face="italic")) + 
  scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1), name='time [days]')  + 
  scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.5,1), name='norm. B6 expression') 
dev.off()

# Test difference in halftime of not silenced and silenced predicted genes
wil = wilcox.test(halftime$halftime[halftime$gene %in% not.silenced.genes],halftime$halftime[halftime$gene %in% silenced.genes])

# Plot halftime of not silenced and silenced predicted genes as a dotplot

halftime$class = "none"
halftime$class[halftime$gene %in% silenced.genes] = "silenced"
halftime$class[halftime$gene %in% not.silenced.genes] = "not silenced"
halftime$class = factor(halftime$class,levels = c("silenced","not silenced"),ordered = TRUE)

cairo_pdf(file='/project/lncrna/Xist/plots/silencing_halftimes/paper_figures_pyro_statistics.pdf',height=3, width=2)
ggplot(halftime,aes(x=class, y=halftime)) + 
  geom_jitter(position=position_jitter(0.2), size=0.6) + 
  stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean,geom="crossbar", color='grey', width=0.7, size=0.2)+
  geom_text(aes(x=1.5,y=4.5), label=paste('p=',round(wil$p.value,4),sep=""), size=2.5) +
  geom_segment(aes(x=1,xend=2,y=4,yend=4)) +
  theme_minimal(base_family = "Source Sans Pro") + 
  theme(panel.grid.minor = element_blank(), axis.text=element_text(size=8), axis.title=element_text(size=8)) + 
  scale_y_continuous(limits=c(0,5),breaks=c(0,2,4,6), name='half-time [days]') + 
  scale_x_discrete (name='predicted class')
dev.off()



