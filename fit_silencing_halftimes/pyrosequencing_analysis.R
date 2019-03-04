
library('ggplot2')
library('gridExtra')
library('tidyr')
library('plyr')
library('dplyr')
control=nls.control(warnOnly=TRUE) # prevent stopping the fit
slow.genes = c('Maoa','Kcnd1','Pim2','Pqbp1','Wdr13')
fast.genes = c('Taf9b', 'Rnf128', 'Ripply1', 'Foxo4', 'Sat1','Nxt2')
# read in data
data = read.table('/project/lncrna/Xist/data/silencing_halftimes/raw_data/pyrosequencing.txt',header=T,dec=',')

# get mean fraction B6 at t=0
t0 = data %>% filter(time==0) %>% group_by(gene) %>% summarize(t0=mean(b6)) %>% ungroup()

# normalize to t=0
data.norm = data %>% left_join(t0) %>% mutate(norm=b6/(100-b6)*(100-t0)/t0)

# define the function used to fit the data
func <- function(t, start) {
  ret=exp(-t/start)
}

plot.genes = c(fast.genes, slow.genes)
halftime=data.frame(gene=plot.genes,k=NA);

plot.list=list()
# Fit each gene with an exponetial and plot
for (g in c(1:length(plot.genes))) {
  expression = data.norm %>% filter(gene==plot.genes[g]) %>% mutate(time_d=time/24) %>% select(time_d,norm)
  fitnl=nls(norm ~ exp(-time_d/k), data=expression, start=list(k=0.5),control=control, algorithm = "port",upper=c(k=5))
  k=abs(coef(fitnl)[1])
  halftime$k[g]=log(2)*k
  p1  = expression %>%  ggplot(aes(x=time_d, y=norm)) + geom_point(size=0.6) + 
    stat_function(fun=func, args = list(k), size=0.3) + 
    geom_text(aes(x=0.5,y=0.1),label=paste('half-time=',round(log(2)*k,2),'h'),size=2) +
    theme_minimal() + theme(text=element_text(size=8), plot.title=element_text(size=8), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) + 
    scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1), name='Time [days]')  + 
    scale_y_continuous(limits=c(0,1.2),breaks=c(0,0.5,1), name='norm. B6 Fraction') + ggtitle(plot.genes[g])
  plot.list[[g]]=p1
}

# Test difference in halftime of slow and fast predicted genes
wil = wilcox.test(halftime$k[halftime$gene %in% slow.genes],halftime$k[halftime$gene %in% fast.genes])

# Plot halftime of slow and fast predicted genes as a dotplot
set.seed(2)
p1 = halftime %>% mutate(pred=ifelse(gene %in% fast.genes,'fast','slow')) %>% 
  ggplot(aes(x=pred, y=k)) + 
  geom_jitter(position=position_jitter(0.2), size=0.6) + 
  stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean,geom="crossbar", color='grey', width=0.7, size=0.2)+
  geom_text(aes(x=1.5,y=4.5), label=paste('p=',round(wil$p.value,4),sep=""), size=2) +
  geom_segment(aes(x=1,xend=2,y=4,yend=4)) +
  theme_minimal() + theme(text=element_text(size=8), panel.grid.minor = element_blank())+ scale_y_continuous(limits=c(0,5),breaks=c(0,2,4,6), name='Halftime [days]') + 
  scale_x_discrete (name='Predicted Class')
plot.list[[g+1]]=p1

pdf(file='/project/lncrna/Xist/plots/silencing_halftimes/pyrosequencing.pdf',height=3, width=7)
grid.arrange(grobs=plot.list, ncol=6,nrow=2)
dev.off()



