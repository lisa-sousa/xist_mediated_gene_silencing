library(Hmisc)

pdf("/project/lncrna/Xist/plots_lisa/silencing_halftimes/sat1_wdr13_fits.pdf")

wdr13 = read.table("/project/lncrna/Xist/data_lisa/silencing_halftimes/raw_data/wdr13.txt",header=T)


GS_XX = data.frame(tp = c(0,0.5,1,2,4,8,24), f_b6_reads = wdr13$X129[1:7])
GS_XX = GS_XX/100


pseudo_count = 0.001
data_trans = ((1-GS_XX$f_b6_reads[1])/GS_XX$f_b6_reads[1]) * (GS_XX$f_b6_reads/(1-GS_XX$f_b6_reads+pseudo_count))
data_trans = data.frame(time = c(0,0.5,1,2,4,8,24)/24, expression = data_trans)


#fitting parameters
fitting_function <- function(t, start) {
  #ret=exp(-t/start)
  ret=exp(-start*t)
}

control=nls.control(warnOnly=TRUE) # prevent stopping the fit


#fitnl=nls(expression ~ fitting_function(time,k), data=data_trans, start=list(k=1),control=control, algorithm="port",upper=c(k=5))
fitnl=nls(expression ~ fitting_function(time,k), data=data_trans, start=list(k=1),control=control, algorithm="port",lower=c(k=1/5))

k = abs(coef(fitnl)[1])
#fitted_data$halftime[i] = k * log(2)
halftime = log(2)/k
RSS = sqrt(sum(resid(fitnl)^2))

plot(data_trans$time,data_trans$expression,xlim=c(-0.2, 1.2),ylim=c(0,1.2),xlab='time [days]', ylab='BL6',main="wdr13 GS XX")
lines(seq.int(0,1,0.01), fitting_function(seq.int(0,1,0.01),k),lwd=2)
legend("topright",legend = paste('square root RSS=',round(sqrt(sum(resid(fitnl)^2)),2),'\nhalftime=',round(log(2)/k,2)))


TX_XX = data.frame(tp = c(0,0.5,1,2,4,8), f_b6_reads = c(mean(wdr13$X129[c(8,14)]),mean(wdr13$X129[c(9,15,20)]),
                                                         mean(wdr13$X129[c(10,16,21)]),mean(wdr13$X129[c(11,17,22)]),
                                                         mean(wdr13$X129[c(12,18,23)]),mean(wdr13$X129[c(13,19,24)])))
TX_XX = TX_XX/100

pseudo_count = 0.001
data_trans = ((1-TX_XX$f_b6_reads[1])/TX_XX$f_b6_reads[1]) * (TX_XX$f_b6_reads/(1-TX_XX$f_b6_reads+pseudo_count))
data_trans = data.frame(time = c(0,0.5,1,2,4,8)/24, expression = data_trans)


#fitting parameters
fitting_function <- function(t, start) {
  #ret=exp(-t/start)
  ret=exp(-start*t)
}

control=nls.control(warnOnly=TRUE) # prevent stopping the fit


fitnl=nls(expression ~ fitting_function(time,k), data=data_trans, start=list(k=1),control=control, algorithm="port",lower=c(k=1/5))

k = abs(coef(fitnl)[1])
halftime = log(2)/k
RSS = sqrt(sum(resid(fitnl)^2))

plot(data_trans$time,data_trans$expression,xlim=c(-0.2, 1.2),ylim=c(0,1.2),xlab='time [days]', ylab='BL6',main="wdr13 TX XX")
lines(seq.int(0,1,0.01), fitting_function(seq.int(0,1,0.01),k),lwd=2)
legend("topright",legend = paste('square root RSS=',round(sqrt(sum(resid(fitnl)^2)),2),'\nhalftime=',round(log(2)/k,2)))

######plot with error bars

normalize_expression = function(data_set){
  data_set = data_set/100
  pseudo_count = 0.001
  data_trans = ((1-data_set$f_b6_reads[1])/data_set$f_b6_reads[1]) * (data_set$f_b6_reads/(1-data_set$f_b6_reads+pseudo_count))
  data_trans = data.frame(time = c(0,0.5,1,2,4,8)/24, expression = data_trans)
  return(data_trans)
}

TX_XX = data.frame(tp = c(0,0.5,1,2,4,8), f_b6_reads = c(mean(wdr13$X129[c(8,14)]),mean(wdr13$X129[c(9,15,20)]),
                                                         mean(wdr13$X129[c(10,16,21)]),mean(wdr13$X129[c(11,17,22)]),
                                                         mean(wdr13$X129[c(12,18,23)]),mean(wdr13$X129[c(13,19,24)])))
TX_XX = normalize_expression(TX_XX)

data_A = data.frame(tp = c(0,0.5,1,2,4,8), f_b6_reads = wdr13$X129[c(8:13)])
data_A = normalize_expression(data_A)

data_B = data.frame(tp = c(0,0.5,1,2,4,8), f_b6_reads = wdr13$X129[c(14:19)])
data_B = normalize_expression(data_B)

data_C = data.frame(tp = c(0,0.5,1,2,4,8), f_b6_reads = c(mean(wdr13$X129[c(14,8)]),wdr13$X129[20:24]))
data_C = normalize_expression(data_C)

table = data.frame(time=data_A$time, A=data_A$expression, B=data_B$expression, C=data_C$expression)

delta_lwr = min(table[1,2:4])
delta_lwr = c(delta_lwr,min(table[2,2:4]))
delta_lwr = c(delta_lwr,min(table[3,2:4]))
delta_lwr = c(delta_lwr,min(table[4,2:4]))
delta_lwr = c(delta_lwr,min(table[5,2:4]))
delta_lwr = c(delta_lwr,min(table[6,2:4]))

delta_upr = max(table[1,2:4])
delta_upr = c(delta_upr,max(table[2,2:4]))
delta_upr = c(delta_upr,max(table[3,2:4]))
delta_upr = c(delta_upr,max(table[4,2:4]))
delta_upr = c(delta_upr,max(table[5,2:4]))
delta_upr = c(delta_upr,max(table[6,2:4]))
errbar(TX_XX$time,TX_XX$expression, delta_lwr, delta_upr , ylim = c(0,1.4), xlim=c(0,0.5), xlab = "time",ylab = "B6 expression")


###############Sat1
sat1 = read.table("/project/lncrna/Xist/data_lisa/silencing_halftimes/raw_data/sat1.txt",header=T)


GS_XX = data.frame(tp = c(0,0.5,1,2,4,8,24), f_b6_reads = sat1$X129[1:7])
GS_XX = GS_XX/100


pseudo_count = 0.001
data_trans = ((1-GS_XX$f_b6_reads[1])/GS_XX$f_b6_reads[1]) * (GS_XX$f_b6_reads/(1-GS_XX$f_b6_reads+pseudo_count))
data_trans = data.frame(time = c(0,0.5,1,2,4,8,24)/24, expression = data_trans)


#fitting parameters
fitting_function <- function(t, start) {
  #ret=exp(-t/start)
  ret=exp(-start*t)
}

control=nls.control(warnOnly=TRUE) # prevent stopping the fit


#fitnl=nls(expression ~ fitting_function(time,k), data=data_trans, start=list(k=1),control=control, algorithm="port",upper=c(k=5))
fitnl=nls(expression ~ fitting_function(time,k), data=data_trans, start=list(k=1),control=control, algorithm="port",lower=c(k=1/5))

k = abs(coef(fitnl)[1])
#fitted_data$halftime[i] = k * log(2)
halftime = log(2)/k
RSS = sqrt(sum(resid(fitnl)^2))

plot(data_trans$time,data_trans$expression,xlim=c(-0.2, 1.2),ylim=c(0,1.2),xlab='time [days]', ylab='BL6',main="sat1 GS XX")
lines(seq.int(0,1,0.01), fitting_function(seq.int(0,1,0.01),k),lwd=2)
legend("topright",legend = paste('square root RSS=',round(sqrt(sum(resid(fitnl)^2)),2),'\nhalftime=',round(log(2)/k,2)))


TX_XX = data.frame(tp = c(0,0.5,1,2,4,8), f_b6_reads = c(mean(sat1$X129[c(8,14)]),mean(sat1$X129[c(9,15,20)]),
                                                         mean(sat1$X129[c(10,16,21)]),mean(sat1$X129[c(11,17,22)]),
                                                         mean(sat1$X129[c(12,18,23)]),mean(sat1$X129[c(13,19,24)])))
TX_XX = TX_XX/100

pseudo_count = 0.001
data_trans = ((1-TX_XX$f_b6_reads[1])/TX_XX$f_b6_reads[1]) * (TX_XX$f_b6_reads/(1-TX_XX$f_b6_reads+pseudo_count))
data_trans = data.frame(time = c(0,0.5,1,2,4,8)/24, expression = data_trans)


#fitting parameters
fitting_function <- function(t, start) {
  #ret=exp(-t/start)
  ret=exp(-start*t)
}

control=nls.control(warnOnly=TRUE) # prevent stopping the fit


#fitnl=nls(expression ~ fitting_function(time,k), data=data_trans, start=list(k=1),control=control, algorithm="port",upper=c(k=5))
fitnl=nls(expression ~ fitting_function(time,k), data=data_trans, start=list(k=1),control=control, algorithm="port",lower=c(k=1/5))

k = abs(coef(fitnl)[1])
#fitted_data$halftime[i] = k * log(2)
halftime = log(2)/k
RSS = sqrt(sum(resid(fitnl)^2))

plot(data_trans$time,data_trans$expression,xlim=c(-0.2, 1.2),ylim=c(0,1.2),xlab='time [days]', ylab='BL6',main="sat1 TX XX")
lines(seq.int(0,1,0.01), fitting_function(seq.int(0,1,0.01),k),lwd=2)
legend("topright",legend = paste('square root RSS=',round(sqrt(sum(resid(fitnl)^2)),2),'\nhalftime=',round(log(2)/k,2)))



######error bar plots
normalize_expression = function(data_set){
  data_set = data_set/100
  pseudo_count = 0.001
  data_trans = ((1-data_set$f_b6_reads[1])/data_set$f_b6_reads[1]) * (data_set$f_b6_reads/(1-data_set$f_b6_reads+pseudo_count))
  data_trans = data.frame(time = c(0,0.5,1,2,4,8)/24, expression = data_trans)
  return(data_trans)
}


TX_XX = data.frame(tp = c(0,0.5,1,2,4,8), f_b6_reads = c(mean(sat1$X129[c(8,14)]),mean(sat1$X129[c(9,15,20)]),
                                                         mean(sat1$X129[c(10,16,21)]),mean(sat1$X129[c(11,17,22)]),
                                                         mean(sat1$X129[c(12,18,23)]),mean(sat1$X129[c(13,19,24)])))
TX_XX = normalize_expression(TX_XX)

data_A = data.frame(tp = c(0,0.5,1,2,4,8), f_b6_reads = sat1$X129[c(8:13)])
data_A = normalize_expression(data_A)

data_B = data.frame(tp = c(0,0.5,1,2,4,8), f_b6_reads = sat1$X129[c(14:19)])
data_B = normalize_expression(data_B)

data_C = data.frame(tp = c(0,0.5,1,2,4,8), f_b6_reads = c(mean(sat1$X129[c(14,8)]),sat1$X129[20:24]))
data_C = normalize_expression(data_C)

table = data.frame(time=data_A$time, A=data_A$expression, B=data_B$expression, C=data_C$expression)

delta_lwr = min(table[1,2:4])
delta_lwr = c(delta_lwr,min(table[2,2:4]))
delta_lwr = c(delta_lwr,min(table[3,2:4]))
delta_lwr = c(delta_lwr,min(table[4,2:4]))
delta_lwr = c(delta_lwr,min(table[5,2:4]))
delta_lwr = c(delta_lwr,min(table[6,2:4]))

delta_upr = max(table[1,2:4])
delta_upr = c(delta_upr,max(table[2,2:4]))
delta_upr = c(delta_upr,max(table[3,2:4]))
delta_upr = c(delta_upr,max(table[4,2:4]))
delta_upr = c(delta_upr,max(table[5,2:4]))
delta_upr = c(delta_upr,max(table[6,2:4]))
errbar(TX_XX$time,TX_XX$expression, delta_lwr, delta_upr , ylim = c(0,1.4), xlim=c(0,0.5), xlab = "time",ylab = "B6 expression")

dev.off()

