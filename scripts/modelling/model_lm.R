########################################
#libraries
########################################
library(glmnet)
library(foreach)
library(doParallel)
library(energy)
library(Cairo)
library(gplots)
library(here)
library(matrixStats)
library(ggplot2)
library(reshape2)
source(here("scripts/modelling","model_functions_lm.R"))

########################################
#directories and files
########################################

output_directory_plots = 'plots/modelling/xci_escape_linear_model'

file_feature_matrix = here("data/modelling/feature_matrix","promoter_matrix_reannotated_normRAdjusted_pro_seq_genes.RData")
file_high_confidence_escapees = here("data/annotation_files/escapees","escapees.RData")
file_results = file(here(output_directory_plots,'log_results.txt'), open = "a")

########################################
#setup parameters
########################################

#threshold combination for "silencing dynamics model"
thr_silencing_lower_seq = c(0.5,0.6,0.7)
thr_silencing_middle_seq = c(0.7,0.8,0.9,1)
thr_silencing_upper_seq = c(1,1.1,1.2,1.3,1.4)
class0_label = "early silenced genes"
class1_label = "late silenced genes"

#threshold combination for "XCI/escape" model
thr_silencing_lower_seq = c(0.9,1.0,1.1,1.2,1.3,1.4)
thr_silencing_middle_seq = "-"
thr_silencing_upper_seq = c(1.4,1.5,1.6,1.7,1.8,1.9,2.0)
class0_label = "silenced genes"
class1_label = "not silenced genes"

ncores = 20 # number of cores for parallelization
test_set_min_size = 10 # minimum number of genes in OOB set
runs = 200 # number of runs for stability test of random forest model
runs_lambda = 50
alpha = 0.5 #1:Lasso, 0.5:ElasticNet, 0:Ridge

########################################
#cluster setup for parallel computing
########################################

comb <- function(...) {
  mapply('cbind', ..., SIMPLIFY=FALSE)
}

my_cluster = makeCluster(ncores)
registerDoParallel(my_cluster)

########################################
#run linear regression model
########################################

load(file_feature_matrix)
data_set$overlap_with_CpG_islands[as.numeric(data_set$overlap_with_CpG_islands) > 2] = 1
cairo_pdf(file = paste(output_directory_plots,"/linear_model.pdf",sep=''), width = 9, height = 5,onefile = T)
linear_regression(data_set,halftime,runs)
dev.off()

########################################
#run pipeline on all combination of thresholds
########################################

for(thr_silencing_lower in thr_silencing_lower_seq){
  for(thr_silencing_middle in thr_silencing_middle_seq){
    for(thr_silencing_upper in thr_silencing_upper_seq){
      
      ########################################
      #load data
      ########################################
      
      data = load_data(file_feature_matrix, thr_silencing_lower, thr_silencing_middle, thr_silencing_upper, file_high_confidence_escapees)
      data[[1]]$overlap_with_CpG_islands[as.numeric(data[[1]]$overlap_with_CpG_islands) > 2] = 1
      
      smaller_set = min(c(sum(data[[2]]==0),sum(data[[2]]==1)))
      training = round(min(c(smaller_set/100*80,smaller_set-test_set_min_size)))
      test = smaller_set - training
      
      if(training >= 25){
        ########################################
        #setup for new threshold
        ########################################
        
        print(paste('thresholds:',thr_silencing_lower,thr_silencing_middle,thr_silencing_upper))
        
        #new entry in log file
        cat("#########################\n#########################\n#########################", file = file_results, sep = "\n")
        cat(paste('thresholds:',thr_silencing_lower,thr_silencing_middle,thr_silencing_upper), file = file_results, sep = "\n")
        cat(paste('size class 0:',sum(data[[2]]==0)), file = file_results, sep = "\n")
        cat(paste('size class 1:',sum(data[[2]]==1)), file = file_results, sep = "\n")
        cat(paste('training set:',training), file = file_results, sep = "\n")
        
        ########################################
        #train linear model on all features
        ########################################
        
        cairo_pdf(file = paste(output_directory_plots,"/",thr_silencing_lower,"_",thr_silencing_middle,"_",thr_silencing_upper,".pdf",sep=''), width = 9, height = 5,onefile = T)
        
        error_log_reg = logistic_regression(data,training,runs)
        cat("#########################\nlogistic regression", file = file_results, sep = "\n")
        cat(paste('total error (sd):',round(mean(error_log_reg)*100,2), "(", round(sd(error_log_reg)*100,2), ")"), file = file_results, sep = "\n")
        
        cat("#########################\nregularized logistic regression", file = file_results, sep = "\n")
        error_reg_log_reg = regularized_logistic_regression(data,training,alpha,runs,runs_lambda,class0_label,class1_label)
        
        cat(paste('total error (sd):',round(mean(error_reg_log_reg[,1])*100,2), "(", round(sd(error_reg_log_reg[,1])*100,2), ")", 
                  '; class 0 error (sd):',round(mean(error_reg_log_reg[,2])*100,2), "(", round(sd(error_reg_log_reg[,2])*100,2), ")",
                  '; class 1 error (sd):',round(mean(error_reg_log_reg[,3])*100,2), "(", round(sd(error_reg_log_reg[,3])*100,2), ")"), 
            file = file_results, sep = "\n")
        dev.off()
      }
    }
  }
}

########################################
#close cluster and file
########################################

close(file_results)
stopCluster(my_cluster)

