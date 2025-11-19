library(MendelianRandomization)
library(TwoSampleMR)
library(mr.raps)
library(MRPRESSO)
library(ggplot2)
library(dplyr)
library(R2jags)

source("D:/PhD Research/uch/research/autoMR/source/MR_check.R")
source("D:/PhD Research/uch/research/autoMR/source/MR_plots.R")
source("D:/PhD Research/uch/research/autoMR/source/MR_forest.R")
source("D:/PhD Research/uch/research/autoMR/source/MR_validoutput.R")
source("D:/PhD Research/uch/research/autoMR/source/GWAS_forest.R")
source("D:/PhD Research/uch/research/autoMR/source/mr_horse.R")

setwd("D:/PhD Research/uch/research/GLP1R MR/data/")

### data loading

heart_atrial_appendage_raw <- read.csv(file = "mr_input_glp1r_heart_atrial_appendage.csv")
heart_left_ventricle_raw <- read.csv(file = "mr_input_glp1r_heart_left_ventricle.csv")
pancreas_raw <- read.csv(file = "mr_input_glp1r_pancreas.csv")

heart_atrial_appendage_raw$Exposure = 'heart_atrial_appendage'
heart_left_ventricle_raw$Exposure = 'heart_left_ventricle'
pancreas_raw$Exposure = 'pancreas'

### rename columns

heart_atrial_appendage_raw <- heart_atrial_appendage_raw %>%
  rename(
    NEA_outcome = NEA_Outcome,	
    EA_outcome = EA_Outcome,
    A1FREQ_outcome = EAFreq,
    beta_outcome = Beta_Outcome,
    se_outcome = SE_Outcome,
    p_outcome	= P_Outcome,
    NEA_exposure = NEA_Exposure,
    EA_exposure	= EA_Exposure,
    A1FREQ_exposure	= EAFreq_Outcome,
    beta_exposure	= Beta_Exposure,
    se_exposure = SE_Exposure
  )

heart_left_ventricle_raw <- heart_left_ventricle_raw %>%
  rename(
    NEA_outcome = NEA_Outcome,	
    EA_outcome = EA_Outcome,
    A1FREQ_outcome = EAFreq,
    beta_outcome = Beta_Outcome,
    se_outcome = SE_Outcome,
    p_outcome	= P_Outcome,
    NEA_exposure = NEA_Exposure,
    EA_exposure	= EA_Exposure,
    A1FREQ_exposure	= EAFreq_Outcome,
    beta_exposure	= Beta_Exposure,
    se_exposure = SE_Exposure
  )

pancreas_raw <- pancreas_raw %>%
  rename(
    NEA_outcome = NEA_Outcome,	
    EA_outcome = EA_Outcome,
    A1FREQ_outcome = EAFreq,
    beta_outcome = Beta_Outcome,
    se_outcome = SE_Outcome,
    p_outcome	= P_Outcome,
    NEA_exposure = NEA_Exposure,
    EA_exposure	= EA_Exposure,
    A1FREQ_exposure	= EAFreq_Outcome,
    beta_exposure	= Beta_Exposure,
    se_exposure = SE_Exposure
  )

# comment the following lines for whole analysis
cols_select <- c("glucose","bmi","hba1c")
heart_atrial_appendage_raw <- subset(heart_atrial_appendage_raw, Outcome %in% cols_select)
heart_left_ventricle_raw <- subset(heart_left_ventricle_raw, Outcome %in% cols_select)
pancreas_raw <- subset(pancreas_raw, Outcome %in% cols_select)

### heart_atrial_appendage

# setwd("D:/PhD Research/uch/research/GLP1R MR/heart_atrial_appendage/")

setwd("D:/PhD Research/uch/research/GLP1R MR/testing/")

heart_atrial_appendage_MR_input <- process_snp_data(heart_atrial_appendage_raw)$input_df

heart_atrial_appendage <- valid.multi.output(MR_input_data = heart_atrial_appendage_MR_input,
                                                    outcome.form = NULL,
                                                    use_ivw = TRUE, 
                                                    use_raps = TRUE, 
                                                    use_median = TRUE, 
                                                    use_egger = TRUE, 
                                                    use_mr_presso = TRUE,
                                                    use_mr_horse = TRUE,
                                                    use_mr_grip  = TRUE,        
                                                    NbDistribution    = 1000,
                                                    SignifThreshold   = 0.05,
                                                    save_all_in_one   = TRUE,
                                                    save_csv          = TRUE,
                                                    mr_horse_n_iter   = 5000,
                                                    mr_horse_n_burnin = 1000,
                                                    mr_grip_parameters = NULL,
                                                    combined_file = "heart atrial appendage.csv")

MRplots.multi(MR_input_data = heart_atrial_appendage_MR_input,
              methods.plot          = c('IVW','RAPS','Egger','PRESSO','Horse','GRIP'),
              NbDistribution_presso = 1000,
              SignifThreshold_presso= 0.05,
              mr_horse_n_iter       = 5000,
              mr_horse_n_burnin     = 1000,
              plot.xlab = "Per allele assoc. with",
              plot.ylab = "Per allele assoc. with",
              show.legend = TRUE)

# MRplots.multi(MR_input_data = heart_atrial_appendage_MR_input,
#               methods.plot          = c('IVW','RAPS','Egger','PRESSO','Horse','GRIP'),
#               NbDistribution_presso = 1000,
#               SignifThreshold_presso= 0.05,
#               mr_horse_n_iter       = 5000,
#               mr_horse_n_burnin     = 1000,
#               plot.xlab = "Per allele assoc. with",
#               plot.ylab = "Per allele assoc. with",
#               show.legend = FALSE)

MR_forest(heart_atrial_appendage,effect = "Beta", 
          # custom_xlim = c(-0.01,0.01), 
          save_plot = TRUE, file_type = "jpeg",
          dot_size = 8,
          axis_text_size = 20,
          axis_title_size = 20,
          pval_text_size = 8,
          plot_width = 20,
          plot_height = 8)

MR_forest_SNP(heart_atrial_appendage_MR_input, report_form = "Beta", 
              # xlim_custom = c(-0.02,0.02),
              save_plot = TRUE, file_type = "jpeg",
              dot_size = 8,
              axis_text_size = 20,
              axis_title_size = 20,
              digits = 3,
              plot_width = 20,
              plot_height = 10,
              label_text_size = 8)


### heart_left_ventricle

setwd("D:/PhD Research/uch/research/GLP1R MR/heart_left_ventricle/")

heart_left_ventricle_MR_input <- process_snp_data(heart_left_ventricle_raw)$input_df

heart_left_ventricle <- valid.multi.output(MR_input_data = heart_left_ventricle_MR_input,
                                           outcome.form = NULL,
                                           use_ivw = TRUE, 
                                           use_raps = TRUE, 
                                           use_median = TRUE, 
                                           use_egger = TRUE, 
                                           use_mr_presso = TRUE,
                                           use_mr_horse = TRUE,
                                           use_mr_grip  = TRUE,        
                                           NbDistribution    = 1000,
                                           SignifThreshold   = 0.05,
                                           save_all_in_one   = TRUE,
                                           save_csv          = TRUE,
                                           mr_horse_n_iter   = 5000,
                                           mr_horse_n_burnin = 1000,
                                           mr_grip_parameters = NULL,
                                           combined_file = "heart left ventricle.csv")

MRplots.multi(MR_input_data = heart_left_ventricle_MR_input,
              methods.plot          = c('IVW','RAPS','Egger','PRESSO','Horse','GRIP'),
              NbDistribution_presso = 1000,
              SignifThreshold_presso= 0.05,
              mr_horse_n_iter       = 5000,
              mr_horse_n_burnin     = 1000,
              plot.xlab = "Per allele assoc. with",
              plot.ylab = "Per allele assoc. with",
              show.legend = TRUE)

# MRplots.multi(MR_input_data = heart_left_ventricle_MR_input,
#               methods.plot          = c('IVW','RAPS','Egger','PRESSO','Horse','GRIP'),
#               NbDistribution_presso = 1000,
#               SignifThreshold_presso= 0.05,
#               mr_horse_n_iter       = 5000,
#               mr_horse_n_burnin     = 1000,
#               plot.xlab = "Per allele assoc. with",
#               plot.ylab = "Per allele assoc. with",
#               show.legend = FALSE)

MR_forest(heart_left_ventricle,effect = "Beta", 
          # custom_xlim = c(-0.01,0.01), 
          save_plot = TRUE, file_type = "jpeg",
          dot_size = 8,
          axis_text_size = 20,
          axis_title_size = 20,
          pval_text_size = 8,
          plot_width = 20,
          plot_height = 8)

MR_forest_SNP(heart_left_ventricle_MR_input, report_form = "Beta", 
              # xlim_custom = c(-0.02,0.02),
              save_plot = TRUE, file_type = "jpeg",
              dot_size = 8,
              axis_text_size = 20,
              axis_title_size = 20,
              digits = 3,
              plot_width = 20,
              plot_height = 10,
              label_text_size = 8)

### pancreas

setwd("D:/PhD Research/uch/research/GLP1R MR/pancreas/")

pancreas_MR_input <- process_snp_data(pancreas_raw)$input_df

pancreas <- valid.multi.output(MR_input_data = pancreas_MR_input,
                               outcome.form = NULL,
                               use_ivw = TRUE, 
                               use_raps = TRUE, 
                               use_median = TRUE, 
                               use_egger = TRUE, 
                               use_mr_presso = TRUE,
                               use_mr_horse = TRUE,
                               use_mr_grip  = TRUE,        
                               NbDistribution    = 1000,
                               SignifThreshold   = 0.05,
                               save_all_in_one   = TRUE,
                               save_csv          = TRUE,
                               mr_horse_n_iter   = 5000,
                               mr_horse_n_burnin = 1000,
                               mr_grip_parameters = NULL,
                               combined_file = "heart atrial appendage.csv")

MRplots.multi(MR_input_data = pancreas_MR_input,
              methods.plot          = c('IVW','RAPS','Egger','PRESSO','Horse','GRIP'),
              NbDistribution_presso = 1000,
              SignifThreshold_presso= 0.05,
              mr_horse_n_iter       = 5000,
              mr_horse_n_burnin     = 1000,
              plot.xlab = "Per allele assoc. with",
              plot.ylab = "Per allele assoc. with",
              show.legend = TRUE)

# MRplots.multi(MR_input_data = pancreas_MR_input,
#               methods.plot          = c('IVW','RAPS','Egger','PRESSO','Horse','GRIP'),
#               NbDistribution_presso = 1000,
#               SignifThreshold_presso= 0.05,
#               mr_horse_n_iter       = 5000,
#               mr_horse_n_burnin     = 1000,
#               plot.xlab = "Per allele assoc. with",
#               plot.ylab = "Per allele assoc. with",
#               show.legend = FALSE)

MR_forest(pancreas,effect = "Beta", 
          # custom_xlim = c(-0.01,0.01), 
          save_plot = TRUE, file_type = "jpeg",
          dot_size = 8,
          axis_text_size = 20,
          axis_title_size = 20,
          pval_text_size = 8,
          plot_width = 20,
          plot_height = 8)

MR_forest_SNP(pancreas_MR_input, report_form = "Beta", 
              # xlim_custom = c(-0.02,0.02),
              save_plot = TRUE, file_type = "jpeg",
              dot_size = 8,
              axis_text_size = 20,
              axis_title_size = 20,
              digits = 3,
              plot_width = 20,
              plot_height = 10,
              label_text_size = 8)
