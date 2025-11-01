# valid_output.R
# ---------------------------------------------------
# Summarise MR methods (IVW, RAPS, Median, Egger, PRESSO, Horse, GRIP)
# Outputs columns compatible with MR_forest() and MRplots() below.
# ---------------------------------------------------

valid.output <- function(MR_input_data, 
                         outcome.form = "Beta",
                         use_ivw = TRUE, 
                         use_raps = TRUE, 
                         use_median = TRUE, 
                         use_egger = TRUE, 
                         use_mr_presso = TRUE,
                         use_mr_horse = TRUE,
                         use_mr_grip  = TRUE,          # <-- GRIP on/off
                         NbDistribution = 1000,
                         SignifThreshold = 0.05,
                         custom_label,
                         return.result = FALSE,
                         save_csv = TRUE,
                         mr_horse_n_iter = 5000,
                         mr_horse_n_burnin = 1000,
                         mr_grip_parameters = NULL) {  # <-- optional GRIP params
  require(MendelianRandomization)
  require(TwoSampleMR)
  require(mr.raps)
  require(MRPRESSO)
  
  EXP <- unique(MR_input_data$Exposure)
  results_list <- list()
  
  for(i in seq_along(EXP)){
    Exposure.i <- as.data.frame(MR_input_data[MR_input_data$Exposure == EXP[i], ])
    clean_Exposure.i <- Exposure.i
    
    exposure_beta <- clean_Exposure.i$beta_exposure
    outcome_beta  <- clean_Exposure.i$beta_outcome
    exposure_se   <- clean_Exposure.i$se_exposure
    outcome_se    <- clean_Exposure.i$se_outcome
    
    MRInputObject <- mr_input(
      bx = exposure_beta, bxse = exposure_se,
      by = outcome_beta,  byse = outcome_se
    )
    
    result_row <- list()
    outcome_value <- unique(MR_input_data$Outcome)[1]
    result_row[["Outcome"]]  <- outcome_value
    result_row[["Exposure"]] <- EXP[i]
    
    # --- IVW ---
    if(use_ivw){
      tryCatch({
        ivw_obj <- MendelianRandomization::mr_ivw(MRInputObject)
        result_row[["SNPs"]]      <- ivw_obj@SNPs
        result_row[["IVW_Beta"]]  <- ivw_obj@Estimate
        result_row[["IVW_Lower"]] <- ivw_obj@CILower
        result_row[["IVW_Upper"]] <- ivw_obj@CIUpper
        result_row[["IVW_P"]]     <- ivw_obj@Pvalue
        heter <- ivw_obj@Heter.Stat; if(length(heter)==1) heter <- c(heter, NA)
        result_row[["IVW_Q"]]     <- heter[1]
        result_row[["IVW_Q_P"]]   <- heter[2]
        result_row[["Fstat"]]     <- ivw_obj@Fstat
      }, error = function(e){
        warning(sprintf("IVW failed for %s: %s", EXP[i], e$message))
        for(nm in c("SNPs","IVW_Beta","IVW_Lower","IVW_Upper","IVW_P","IVW_Q","IVW_Q_P","Fstat")) result_row[[nm]] <- NA
      })
    } else {
      for(nm in c("SNPs","IVW_Beta","IVW_Lower","IVW_Upper","IVW_P","IVW_Q","IVW_Q_P","Fstat")) result_row[[nm]] <- NA
    }
    
    # --- MR-RAPS ---
    if(use_raps){
      if(length(exposure_beta)>0 && length(outcome_beta)>0 &&
         length(exposure_se)>0   && length(outcome_se)>0 &&
         all(is.finite(exposure_beta), is.finite(outcome_beta),
             is.finite(exposure_se),  is.finite(outcome_se)) &&
         all(exposure_se>0, outcome_se>0)){
        tryCatch({
          r <- mr.raps(exposure_beta, outcome_beta, exposure_se, outcome_se)
          result_row[["RAPS_Beta"]]  <- r$beta.hat
          result_row[["RAPS_Lower"]] <- r$beta.hat - qnorm(0.975)*r$beta.se
          result_row[["RAPS_Upper"]] <- r$beta.hat + qnorm(0.975)*r$beta.se
          result_row[["RAPS_P"]]     <- r$beta.p.value
        }, error = function(e){
          warning(sprintf("MR-RAPS failed for %s: %s", EXP[i], e$message))
          for(nm in c("RAPS_Beta","RAPS_Lower","RAPS_Upper","RAPS_P")) result_row[[nm]] <- NA
        })
      } else {
        warning(sprintf("Invalid data for MR-RAPS: %s", EXP[i]))
        for(nm in c("RAPS_Beta","RAPS_Lower","RAPS_Upper","RAPS_P")) result_row[[nm]] <- NA
      }
    } else {
      for(nm in c("RAPS_Beta","RAPS_Lower","RAPS_Upper","RAPS_P")) result_row[[nm]] <- NA
    }
    
    if(nrow(Exposure.i) > 2){
      # --- Median ---
      if(use_median){
        tryCatch({
          m <- MendelianRandomization::mr_median(MRInputObject)
          result_row[["Med_Beta"]]  <- m@Estimate
          result_row[["Med_Lower"]] <- m@CILower
          result_row[["Med_Upper"]] <- m@CIUpper
          result_row[["Med_P"]]     <- m@Pvalue
        }, error = function(e){
          warning(sprintf("Median failed for %s: %s", EXP[i], e$message))
          for(nm in c("Med_Beta","Med_Lower","Med_Upper","Med_P")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Med_Beta","Med_Lower","Med_Upper","Med_P")) result_row[[nm]] <- NA
      }
      
      # --- Egger ---
      if(use_egger){
        tryCatch({
          eg <- MendelianRandomization::mr_egger(MRInputObject)
          heter <- eg@Heter.Stat; if(length(heter)==1) heter <- c(heter, NA)
          result_row[["Egger_Beta"]]       <- eg@Estimate
          result_row[["Egger_Lower"]]      <- eg@CILower.Est
          result_row[["Egger_Upper"]]      <- eg@CIUpper.Est
          result_row[["Egger_P_value"]]    <- eg@Pvalue.Est
          result_row[["Egger_Q"]]          <- heter[1]
          result_row[["Egger_Q_P"]]        <- heter[2]
          result_row[["I_sq"]]             <- eg@I.sq
          result_row[["Intercept_Est"]]    <- eg@Intercept
          result_row[["Intercept_Lower"]]  <- eg@CILower.Int
          result_row[["Intercept_Upper"]]  <- eg@CIUpper.Int
          result_row[["Intercept_P"]]      <- eg@Pvalue.Int
        }, error = function(e){
          warning(sprintf("Egger failed for %s: %s", EXP[i], e$message))
          for(nm in c("Egger_Beta","Egger_Lower","Egger_Upper","Egger_P_value",
                      "Egger_Q","Egger_Q_P","I_sq","Intercept_Est",
                      "Intercept_Lower","Intercept_Upper","Intercept_P")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Egger_Beta","Egger_Lower","Egger_Upper","Egger_P_value",
                    "Egger_Q","Egger_Q_P","I_sq","Intercept_Est",
                    "Intercept_Lower","Intercept_Upper","Intercept_P")) result_row[[nm]] <- NA
      }
      
      # --- PRESSO (>=4 SNPs) ---
      if(use_mr_presso){
        if(nrow(Exposure.i) >= 4){
          tryCatch({
            pr <- mr_presso(
              BetaOutcome = "beta_outcome", BetaExposure = "beta_exposure",
              SdOutcome   = "se_outcome",   SdExposure   = "se_exposure",
              data = clean_Exposure.i,
              OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
              NbDistribution = NbDistribution, SignifThreshold = SignifThreshold
            )
            b  <- pr$`Main MR results`[1,3]
            se <- pr$`Main MR results`[1,4]
            p  <- pr$`Main MR results`[1,6]
            result_row[["Presso_Beta"]]  <- b
            result_row[["Presso_lower"]] <- b - 1.96*se
            result_row[["Presso_upper"]] <- b + 1.96*se
            result_row[["Presso_p"]]     <- p
            # Diagnostics
            orig <- nrow(clean_Exposure.i)
            if(!is.null(pr$`Outlier Test`) && nrow(pr$`Outlier Test`) > 0){
              nout <- nrow(pr$`Outlier Test`)
              result_row[["Presso_SNPs"]]  <- orig - nout
              result_row[["outlier_SNPs"]] <- paste(pr$`Outlier Test`$SNP, collapse = ",")
            } else {
              result_row[["Presso_SNPs"]]  <- orig
              result_row[["outlier_SNPs"]] <- NA
            }
          }, error = function(e){
            warning(sprintf("PRESSO failed for %s: %s", EXP[i], e$message))
            for(nm in c("Presso_Beta","Presso_lower","Presso_upper","Presso_p","Presso_SNPs","outlier_SNPs")) result_row[[nm]] <- NA
          })
        } else {
          for(nm in c("Presso_Beta","Presso_lower","Presso_upper","Presso_p","Presso_SNPs","outlier_SNPs")) result_row[[nm]] <- NA
        }
      } else {
        for(nm in c("Presso_Beta","Presso_lower","Presso_upper","Presso_p","Presso_SNPs","outlier_SNPs")) result_row[[nm]] <- NA
      }
      
      # --- Horse ---
      if(use_mr_horse){
        tryCatch({
          D <- data.frame(betaY = outcome_beta, betaX = exposure_beta,
                          betaYse = outcome_se, betaXse = exposure_se)
          ho <- mr_horse(D, no_ini = 3, variable.names = "theta",
                         n.iter = mr_horse_n_iter, n.burnin = mr_horse_n_burnin)
          he <- ho$MR_Estimate
          b  <- he$Estimate[1]; se <- he$SD[1]
          l  <- b - qnorm(0.975)*se
          u  <- b + qnorm(0.975)*se
          p  <- 2*(1 - pnorm(abs(b/se)))
          result_row[["Horse_Beta"]]  <- b
          result_row[["Horse_Lower"]] <- l
          result_row[["Horse_Upper"]] <- u
          result_row[["Horse_P"]]     <- p
        }, error = function(e){
          warning(sprintf("MR-Horse failed for %s: %s", EXP[i], e$message))
          for(nm in c("Horse_Beta","Horse_Lower","Horse_Upper","Horse_P")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Horse_Beta","Horse_Lower","Horse_Upper","Horse_P")) result_row[[nm]] <- NA
      }
      
      # --- GRIP (TwoSampleMR) ---
      if(use_mr_grip){
        tryCatch({
          params <- if(is.null(mr_grip_parameters)) TwoSampleMR::default_parameters() else mr_grip_parameters
          gr <- TwoSampleMR::mr_grip(
            b_exp = exposure_beta, b_out = outcome_beta,
            se_exp = exposure_se,  se_out = outcome_se,
            parameters = params
          )
          # primary estimate
          result_row[["Grip_Beta"]]  <- as.numeric(gr$b)
          result_row[["Grip_Lower"]] <- as.numeric(gr$b - 1.96*gr$se)
          result_row[["Grip_Upper"]] <- as.numeric(gr$b + 1.96*gr$se)
          result_row[["Grip_P"]]     <- as.numeric(gr$pval)
          # (Optionally available adjusted estimates could be stored as *_adj if desired)
        }, error = function(e){
          warning(sprintf("MR-GRIP failed for %s: %s", EXP[i], e$message))
          for(nm in c("Grip_Beta","Grip_Lower","Grip_Upper","Grip_P")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Grip_Beta","Grip_Lower","Grip_Upper","Grip_P")) result_row[[nm]] <- NA
      }
      
    } else {
      # <3 SNPs: set NA for methods that need >2
      for(nm in c("Med_Beta","Med_Lower","Med_Upper","Med_P",
                  "Egger_Beta","Egger_Lower","Egger_Upper","Egger_P_value",
                  "Egger_Q","Egger_Q_P","I_sq","Intercept_Est",
                  "Intercept_Lower","Intercept_Upper","Intercept_P",
                  "Presso_Beta","Presso_lower","Presso_upper","Presso_p","Presso_SNPs","outlier_SNPs",
                  "Horse_Beta","Horse_Lower","Horse_Upper","Horse_P",
                  "Grip_Beta","Grip_Lower","Grip_Upper","Grip_P")) result_row[[nm]] <- NA
    }
    
    results_list[[i]] <- as.data.frame(result_row, stringsAsFactors = FALSE)
  }
  
  mr_res <- do.call(rbind, results_list)
  
  # Scale handling for OR/HR
  if(outcome.form %in% c("OR","HR")){
    exp_cols <- c("IVW_Beta","RAPS_Beta","Med_Beta","Egger_Beta","Presso_Beta","Horse_Beta","Grip_Beta",
                  "IVW_Lower","RAPS_Lower","Med_Lower","Egger_Lower","Presso_lower","Horse_Lower","Grip_Lower",
                  "IVW_Upper","RAPS_Upper","Med_Upper","Egger_Upper","Presso_upper","Horse_Upper","Grip_Upper")
    for(col in exp_cols) if(col %in% colnames(mr_res)) mr_res[[col]] <- exp(mr_res[[col]])
    # Rename Beta column to effect name for convenience in plotting; CIs stay as *_Lower/Upper (on effect scale)
    if(outcome.form == "OR"){
      for(pat in c("IVW","RAPS","Med","Egger","Presso","Horse","Grip")){
        names(mr_res) <- sub(paste0("^",pat,"_Beta$"), paste0(pat,"_OR"), names(mr_res))
      }
    }
    if(outcome.form == "HR"){
      for(pat in c("IVW","RAPS","Med","Egger","Presso","Horse","Grip")){
        names(mr_res) <- sub(paste0("^",pat,"_Beta$"), paste0(pat,"_HR"), names(mr_res))
      }
    }
  }
  
  # Diagnostics flags
  if(all(c("Fstat","IVW_Q_P") %in% colnames(mr_res))){
    if(use_ivw){
      mr_res$FLe10 <- as.integer(mr_res$Fstat < 10)
      mr_res$SigQ  <- as.integer(mr_res$IVW_Q_P < 0.01)
    } else {
      mr_res$FLe10 <- NA; mr_res$SigQ <- NA
    }
  }
  if("Intercept_P" %in% colnames(mr_res)){
    if(use_egger){
      mr_res$I2Le90   <- as.integer(mr_res$I_sq < 0.9)
      mr_res$Pleiotropy <- as.integer(mr_res$Intercept_P < 0.01)
    } else {
      mr_res$I2Le90 <- NA; mr_res$Pleiotropy <- NA
    }
  }
  
  new_order <- c("Outcome","Exposure","SNPs","Presso_SNPs","outlier_SNPs","FLe10","SigQ","I2Le90","Pleiotropy",
                 setdiff(colnames(mr_res), c("Outcome","Exposure","SNPs","Presso_SNPs","outlier_SNPs","FLe10","SigQ","I2Le90","Pleiotropy")))
  mr_res <- mr_res[, new_order]
  
  if(save_csv){
    if(missing(custom_label) || is.null(custom_label) || !nzchar(custom_label)){
      custom_label <- paste0(gsub("[^A-Za-z0-9._-]","_", outcome_value), "_",
                             gsub("[^A-Za-z0-9._-]","_", EXP[1]), "_results.csv")
    }
    write.csv(mr_res, file = custom_label, row.names = FALSE)
  }
  
  if(return.result) mr_res else invisible(NULL)
}


valid.multi.output <- function(MR_input_data,
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
                               save_all_in_one   = FALSE,
                               save_csv          = TRUE,
                               combined_file     = "combined_results.csv",
                               mr_horse_n_iter   = 5000,
                               mr_horse_n_burnin = 1000,
                               mr_grip_parameters = NULL) {
  outcome <- unique(MR_input_data$Outcome)
  if(is.null(outcome.form)) outcome.form <- rep("Beta", length(outcome))
  
  results_list <- list()
  for(i in seq_along(outcome)){
    Outcome.i <- as.data.frame(MR_input_data[MR_input_data$Outcome == outcome[i], ])
    res <- valid.output(
      MR_input_data = Outcome.i,
      outcome.form  = outcome.form[i],
      use_ivw = use_ivw, use_raps = use_raps, use_median = use_median,
      use_egger = use_egger, use_mr_presso = use_mr_presso,
      use_mr_horse = use_mr_horse, use_mr_grip = use_mr_grip,
      NbDistribution = NbDistribution, SignifThreshold = SignifThreshold,
      custom_label = paste(outcome[i], "filter_results.csv"),
      return.result = TRUE,
      save_csv = FALSE,  # keep per-outcome off here; control with save_all_in_one below
      mr_horse_n_iter = mr_horse_n_iter, mr_horse_n_burnin = mr_horse_n_burnin,
      mr_grip_parameters = mr_grip_parameters
    )
    results_list[[i]] <- res
  }
  combined_results <- do.call(rbind, results_list)
  
  if(save_all_in_one && save_csv){
    write.csv(combined_results, file = combined_file, row.names = FALSE)
  }
  invisible(combined_results)
}
