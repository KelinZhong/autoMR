# ==============================================================================
# Global variables declaration for R CMD check
# ==============================================================================
utils::globalVariables(c(
  "Outcome", "Exposure", "Instrument", "ALLELE0", "ALLELE1", "A1FREQ",
  "beta_exposure", "beta_outcome", "se_exposure", "se_outcome"
))

# ==============================================================================
# Internal Helper Functions (Calculation & Lookup Engines)
# ==============================================================================

#' Extract method slope from summary table
#' @noRd
get_slope_from_summary <- function(df_row, m, effect_scale) {
  if (is.null(df_row) || nrow(df_row) == 0) return(NA_real_)

  col <- if (m == "PRESSO") paste0("Presso_", effect_scale)
  else if (m == "Horse") paste0("Horse_", effect_scale)
  else if (m == "GRIP") paste0("Grip_", effect_scale)
  else paste0(m, "_", effect_scale)

  if (!col %in% names(df_row)) return(NA_real_)
  val <- suppressWarnings(as.numeric(df_row[[col]]))

  if (effect_scale != "Beta") {
    if (is.na(val) || !is.finite(val) || val <= 0) return(NA_real_)
    return(log(val))
  }
  val
}

#' Extract p-value from summary table
#' @noRd
get_p_from_summary <- function(df_row, m) {
  if (is.null(df_row) || nrow(df_row) == 0) return(NA_real_)
  col <- if (m == "Egger") "Egger_P_value"
  else if (m == "PRESSO") "Presso_p"
  else if (m == "Horse") "Horse_P"
  else if (m == "GRIP") "Grip_P"
  else paste0(m, "_P")
  if (!col %in% names(df_row)) return(NA_real_)
  suppressWarnings(as.numeric(df_row[[col]]))
}

#' Internal single-outcome analysis engine
#' @noRd
valid.output <- function(MR_input_data,
                         outcome.form = "Beta",
                         use_ivw = TRUE,
                         use_raps = TRUE,
                         use_median = TRUE,
                         use_egger = TRUE,
                         use_mr_presso = TRUE,
                         use_mr_horse = TRUE,
                         use_mr_grip  = TRUE,
                         NbDistribution = 1000,
                         SignifThreshold = 0.05,
                         custom_label = NULL,
                         return.result = TRUE,
                         save_csv = FALSE,
                         save_dir = ".",
                         mr_horse_n_iter = 5000,
                         mr_horse_n_burnin = 1000,
                         mr_grip_parameters = NULL) {

  MR_input_data <- ensure_dummy_vars(MR_input_data)
  EXP <- unique(MR_input_data$Exposure)
  results_list <- list()

  for(i in seq_along(EXP)){
    clean_Exposure.i <- as.data.frame(MR_input_data[MR_input_data$Exposure == EXP[i], ])
    exposure_beta <- clean_Exposure.i$beta_exposure
    outcome_beta  <- clean_Exposure.i$beta_outcome
    exposure_se   <- clean_Exposure.i$se_exposure
    outcome_se    <- clean_Exposure.i$se_outcome

    MRInputObject <- MendelianRandomization::mr_input(
      bx = exposure_beta, bxse = exposure_se,
      by = outcome_beta,  byse = outcome_se,
      snps = clean_Exposure.i$Instrument,
      effect_allele = clean_Exposure.i$ALLELE0,
      other_allele = clean_Exposure.i$ALLELE1,
      eaf = clean_Exposure.i$A1FREQ
    )

    result_row <- list()
    result_row[["Outcome"]]  <- unique(MR_input_data$Outcome)[1]
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
        for(nm in c("SNPs","IVW_Beta","IVW_Lower","IVW_Upper","IVW_P","IVW_Q","IVW_Q_P","Fstat")) result_row[[nm]] <- NA
      })
    } else {
      for(nm in c("SNPs","IVW_Beta","IVW_Lower","IVW_Upper","IVW_P","IVW_Q","IVW_Q_P","Fstat")) result_row[[nm]] <- NA
    }

    # --- MR-RAPS ---
    if(use_raps){
      tryCatch({
        r <- mr.raps(exposure_beta, outcome_beta, exposure_se, outcome_se)
        result_row[["RAPS_Beta"]]  <- r$beta.hat
        result_row[["RAPS_Lower"]] <- r$beta.hat - stats::qnorm(0.975)*r$beta.se
        result_row[["RAPS_Upper"]] <- r$beta.hat + stats::qnorm(0.975)*r$beta.se
        result_row[["RAPS_P"]]     <- r$beta.p.value
      }, error = function(e){
        for(nm in c("RAPS_Beta","RAPS_Lower","RAPS_Upper","RAPS_P")) result_row[[nm]] <- NA
      })
    } else {
      for(nm in c("RAPS_Beta","RAPS_Lower","RAPS_Upper","RAPS_P")) result_row[[nm]] <- NA
    }

    if(nrow(clean_Exposure.i) > 2){

      # --- Median ---
      if(use_median){
        tryCatch({
          m <- MendelianRandomization::mr_median(MRInputObject)
          result_row[["Med_Beta"]]  <- m@Estimate
          result_row[["Med_Lower"]] <- m@CILower
          result_row[["Med_Upper"]] <- m@CIUpper
          result_row[["Med_P"]]     <- m@Pvalue
        }, error = function(e){
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
          for(nm in c("Egger_Beta","Egger_Lower","Egger_Upper","Egger_P_value",
                      "Egger_Q","Egger_Q_P","I_sq","Intercept_Est",
                      "Intercept_Lower","Intercept_Upper","Intercept_P")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Egger_Beta","Egger_Lower","Egger_Upper","Egger_P_value",
                    "Egger_Q","Egger_Q_P","I_sq","Intercept_Est",
                    "Intercept_Lower","Intercept_Upper","Intercept_P")) result_row[[nm]] <- NA
      }

      # --- PRESSO ---
      if(use_mr_presso && nrow(clean_Exposure.i) >= 4){
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

          orig <- nrow(clean_Exposure.i)
          if(!is.null(pr$`MR-PRESSO results`$`Outlier Test`) && nrow(pr$`MR-PRESSO results`$`Outlier Test`) > 0){
            nout <- nrow(pr$`MR-PRESSO results`$`Outlier Test`)
            result_row[["Presso_SNPs"]]  <- orig - nout
            result_row[["outlier_SNPs"]] <- paste(row.names(pr$`MR-PRESSO results`$`Outlier Test`), collapse = ",")
          } else {
            result_row[["Presso_SNPs"]]  <- orig
            result_row[["outlier_SNPs"]] <- NA
          }
        }, error = function(e){
          for(nm in c("Presso_Beta","Presso_lower","Presso_upper","Presso_p","Presso_SNPs","outlier_SNPs")) result_row[[nm]] <- NA
        })
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
          b <- as.numeric(ho$MR_Estimate$Estimate[1])
          se <- as.numeric(ho$MR_Estimate$SD[1])
          result_row[["Horse_Beta"]]  <- b
          result_row[["Horse_Lower"]] <- b - stats::qnorm(0.975)*se
          result_row[["Horse_Upper"]] <- b + stats::qnorm(0.975)*se
          result_row[["Horse_P"]]     <- 2*(1 - stats::pnorm(abs(b/se)))
        }, error = function(e){
          for(nm in c("Horse_Beta","Horse_Lower","Horse_Upper","Horse_P")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Horse_Beta","Horse_Lower","Horse_Upper","Horse_P")) result_row[[nm]] <- NA
      }

      # --- GRIP ---
      if(use_mr_grip){
        tryCatch({
          params <- if(is.null(mr_grip_parameters)) default_parameters_mr_grip() else mr_grip_parameters
          gr <- mr_grip(
            b_exp = exposure_beta, b_out = outcome_beta,
            se_exp = exposure_se,  se_out = outcome_se,
            parameters = params
          )
          result_row[["Grip_Beta"]]  <- as.numeric(gr$b)
          result_row[["Grip_Lower"]] <- as.numeric(gr$b - 1.96*gr$se)
          result_row[["Grip_Upper"]] <- as.numeric(gr$b + 1.96*gr$se)
          result_row[["Grip_P"]]     <- as.numeric(gr$pval)
        }, error = function(e){
          for(nm in c("Grip_Beta","Grip_Lower","Grip_Upper","Grip_P")) result_row[[nm]] <- NA
        })
      } else {
        for(nm in c("Grip_Beta","Grip_Lower","Grip_Upper","Grip_P")) result_row[[nm]] <- NA
      }

    } else {
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

  if(outcome.form %in% c("OR","HR")){
    exp_cols <- c("IVW_Beta","RAPS_Beta","Med_Beta","Egger_Beta","Presso_Beta","Horse_Beta","Grip_Beta",
                  "IVW_Lower","RAPS_Lower","Med_Lower","Egger_Lower","Presso_lower","Horse_Lower","Grip_Lower",
                  "IVW_Upper","RAPS_Upper","Med_Upper","Egger_Upper","Presso_upper","Horse_Upper","Grip_Upper")
    for(col in exp_cols) if(col %in% colnames(mr_res)) mr_res[[col]] <- exp(mr_res[[col]])

    prefix <- if (outcome.form == "OR") "_OR" else "_HR"
    for(pat in c("IVW","RAPS","Med","Egger","Presso","Horse","Grip")){
      names(mr_res) <- sub(paste0("^",pat,"_Beta$"), paste0(pat, prefix), names(mr_res))
    }
  }

  if(all(c("Fstat","IVW_Q_P") %in% colnames(mr_res))){
    mr_res$FLe10 <- if(use_ivw) as.integer(mr_res$Fstat < 10) else NA
    mr_res$SigQ  <- if(use_ivw) as.integer(mr_res$IVW_Q_P < 0.01) else NA
  }
  if("Intercept_P" %in% colnames(mr_res)){
    mr_res$I2Le90     <- if(use_egger) as.integer(mr_res$I_sq < 0.9) else NA
    mr_res$Pleiotropy <- if(use_egger) as.integer(mr_res$Intercept_P < 0.01) else NA
  }

  new_order <- c("Outcome","Exposure","SNPs","Presso_SNPs","outlier_SNPs","FLe10","SigQ","I2Le90","Pleiotropy",
                 setdiff(colnames(mr_res), c("Outcome","Exposure","SNPs","Presso_SNPs","outlier_SNPs","FLe10","SigQ","I2Le90","Pleiotropy")))
  mr_res <- mr_res[, new_order]

  if(save_csv){
    if(is.null(custom_label)) custom_label <- paste0(unique(MR_input_data$Outcome)[1], "_results.csv")
    utils::write.csv(mr_res, file = file.path(save_dir, custom_label), row.names = FALSE)
  }

  mr_res
}

#' Internal single-outcome scatter plot engine
#' @noRd
#' @importFrom grDevices png dev.off dev.cur
MRplots <- function(MR_input_data,
                    d.title = NULL,
                    d.subtitle = NULL,
                    plot.xlab,
                    plot.ylab,
                    outcome_label,
                    methods.plot = c("IVW","RAPS","Egger","PRESSO","Horse","GRIP"),
                    NbDistribution_presso = 1000,
                    SignifThreshold_presso = 0.05,
                    mr_horse_n_iter = 5000,
                    mr_horse_n_burnin = 1000,
                    show.legend = TRUE,
                    summary_df = NULL,
                    effect_scale = "Beta",
                    use_df_results = TRUE,
                    save_plot = FALSE,
                    save_dir = ".") {

  MR_input_data <- ensure_dummy_vars(MR_input_data)

  if (is.null(summary_df) || !use_df_results) {
    summary_df <- valid.output(
      MR_input_data, outcome.form = effect_scale, return.result = TRUE, save_csv = FALSE,
      NbDistribution = NbDistribution_presso, SignifThreshold = SignifThreshold_presso,
      mr_horse_n_iter = mr_horse_n_iter, mr_horse_n_burnin = mr_horse_n_burnin
    )
  }

  colors <- list(IVW = "chartreuse3", RAPS = "turquoise3", Egger = "cornflowerblue",
                 PRESSO = "red", Horse = "purple", GRIP = "darkorange")

  for (ex in unique(MR_input_data$Exposure)) {
    dat <- subset(MR_input_data, Exposure == ex)
    df_row <- subset(summary_df, Outcome == outcome_label & Exposure == ex)
    if (nrow(df_row) == 0) next

    d.x <- dat$beta_exposure; d.y <- dat$beta_outcome
    d.x.se <- dat$se_exposure; d.y.se <- dat$se_outcome

    x_ext <- max(abs(d.x - d.x.se), abs(d.x + d.x.se), na.rm = TRUE)
    y_ext <- max(abs(d.y - d.y.se), abs(d.y + d.y.se), na.rm = TRUE)
    limx <- c(0, x_ext)
    limy <- c(-y_ext, y_ext)

    if (save_plot) {
      fname <- file.path(save_dir, paste0("Scatter_", gsub("\\s+", "_", ex), "_", gsub("\\s+", "_", outcome_label), ".png"))
      grDevices::png(fname, width = 8, height = 6, units = "in", res = 300)
    }

    graphics::plot(d.x, d.y, xlim = limx, ylim = limy,
                   xlab = paste(plot.xlab, ex), ylab = paste(plot.ylab, outcome_label),
                   main = d.title, pch = 16, bty = "L")
    graphics::mtext(d.subtitle)
    graphics::abline(h = 0, lty = 2); graphics::abline(v = 0, lty = 2)
    graphics::arrows(d.x, d.y - d.y.se, d.x, d.y + d.y.se, length = 0, angle = 90, code = 3, col = "grey")
    graphics::arrows(d.x - d.x.se, d.y, d.x + d.x.se, d.y, length = 0, angle = 90, code = 3, col = "grey")

    legend_txt <- character(); legend_cols <- character(); legend_lty <- integer()

    for (m_info in list(
      list(name="IVW", lty=3), list(name="RAPS", lty=3), list(name="Egger", lty=4),
      list(name="PRESSO", lty=5), list(name="Horse", lty=6), list(name="GRIP", lty=2)
    )) {
      if (m_info$name %in% methods.plot) {
        b <- get_slope_from_summary(df_row, m_info$name, effect_scale)
        p <- get_p_from_summary(df_row, m_info$name)

        if (!is.na(b) && is.finite(b)) {
          graphics::abline(a = 0, b = b, lty = m_info$lty, lwd = 2, col = colors[[m_info$name]])
          legend_txt <- c(legend_txt, sprintf("%s: beta=%.3f, p=%.3g", m_info$name, b, p))
        } else {
          legend_txt <- c(legend_txt, sprintf("%s: failed", m_info$name))
        }
        legend_cols <- c(legend_cols, colors[[m_info$name]])
        legend_lty <- c(legend_lty, m_info$lty)
      }
    }

    if (show.legend && length(legend_txt) > 0) {
      graphics::legend("topright", legend = legend_txt, col = legend_cols, lty = legend_lty, lwd = 2, bg = "white", cex = 0.9)
    }

    if (save_plot && grDevices::dev.cur() != 1) grDevices::dev.off()
  }
}

# ==============================================================================
# Exported User-Facing Functions
# ==============================================================================

#' Run MR Analysis for Multiple Outcomes
#'
#' Performs causal inference analysis using multiple Mendelian Randomization methods
#' across one or more outcomes and exposures.
#'
#' @param MR_input_data Harmonised MR input data frame. Must contain Outcome and Exposure columns.
#' @param outcome.form Character vector indicating standard output scale for each outcome (e.g., "Beta", "OR", "HR"). Defaults to "Beta".
#' @param use_ivw Logical; whether to use Inverse Variance Weighted method.
#' @param use_raps Logical; whether to use Robust Adjusted Profile Score method.
#' @param use_median Logical; whether to use Weighted Median method.
#' @param use_egger Logical; whether to use MR-Egger regression.
#' @param use_mr_presso Logical; whether to use MR-PRESSO method.
#' @param use_mr_horse Logical; whether to use MR-Horse method.
#' @param use_mr_grip Logical; whether to use MR-GRIP method.
#' @param NbDistribution Integer; number of distributions for MR-PRESSO. Default is 1000.
#' @param SignifThreshold Numeric; significance threshold for MR-PRESSO. Default is 0.05.
#' @param save_all_in_one Logical; whether to save combined results of all outcomes into a single file.
#' @param save_csv Logical; whether to save results as CSV files. Default is TRUE.
#' @param combined_file Character string for the combined output filename. Default is "combined_results.csv".
#' @param save_dir Directory to save output files. Default is current directory (".").
#' @param mr_horse_n_iter Integer; number of iterations for MR-Horse. Default is 5000.
#' @param mr_horse_n_burnin Integer; number of burn-in samples for MR-Horse. Default is 1000.
#' @param mr_grip_parameters List; additional parameters for the MR-GRIP method.
#' @return Combined results data frame.
#' @examples
#' \donttest{
#' data("fi_49item")
#' data("fried_frailty")
#' input1 <- harmonize_mr_data(df = fi_49item)
#' input2 <- harmonize_mr_data(df = fried_frailty)
#' input1 <- harmonize_mr_data(df = fi_49item)
#' input2 <- harmonize_mr_data(df = fried_frailty)
#'
#' outcome1 <- run_mr_analysis(
#'   MR_input_data = input1,
#'   outcome.form = NULL, # defaul is "Beta"
#'   use_ivw = TRUE,
#'   use_raps = TRUE,
#'   use_median = TRUE,
#'   use_egger = TRUE,
#'   use_mr_presso = TRUE,
#'   use_mr_horse = TRUE,
#'   use_mr_grip  = TRUE,
#'   NbDistribution    = 1000,
#'   SignifThreshold   = 0.05,
#'   save_all_in_one   = TRUE,
#'   save_csv          = FALSE,
#'   save_dir          = ".",
#'   mr_horse_n_iter   = 5000,
#'   mr_horse_n_burnin = 1000,
#'   mr_grip_parameters = NULL,
#'   combined_file = "fi_49item results.csv" # you can customize the output csv file name
#'   )
#'
#' outcome2 <- run_mr_analysis(
#'   MR_input_data = input2,
#'   outcome.form = "OR", # defaul is "Beta"
#'   use_ivw = TRUE,
#'   use_raps = TRUE,
#'   use_median = TRUE,
#'   use_egger = TRUE,
#'   use_mr_presso = TRUE,
#'   use_mr_horse = TRUE,
#'   use_mr_grip  = TRUE,
#'   NbDistribution    = 1000,
#'   SignifThreshold   = 0.05,
#'   save_all_in_one   = TRUE,
#'   save_csv          = FALSE,
#'   save_dir          = ".",
#'   mr_horse_n_iter   = 5000,
#'   mr_horse_n_burnin = 1000,
#'   mr_grip_parameters = NULL,
#'   combined_file = "fried_frailty results.csv" # you can customize the output csv file name
#' )
#' }
#' @export
run_mr_analysis <- function(MR_input_data,
                            outcome.form = NULL,
                            use_ivw = TRUE,
                            use_raps = TRUE,
                            use_median = TRUE,
                            use_egger = TRUE,
                            use_mr_presso = TRUE,
                            use_mr_horse = TRUE,
                            use_mr_grip  = TRUE,
                            NbDistribution = 1000,
                            SignifThreshold = 0.05,
                            save_all_in_one = FALSE,
                            save_csv = TRUE,
                            combined_file = "combined_results.csv",
                            save_dir = ".",
                            mr_horse_n_iter = 5000,
                            mr_horse_n_burnin = 1000,
                            mr_grip_parameters = NULL) {

  outcomes <- unique(MR_input_data$Outcome)
  if(is.null(outcome.form)) outcome.form <- rep("Beta", length(outcomes))
  if(length(outcome.form) == 1) outcome.form <- rep(outcome.form, length(outcomes))

  results_list <- list()
  for(i in seq_along(outcomes)){
    Outcome.i <- as.data.frame(MR_input_data[MR_input_data$Outcome == outcomes[i], ])
    res <- valid.output(
      MR_input_data = Outcome.i,
      outcome.form = outcome.form[i],
      use_ivw = use_ivw,
      use_raps = use_raps,
      use_median = use_median,
      use_egger = use_egger,
      use_mr_presso = use_mr_presso,
      use_mr_horse = use_mr_horse,
      use_mr_grip = use_mr_grip,
      NbDistribution = NbDistribution,
      SignifThreshold = SignifThreshold,
      save_csv = save_csv,
      save_dir = save_dir,
      mr_horse_n_iter = mr_horse_n_iter,
      mr_horse_n_burnin = mr_horse_n_burnin,
      mr_grip_parameters = mr_grip_parameters
    )
    results_list[[i]] <- res
  }
  combined_results <- do.call(rbind, results_list)

  if(save_all_in_one && save_csv){
    utils::write.csv(combined_results, file = file.path(save_dir, combined_file), row.names = FALSE)
  }
  return(combined_results)
}

#' Plot MR Scatter plots for Multiple Outcomes
#'
#' Generates scatter plots with fitted lines for each MR method. Efficiently
#' uses pre-calculated results if provided via summary_df.
#'
#' @param MR_input_data Harmonised MR input data frame. Must contain Outcome and Exposure columns.
#' @param plot.xlab Character string; prefix for the X-axis label. Default is "Exposure".
#' @param plot.ylab Character string; prefix for the Y-axis label. Default is "Outcome".
#' @param methods.plot Character vector; methods to overlay on the plot. Default includes all supported methods.
#' @param NbDistribution_presso Integer; distributions for on-the-fly MR-PRESSO calculation. Default 1000.
#' @param SignifThreshold_presso Numeric; significance threshold for on-the-fly MR-PRESSO calculation. Default 0.05.
#' @param mr_horse_n_iter Integer; iterations for on-the-fly MR-Horse calculation. Default 5000.
#' @param mr_horse_n_burnin Integer; burn-in samples for on-the-fly MR-Horse calculation. Default 1000.
#' @param show.legend Logical; whether to display the method legend. Default is TRUE.
#' @param summary_df Data frame; optional pre-calculated summary results from run_mr_analysis().
#' @param effect_scale Character string; scale of effect to match summary_df ("Beta", "OR", "HR"). Default "Beta".
#' @param use_df_results Logical; if TRUE, uses results from summary_df instead of re-calculating. Default is TRUE.
#' @param save_plot Logical; whether to save plots to disk. Default is FALSE.
#' @param save_dir Directory to save plots. Default is current directory (".").
#' @return Invisibly returns NULL.
#' @examples
#' \donttest{
#' data("fi_49item")
#' input1 <- harmonize_mr_data(df = fi_49item)
#' outcome1 <- run_mr_analysis(
#'   MR_input_data = input1,
#'   outcome.form = NULL, # defaul is "Beta"
#'   use_ivw = TRUE,
#'   use_raps = TRUE,
#'   use_median = TRUE,
#'   use_egger = TRUE,
#'   use_mr_presso = TRUE,
#'   use_mr_horse = TRUE,
#'   use_mr_grip  = TRUE,
#'   NbDistribution    = 1000,
#'   SignifThreshold   = 0.05,
#'   save_all_in_one   = TRUE,
#'   save_csv          = FALSE,
#'   save_dir          = ".",
#'   mr_horse_n_iter   = 5000,
#'   mr_horse_n_burnin = 1000,
#'   mr_grip_parameters = NULL,
#'   combined_file = "fi_49item results.csv" # you can customize the output csv file name
#'   )
#'
#' plot_mr_scatter(
#'   MR_input_data = input1,
#'   plot.xlab = "Exposure",
#'   plot.ylab = "Outcome",
#'   methods.plot = c("IVW", "RAPS", "Egger", "PRESSO", "Horse", "GRIP"),
#'   NbDistribution_presso = 1000,
#'   SignifThreshold_presso = 0.05,
#'   mr_horse_n_iter = 5000,
#'   mr_horse_n_burnin = 1000,
#'   show.legend = TRUE,
#'   summary_df = outcome1, # avoid computing the MR results again if you already run run_mr_analysis()
#'   effect_scale = "Beta",
#'   use_df_results = TRUE,
#'   save_plot = FALSE,
#'   save_dir = "."
#' )
#'
#' plot_mr_scatter(
#'   MR_input_data = input1,
#'   plot.xlab = "Exposure",
#'   plot.ylab = "Outcome",
#'   methods.plot = c("IVW", "RAPS", "Egger", "PRESSO", "Horse", "GRIP"),
#'   NbDistribution_presso = 1000,
#'   SignifThreshold_presso = 0.05,
#'   mr_horse_n_iter = 5000,
#'   mr_horse_n_burnin = 1000,
#'   show.legend = TRUE,
#'   summary_df = NULL,     # you can also generate the plot without the results from run_mr_analysis()
#'   effect_scale = "Beta",
#'   use_df_results = FALSE,
#'   save_plot = FALSE,
#'   save_dir = "."
#' )
#' }
#' @export
plot_mr_scatter <- function(MR_input_data,
                            plot.xlab = "Exposure",
                            plot.ylab = "Outcome",
                            methods.plot = c("IVW","RAPS","Egger","PRESSO","Horse","GRIP"),
                            NbDistribution_presso = 1000,
                            SignifThreshold_presso = 0.05,
                            mr_horse_n_iter = 5000,
                            mr_horse_n_burnin = 1000,
                            show.legend = TRUE,
                            summary_df = NULL,
                            effect_scale = "Beta",
                            use_df_results = TRUE,
                            save_plot = FALSE,
                            save_dir = ".") {
  for (out in unique(MR_input_data$Outcome)) {
    sub_data <- subset(MR_input_data, Outcome == out)
    MRplots(
      MR_input_data = sub_data,
      outcome_label = out,
      plot.xlab = plot.xlab,
      plot.ylab = plot.ylab,
      methods.plot = methods.plot,
      NbDistribution_presso = NbDistribution_presso,
      SignifThreshold_presso = SignifThreshold_presso,
      mr_horse_n_iter = mr_horse_n_iter,
      mr_horse_n_burnin = mr_horse_n_burnin,
      show.legend = show.legend,
      summary_df = summary_df,
      effect_scale = effect_scale,
      use_df_results = use_df_results,
      save_plot = save_plot,
      save_dir = save_dir
    )
  }
  invisible(NULL)
}
