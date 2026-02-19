library(ggplot2)
library(dplyr)
library(MendelianRandomization)

MR_forest_SNP <- function(df, report_form, 
                          save_plot = FALSE, file_type = "jpeg", 
                          save_dir = ".", xlim_custom = NULL,
                          dot_size = 2,
                          axis_text_size = 10,
                          axis_title_size = 12,
                          digits = 2,
                          plot_width = 8,
                          plot_height = 6,
                          label_text_size = 3) {
  # df: must contain columns: Outcome, SNP, ALLELE1, ALLELE0, A1FREQ,
  #     beta_exposure, se_exposure, beta_outcome, se_outcome, Exposure.
  # report_form: vector matching unique outcomes.
  
  unique_outcomes <- unique(df$Outcome)
  if (length(report_form) == 1) {
    report_form <- rep(report_form, length(unique_outcomes))
  }
  if (length(report_form) != length(unique_outcomes)) {
    stop("Length of report_form must be 1 or equal to the number of unique outcomes in the data frame.")
  }
  
  outcome_exposure <- df %>% distinct(Outcome, Exposure)
  outcome_exposure$effect <- sapply(outcome_exposure$Outcome, function(x) {
    report_form[which(unique_outcomes == x)]
  })
  
  plots_list <- list()
  for (i in seq_len(nrow(outcome_exposure))) {
    current_outcome  <- outcome_exposure$Outcome[i]
    current_exposure <- outcome_exposure$Exposure[i]
    current_form     <- as.character(outcome_exposure$effect[i])
    
    df_sub <- df %>% 
      filter(Outcome == current_outcome, Exposure == current_exposure) %>% 
      mutate(
        Estimate = beta_outcome / beta_exposure,
        se_ratio = sqrt((se_outcome / beta_exposure)^2 +
                          ((beta_outcome * se_exposure) / (beta_exposure^2))^2),
        Lower    = Estimate - 1.96 * se_ratio,
        Upper    = Estimate + 1.96 * se_ratio
      )
    
    if (current_form != "Beta") {
      df_sub <- df_sub %>% mutate(
        Estimate = exp(Estimate),
        Lower    = exp(Lower),
        Upper    = exp(Upper)
      )
      ref_line <- 1
    } else {
      ref_line <- 0
    }
    
    MRInput    <- mr_input(bx = df_sub$beta_exposure,
                           bxse = df_sub$se_exposure,
                           by   = df_sub$beta_outcome,
                           byse = df_sub$se_outcome)
    ivw_obj    <- MendelianRandomization::mr_ivw(MRInput)
    pooled_beta  <- ivw_obj@Estimate
    pooled_lower <- ivw_obj@CILower
    pooled_upper <- ivw_obj@CIUpper
    if (current_form != "Beta") {
      pooled_beta  <- exp(pooled_beta)
      pooled_lower <- exp(pooled_lower)
      pooled_upper <- exp(pooled_upper)
    }
    pooled_row <- tibble(
      SNP      = "IVW Pooled",
      Estimate = pooled_beta,
      Lower    = pooled_lower,
      Upper    = pooled_upper
    )
    
    df_plot <- df_sub %>%
      select("SNP", "Estimate", "Lower", "Upper") %>%
      bind_rows(pooled_row)
    
    snp_order <- df_plot %>% filter(SNP != "IVW Pooled") %>% arrange(Estimate) %>% pull(SNP)
    new_levels <- c("IVW Pooled", snp_order)
    df_plot <- df_plot %>% mutate(SNP = factor(SNP, levels = new_levels))
    
    # compute raw data limits and range
    all_vals <- c(df_plot$Estimate, df_plot$Lower, df_plot$Upper)
    x_min <- min(all_vals, na.rm = TRUE)
    x_max <- max(all_vals, na.rm = TRUE)
    range_x <- x_max - x_min
    
    # allow 20% extra space on right for labels
    data_limits <- c(x_min, x_max)
    label_x <- x_max
    
    # override with custom if provided
    if (!is.null(xlim_custom)) {
      if (length(xlim_custom) != 2 || !is.numeric(xlim_custom)) {
        stop("xlim_custom must be a numeric vector of length 2.")
      }
      data_limits <- xlim_custom
      label_x <- data_limits[2]
    }
    
    fmt_str <- paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]")
    p <- ggplot(df_plot, aes(x = Estimate, y = SNP)) +
      geom_point(size = dot_size) +
      geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = dot_size / 20) +
      geom_vline(xintercept = ref_line, linetype = "dashed") +
      # fixed data limits, expanded on right by 20%
      scale_x_continuous(limits = data_limits,
                         expand = expansion(mult = c(0, 0.2))) +
      coord_cartesian(clip = "off") +
      geom_text(aes(label = sprintf(fmt_str, Estimate, Lower, Upper)),
                x = label_x, hjust = 0,
                size = label_text_size) +
      labs(
        title = paste("Exposure:", current_exposure, "\nOutcome:", current_outcome),
        x     = current_form,
        y     = NULL
      ) +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text       = element_text(size = axis_text_size),
        axis.title      = element_text(size = axis_title_size),
        plot.title      = element_text(face = "bold", size = axis_title_size + 2),
        plot.margin     = margin(5, 20, 5, 5)
      )
    
    if (save_plot) {
      file_name <- paste0("SNP_forest_",
                          gsub("\\s+", "_", current_exposure), "_",
                          gsub("\\s+", "_", current_outcome),
                          ".", file_type)
      ggsave(
        filename = file.path(save_dir, file_name),
        plot     = p,
        device   = file_type,
        width    = plot_width,
        height   = plot_height,
        dpi      = 300
      )
    }
    
    plots_list[[i]] <- p
  }
  names(plots_list) <- paste(outcome_exposure$Outcome, outcome_exposure$Exposure, sep = " vs ")
  return(plots_list)
}
