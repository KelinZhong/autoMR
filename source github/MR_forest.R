# MR_forest.R
# ---------------------------------------------------
# Forest plot for MR methods, with robust GRIP handling and CI fixes.
# Reads columns from summary table produced by valid.output/.multi
# ---------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

MR_forest <- function(df, effect, 
                      save_plot = FALSE, file_type = "jpeg", save_dir = ".",
                      custom_xlim = NULL,
                      dot_size = 3,
                      axis_text_size = 10,
                      axis_title_size = 12,
                      pval_text_size = 3,
                      plot_width = 8,
                      plot_height = 6,
                      clamp_nonpositive = FALSE) {
  outcome_exposure <- df %>% distinct(Outcome, Exposure)
  unique_outcomes  <- unique(df$Outcome)
  if(length(effect) == 1) effect <- rep(effect, length(unique_outcomes))
  if(length(effect) != length(unique_outcomes)) {
    stop("Length of effect vector must be 1 or equal to the number of unique outcomes.")
  }
  outcome_exposure$effect <- sapply(outcome_exposure$Outcome, function(x) {
    effect[which(unique_outcomes == x)]
  })
  
  methods <- c("IVW","RAPS","Med","Egger","PRESSO","Horse","GRIP")
  plots_list <- list()
  
  for(i in seq_len(nrow(outcome_exposure))) {
    curr_out <- outcome_exposure$Outcome[i]
    curr_exp <- outcome_exposure$Exposure[i]
    curr_eff <- as.character(outcome_exposure$effect[i])
    xlab_val <- if(curr_eff == "Beta") curr_eff else if(curr_eff %in% c("OR","HR")) paste0("log(",curr_eff,")") else curr_eff
    df_row <- df %>% filter(Outcome == curr_out, Exposure == curr_exp)
    snp_label <- df_row$SNPs
    method_df <- data.frame(Method = methods, stringsAsFactors = FALSE)
    
    # helper to fetch a value with GRIP fallbacks
    .get_col <- function(row, method, what, eff){
      # what in {"Estimate","Lower","Upper","P"}
      if (eff == "Beta") {
        nm <- switch(what,
                     Estimate = if (method=="PRESSO") "Presso_Beta" else if (method=="Horse") "Horse_Beta" else if (method=="GRIP") "Grip_Beta" else paste0(method,"_Beta"),
                     Lower    = if (method=="PRESSO") "Presso_lower" else if (method=="Horse") "Horse_Lower" else if (method=="GRIP") "Grip_Lower" else paste0(method,"_Lower"),
                     Upper    = if (method=="PRESSO") "Presso_upper" else if (method=="Horse") "Horse_Upper" else if (method=="GRIP") "Grip_Upper" else paste0(method,"_Upper"),
                     P        = if (method=="Egger")  "Egger_P_value" else if (method=="PRESSO") "Presso_p" else if (method=="Horse") "Horse_P" else if (method=="GRIP") "Grip_P" else paste0(method,"_P")
        )
      } else {
        nm <- switch(what,
                     Estimate = if (method=="PRESSO") paste0("Presso_", eff) else if (method=="Horse") paste0("Horse_", eff) else if (method=="GRIP") paste0("Grip_", eff) else paste0(method,"_", eff),
                     Lower    = if (method=="PRESSO") "Presso_lower" else if (method=="Horse") "Horse_Lower" else if (method=="GRIP") "Grip_Lower" else paste0(method,"_Lower"),
                     Upper    = if (method=="PRESSO") "Presso_upper" else if (method=="Horse") "Horse_Upper" else if (method=="GRIP") "Grip_Upper" else paste0(method,"_Upper"),
                     P        = if (method=="Egger")  "Egger_P_value" else if (method=="PRESSO") "Presso_p" else if (method=="Horse") "Horse_P" else if (method=="GRIP") "Grip_P" else paste0(method,"_P")
        )
      }
      # GRIP alt naming fallback
      if (method == "GRIP" && !nm %in% names(row)) {
        alt <- sub("^Grip_", "GRIP_", nm)
        if (alt %in% names(row)) nm <- alt
      }
      # GRIP adjusted CI fallback if only *_adj present
      if (method == "GRIP" && what %in% c("Lower","Upper") && !nm %in% names(row)) {
        nm2 <- sub("Grip_Lower$", "Grip_Lower_adj", nm)
        nm2 <- sub("Grip_Upper$", "Grip_Upper_adj", nm2)
        if (nm2 %in% names(row)) nm <- nm2
      }
      if (!nm %in% names(row)) return(NA_real_)
      suppressWarnings(as.numeric(row[[nm]]))
    }
    
    # --- Estimates & CI ---
    if(curr_eff == "Beta") {
      method_df$Estimate <- sapply(methods, function(m) .get_col(df_row, m, "Estimate", curr_eff))
      method_df$Lower    <- sapply(methods, function(m) .get_col(df_row, m, "Lower",    curr_eff))
      method_df$Upper    <- sapply(methods, function(m) .get_col(df_row, m, "Upper",    curr_eff))
    } else {
      est_raw <- sapply(methods, function(m) .get_col(df_row, m, "Estimate", curr_eff))
      lo_raw  <- sapply(methods, function(m) .get_col(df_row, m, "Lower",    curr_eff))
      hi_raw  <- sapply(methods, function(m) .get_col(df_row, m, "Upper",    curr_eff))
      if (clamp_nonpositive){
        eps <- .Machine$double.eps
        est_raw[!is.na(est_raw) & est_raw <= 0] <- eps
        lo_raw [!is.na(lo_raw)  & lo_raw  <= 0] <- eps
        hi_raw [!is.na(hi_raw)  & hi_raw  <= 0] <- eps
      } else {
        bad <- which((!is.na(est_raw) & est_raw <= 0) |
                       (!is.na(lo_raw)  & lo_raw  <= 0) |
                       (!is.na(hi_raw)  & hi_raw  <= 0))
        if (length(bad)) warning("Non-positive effect-scale values found; dropped before log().")
        est_raw[bad] <- NA_real_; lo_raw[bad] <- NA_real_; hi_raw[bad] <- NA_real_
      }
      method_df$Estimate <- ifelse(is.na(est_raw), NA_real_, log(est_raw))
      method_df$Lower    <- ifelse(is.na(lo_raw ), NA_real_, log(lo_raw ))
      method_df$Upper    <- ifelse(is.na(hi_raw ), NA_real_, log(hi_raw ))
    }
    method_df$Estimate <- as.numeric(method_df$Estimate)
    method_df$Lower    <- as.numeric(method_df$Lower)
    method_df$Upper    <- as.numeric(method_df$Upper)
    
    # --- P-values ---
    method_df$p <- as.numeric(sapply(methods, function(m) .get_col(df_row, m, "P", curr_eff)))
    method_df$Method <- factor(method_df$Method, levels = methods)
    
    # reference line
    ref_line <- 0
    
    # x-limits
    all_values <- c(method_df$Estimate, method_df$Lower, method_df$Upper)
    finite_values <- all_values[is.finite(all_values) & !is.na(all_values)]
    if(length(finite_values) == 0) {
      auto_xlim <- c(-1, 1); range_x <- 2
    } else {
      if(curr_eff == "Beta") {
        max_abs <- max(abs(finite_values)); auto_xlim <- c(-max_abs, max_abs)
      } else {
        auto_xlim <- c(min(finite_values), max(finite_values))
      }
      range_x <- diff(auto_xlim); if(range_x <= 0) { range_x <- 2; auto_xlim <- c(-1, 1) }
    }
    offset <- if(curr_eff == "Beta") 0.1 * range_x else 0.2 * range_x
    
    if(!is.null(custom_xlim) && length(custom_xlim)==2 && is.numeric(custom_xlim)) {
      xlim_to_use <- custom_xlim; p_label_x <- custom_xlim[2] + 0.1 * diff(custom_xlim)
    } else {
      xlim_to_use <- auto_xlim;  p_label_x <- auto_xlim[2] + offset
    }
    final_xlim <- c(xlim_to_use[1], p_label_x + 0.05 * range_x)
    final_xlim <- as.numeric(final_xlim)
    if(length(final_xlim)!=2 || any(!is.finite(final_xlim))) { final_xlim <- c(-1, 1); p_label_x <- 1.1 }
    
    # p labels
    p_labels <- ifelse(is.na(method_df$p) | !is.finite(method_df$p), "", paste0("p=", formatC(method_df$p, format="f", digits=3)))
    p_label_x <- as.numeric(p_label_x)[1]; if(!is.finite(p_label_x)) p_label_x <- 1.1
    
    # plot
    p <- ggplot(method_df, aes(x=Estimate, y=Method)) +
      geom_point(size=dot_size, na.rm=TRUE) +
      geom_errorbarh(aes(xmin=Lower, xmax=Upper), height=0.2, linewidth=dot_size*0.1, na.rm=TRUE) +
      geom_vline(xintercept=ref_line, linetype="dashed") +
      scale_x_continuous(limits=final_xlim) +
      coord_cartesian(clip="off") +
      labs(title=paste("Exposure:",curr_exp, "Outcome:",curr_out),
           subtitle=paste("SNPs:",snp_label), x=xlab_val, y="") +
      theme_minimal() +
      theme(
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text       = element_text(size=axis_text_size),
        axis.title.x    = element_text(size=axis_title_size),
        axis.title.y    = element_text(size=axis_title_size),
        plot.title      = element_text(face="bold", size=axis_title_size+2),
        plot.subtitle   = element_text(size=axis_text_size)
      ) +
      annotate("text", x=p_label_x, y=methods, label=p_labels, hjust=0, size=pval_text_size)
    
    if(save_plot) {
      fname <- paste0(gsub("\\s+","_",curr_exp),"_",gsub("\\s+","_",curr_out),".",file_type)
      ggsave(filename=file.path(save_dir,fname), plot=p, device=file_type, width=plot_width, height=plot_height, dpi=300)
    }
    plots_list[[i]] <- p
  }
  names(plots_list) <- paste(outcome_exposure$Outcome, outcome_exposure$Exposure, sep=" - ")
  return(plots_list)
}
