# mr_scatter_doMR.R
# ---------------------------------------------------
# Egger-style scatter with fitted lines for multiple methods.
# NOW supports pulling numbers directly from the summary table (for exact consistency).
# ---------------------------------------------------
library(MendelianRandomization)
library(TwoSampleMR)
library(mr.raps)
library(MRPRESSO)

doMR <- function(d.x, d.y, d.x.se, d.y.se,
                 corr.mat,
                 n.it = 1000,
                 d.title, d.subtitle, d.xlab, d.ylab,
                 d.snpinfo,
                 d.out = NA,
                 methods.plot = c('IVW','RAPS','Egger','PRESSO','Horse','GRIP'),
                 NbDistribution_presso = 1000,
                 SignifThreshold_presso = 0.05,
                 mr_horse_n_iter = 5000,
                 mr_horse_n_burnin = 1000,
                 show.legend = TRUE,
                 df_row = NULL,                 # <- row from summary table (same as forest)
                 effect_scale = "Beta",         # 'Beta','OR','HR' (needed if using df_row)
                 use_df_results = TRUE) {       # <- when TRUE, read betas/p from df_row
  # determine file name
  suffix <- if (show.legend) "" else "_no_legend"
  if (!is.na(d.out)) file_name <- paste0(d.out, suffix, ".png")
  
  # axis limits based on raw points ± SE
  x_ext <- max(abs(d.x - d.x.se), abs(d.x + d.x.se), na.rm = TRUE)
  y_ext <- max(abs(d.y - d.y.se), abs(d.y + d.y.se), na.rm = TRUE)
  limx <- c(0, x_ext)
  limy <- c(-y_ext, y_ext)
  
  # plotting base
  if (!is.na(d.out)) png(file_name, width = 8, height = 6, units = "in", res = 300)
  plot(d.x, d.y, xlim = limx, ylim = limy,
       xlab = d.xlab, ylab = d.ylab,
       main = d.title, pch = 16, bty = "L")
  mtext(d.subtitle)
  abline(h = 0, lty = 2); abline(v = 0, lty = 2)
  arrows(d.x, d.y - d.y.se, d.x, d.y + d.y.se, length = 0, angle = 90, code = 3, col = "grey")
  arrows(d.x - d.x.se, d.y, d.x + d.x.se, d.y, length = 0, angle = 90, code = 3, col = "grey")
  
  # MR input for recomputation fallback
  d.mr <- mr_input(
    bx = d.x, bxse = d.x.se,
    by = d.y, byse = d.y.se,
    correlation   = corr.mat,
    exposure      = d.xlab,
    outcome       = d.ylab,
    snps          = d.snpinfo$SNP,
    effect_allele = d.snpinfo$ALLELE0,
    other_allele  = d.snpinfo$ALLELE1,
    eaf           = d.snpinfo$A1FREQ
  )
  all_res <- if (!use_df_results) MendelianRandomization::mr_allmethods(d.mr, iterations = n.it)$Values else NULL
  
  colors <- list(IVW='chartreuse3', RAPS='turquoise3', Egger='cornflowerblue',
                 PRESSO='red', Horse='purple', GRIP='darkorange')
  legend_txt  <- character(); legend_cols <- character(); legend_lty <- integer()
  
  # helpers to pull from df_row (same mapping as MR_forest)
  .from_df_beta <- function(m){
    stopifnot(!is.null(df_row))
    if (effect_scale == "Beta"){
      col <- if (m=="PRESSO") "Presso_Beta" else if (m=="Horse") "Horse_Beta" else if (m=="GRIP") "Grip_Beta" else paste0(m,"_Beta")
      return(suppressWarnings(as.numeric(df_row[[col]])))
    } else {
      col <- if (m=="PRESSO") paste0("Presso_",effect_scale) else if (m=="Horse") paste0("Horse_",effect_scale) else if (m=="GRIP") paste0("Grip_",effect_scale) else paste0(m,"_",effect_scale)
      val <- suppressWarnings(as.numeric(df_row[[col]]))
      if (is.na(val) || !is.finite(val) || val <= 0) return(NA_real_)
      return(log(val))
    }
  }
  .from_df_p <- function(m){
    stopifnot(!is.null(df_row))
    col <- if (m=="Egger") "Egger_P_value" else if (m=="PRESSO") "Presso_p" else if (m=="Horse") "Horse_P" else if (m=="GRIP") "Grip_P" else paste0(m,"_P")
    suppressWarnings(as.numeric(df_row[[col]]))
  }
  
  add_line <- function(label, m, lty_code){
    # slope
    b <- if (use_df_results && !is.null(df_row)){
      .from_df_beta(m)
    } else {
      switch(m,
             IVW = subset(all_res, Method=='IVW')$Estimate,
             RAPS = {
               rr <- try(mr.raps(d.x, d.y, d.x.se, d.y.se), TRUE)
               if (inherits(rr, "try-error")) NA_real_ else rr$beta.hat
             },
             Egger = { ee <- MendelianRandomization::mr_egger(d.mr); ee@Estimate },
             PRESSO = {
               if (length(d.x) < 4) NA_real_ else {
                 pr <- mr_presso(
                   BetaOutcome='beta_outcome', BetaExposure='beta_exposure',
                   SdOutcome='se_outcome',     SdExposure='se_exposure',
                   data=data.frame(beta_exposure=d.x, se_exposure=d.x.se,
                                   beta_outcome=d.y,  se_outcome=d.y.se,
                                   SNP=d.snpinfo$SNP),
                   OUTLIERtest=TRUE, DISTORTIONtest=TRUE,
                   NbDistribution=NbDistribution_presso, SignifThreshold=SignifThreshold_presso
                 )
                 as.numeric(pr$`Main MR results`[1,'Causal Estimate'])
               }
             },
             Horse = {
               H <- mr_horse(data.frame(betaY=d.y, betaX=d.x, betaYse=d.y.se, betaXse=d.x.se),
                             no_ini=3, variable.names='theta',
                             n.iter=mr_horse_n_iter, n.burnin=mr_horse_n_burnin)
               H$MR_Estimate$Estimate[1]
             },
             GRIP = {
               gg <- try(TwoSampleMR::mr_grip(b_exp=d.x, b_out=d.y, se_exp=d.x.se, se_out=d.y.se,
                                              parameters=try(TwoSampleMR::default_parameters(), TRUE)), TRUE)
               if (inherits(gg,"try-error")) NA_real_ else as.numeric(gg$b)
             }
      )
    }
    # p-value
    p <- if (use_df_results && !is.null(df_row)){
      .from_df_p(m)
    } else {
      switch(m,
             IVW = subset(all_res, Method=='IVW')$`P-value`,
             RAPS = {
               rr <- try(mr.raps(d.x, d.y, d.x.se, d.y.se), TRUE)
               if (inherits(rr,"try-error")) NA_real_ else rr$beta.p.value
             },
             Egger  = MendelianRandomization::mr_egger(d.mr)@Pvalue.Est,
             PRESSO = {
               if (length(d.x) < 4) NA_real_ else {
                 pr <- mr_presso(
                   BetaOutcome='beta_outcome', BetaExposure='beta_exposure',
                   SdOutcome='se_outcome',     SdExposure='se_exposure',
                   data=data.frame(beta_exposure=d.x, se_exposure=d.x.se,
                                   beta_outcome=d.y,  se_outcome=d.y.se,
                                   SNP=d.snpinfo$SNP),
                   OUTLIERtest=TRUE, DISTORTIONtest=TRUE,
                   NbDistribution=NbDistribution_presso, SignifThreshold=SignifThreshold_presso
                 )
                 as.numeric(pr$`Main MR results`[1,'P-value'])
               }
             },
             Horse = {
               H <- mr_horse(data.frame(betaY=d.y, betaX=d.x, betaYse=d.y.se, betaXse=d.x.se),
                             no_ini=3, variable.names='theta',
                             n.iter=mr_horse_n_iter, n.burnin=mr_horse_n_burnin)
               bH <- H$MR_Estimate$Estimate[1]; seH <- H$MR_Estimate$SD[1]
               2*(1 - pnorm(abs(bH/seH)))
             },
             GRIP = {
               gg <- try(TwoSampleMR::mr_grip(b_exp=d.x, b_out=d.y, se_exp=d.x.se, se_out=d.y.se,
                                              parameters=try(TwoSampleMR::default_parameters(), TRUE)), TRUE)
               if (inherits(gg,"try-error")) NA_real_ else as.numeric(gg$pval)
             }
      )
    }
    if (is.finite(b)){
      abline(a=0, b=b, lty=lty_code, lwd=2, col=colors[[m]])
      legend_txt  <<- c(legend_txt, sprintf('%s: beta=%.3f, p=%.3g', label, b, p))
    } else {
      legend_txt  <<- c(legend_txt, sprintf('%s: failed', label))
    }
    legend_cols <<- c(legend_cols, colors[[m]])
    legend_lty  <<- c(legend_lty, lty_code)
  }
  
  if ('IVW'    %in% methods.plot) add_line('IVW',    'IVW',    3)
  if ('RAPS'   %in% methods.plot) add_line('RAPS',   'RAPS',   3)
  if ('Egger'  %in% methods.plot) add_line('Egger',  'Egger',  4)
  if ('PRESSO' %in% methods.plot) add_line('PRESSO', 'PRESSO', 5)
  if ('Horse'  %in% methods.plot) add_line('Horse',  'Horse',  6)
  if ('GRIP'   %in% methods.plot) add_line('GRIP',   'GRIP',   2)
  
  if (show.legend) legend('topright', legend=legend_txt, col=legend_cols, lty=legend_lty, lwd=2, bg='white', cex=0.9)
  if (!is.na(d.out) && dev.cur() != 1) dev.off()
}

# wrappers that pass the summary table row for exact consistency with the forest/table
MRplots <- function(MR_input_data,
                    d.title = NULL,
                    d.subtitle = NULL,
                    plot.xlab,
                    plot.ylab,
                    outcome_label,
                    methods.plot = c('IVW','RAPS','Egger','PRESSO','Horse','GRIP'),
                    NbDistribution_presso = 1000,
                    SignifThreshold_presso= 0.05,
                    mr_horse_n_iter       = 5000,
                    mr_horse_n_burnin     = 1000,
                    show.legend           = TRUE,
                    summary_df            = NULL,    # <- combined results from valid.output/.multi
                    effect_scale          = "Beta",  # "Beta","OR","HR" to match summary_df
                    use_df_results        = TRUE) {  # <- read directly from summary_df
  for (ex in unique(MR_input_data$Exposure)) {
    dat <- subset(MR_input_data, Exposure == ex)
    df_row <- NULL
    if(!is.null(summary_df)){
      hit <- subset(summary_df, Outcome == outcome_label & Exposure == ex)
      if(nrow(hit) > 0) df_row <- hit[1, , drop = FALSE]
    }
    doMR(
      d.x = dat$beta_exposure, d.y = dat$beta_outcome,
      d.x.se = dat$se_exposure, d.y.se = dat$se_outcome,
      corr.mat = diag(1, nrow(dat)),
      d.title = d.title, d.subtitle = d.subtitle,
      d.xlab = paste(plot.xlab, ex), d.ylab = paste(plot.ylab, outcome_label),
      d.snpinfo = dat[, c('SNP','ALLELE1','ALLELE0','A1FREQ')],
      d.out = paste0('Exp_', ex, '_Out_', outcome_label),
      show.legend = show.legend,
      methods.plot = methods.plot,
      NbDistribution_presso = NbDistribution_presso,
      SignifThreshold_presso= SignifThreshold_presso,
      mr_horse_n_iter = mr_horse_n_iter, mr_horse_n_burnin = mr_horse_n_burnin,
      df_row = df_row, effect_scale = effect_scale, use_df_results = use_df_results
    )
  }
}

MRplots.multi <- function(MR_input_data,
                          plot.xlab,
                          plot.ylab,
                          methods.plot          = c('IVW','RAPS','Egger','PRESSO','Horse','GRIP'),
                          NbDistribution_presso = 1000,
                          SignifThreshold_presso= 0.05,
                          mr_horse_n_iter       = 5000,
                          mr_horse_n_burnin     = 1000,
                          show.legend           = TRUE,
                          summary_df            = NULL,
                          effect_scale          = "Beta",
                          use_df_results        = TRUE) {
  for (out in unique(MR_input_data$Outcome)) {
    sub <- subset(MR_input_data, Outcome == out)
    MRplots(sub,
            d.title = NULL, d.subtitle = NULL,
            plot.xlab, plot.ylab,
            outcome_label = out,
            methods.plot = methods.plot,
            NbDistribution_presso = NbDistribution_presso,
            SignifThreshold_presso= SignifThreshold_presso,
            mr_horse_n_iter = mr_horse_n_iter, mr_horse_n_burnin = mr_horse_n_burnin,
            show.legend = show.legend,
            summary_df = summary_df,
            effect_scale = effect_scale,
            use_df_results = use_df_results)
  }
}
