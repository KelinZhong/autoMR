library(dplyr)

process_snp_data <- function(df) {
  
  df <- df %>%
    rowwise() %>%
    mutate(
      # Task 1: Check if exposure and outcome alleles match
      swap_needed = NEA_exposure != NEA_outcome | EA_exposure != EA_outcome,
      
      # Store swapped values for outcome alleles and beta
      temp_NEA_outcome = ifelse(swap_needed, EA_outcome, NEA_outcome),
      temp_EA_outcome = ifelse(swap_needed, NEA_outcome, EA_outcome),
      temp_beta_outcome = ifelse(swap_needed, -beta_outcome, beta_outcome),
      temp_A1FREQ_outcome = ifelse(swap_needed, 1 - A1FREQ_outcome, A1FREQ_outcome)
    ) %>%
    mutate(
      # Assign the swapped values simultaneously
      NEA_outcome = temp_NEA_outcome,
      EA_outcome = temp_EA_outcome,
      beta_outcome = temp_beta_outcome,
      A1FREQ_outcome = temp_A1FREQ_outcome
    ) %>%
    mutate(
      # Task 2: Flip beta exposure if negative
      flip_beta = beta_exposure < 0,
      
      # Store flipped values for exposure
      temp_NEA_exposure = ifelse(flip_beta, EA_exposure, NEA_exposure),
      temp_EA_exposure = ifelse(flip_beta, NEA_exposure, EA_exposure),
      temp_A1FREQ_exposure = ifelse(flip_beta, 1 - A1FREQ_exposure, A1FREQ_exposure),
      
      # Store flipped values for outcome (to remain consistent with exposure)
      temp_NEA_outcome_flip = ifelse(flip_beta, EA_outcome, NEA_outcome),
      temp_EA_outcome_flip = ifelse(flip_beta, NEA_outcome, EA_outcome),
      temp_A1FREQ_outcome_flip = ifelse(flip_beta, 1 - A1FREQ_outcome, A1FREQ_outcome),
      
      # Store flipped beta values
      temp_beta_exposure = ifelse(flip_beta, -beta_exposure, beta_exposure),
      temp_beta_outcome = ifelse(flip_beta, -beta_outcome, beta_outcome)
    ) %>%
    mutate(
      # Assign the flipped values simultaneously
      NEA_exposure = temp_NEA_exposure,
      EA_exposure = temp_EA_exposure,
      A1FREQ_exposure = temp_A1FREQ_exposure,
      beta_exposure = temp_beta_exposure,
      beta_outcome = temp_beta_outcome,
      
      # Apply the second allele flip to the outcome (AFTER exposure flip)
      NEA_outcome = temp_NEA_outcome_flip,
      EA_outcome = temp_EA_outcome_flip,
      A1FREQ_outcome = temp_A1FREQ_outcome_flip
    ) %>%
    ungroup() %>%
    select(-swap_needed, -flip_beta, 
           -temp_NEA_outcome, -temp_EA_outcome, -temp_beta_outcome, -temp_A1FREQ_outcome,
           -temp_NEA_exposure, -temp_EA_exposure, -temp_A1FREQ_exposure, 
           -temp_NEA_outcome_flip, -temp_EA_outcome_flip, -temp_A1FREQ_outcome_flip, 
           -temp_beta_exposure)
  
  # Create transformed SNP data
  
  if("Outcome" %in% names(df) == TRUE) {
    new_df <- df %>%
      select(Outcome, SNP, NEA_exposure, EA_exposure, A1FREQ_outcome, 
             beta_exposure, se_exposure, beta_outcome, se_outcome, Exposure) %>%
      rename(ALLELE1 = NEA_exposure, 
             ALLELE0 = EA_exposure, 
             A1FREQ = A1FREQ_outcome)
  }
  
  else {
    new_df <- df %>%
      select(SNP, NEA_exposure, EA_exposure, A1FREQ_outcome, 
             beta_exposure, se_exposure, beta_outcome, se_outcome, Exposure) %>%
      rename(ALLELE1 = NEA_exposure, 
             ALLELE0 = EA_exposure, 
             A1FREQ = A1FREQ_outcome)
  }
  
  
  return(list(check_df = df, input_df = new_df))
}
