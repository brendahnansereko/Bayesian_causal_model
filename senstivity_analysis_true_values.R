library(rbmi)
library(rstan)
library(dplyr)
library(purrr)
library(mmrm)
#_____________________________________________
# simulate complete data set
#______________________________________________
#set.seed(423456)

samples <- 5000
sigma_k0 <- seq(.0, 1.5, .1) 
iter <- 5000
mar_mean <- j2r_mean <- mmrm_mean  <- cir_mean <- cr_mean  <-  true_value <- 0
j2r_se <- mmrm_se <- cir_se <- cr_se <- 0

j2r_constant_k0  <- cir_constant_k1 <- j2r_constant_k0_se  <- cir_constant_k1_se  <- 0


prop_check <- 0
causal_k_mean0  <- causal_k_mean05 <- causal_k_mean1 <- matrix(,nrow = iter, ncol = length(sigma_k0))
causal_k_mean0_se  <- causal_k_mean05_se <- causal_k_mean1_se <- matrix(,nrow = iter, ncol = length(sigma_k0))
causal_k_mean0_true <- causal_k_mean1_true <- causal_k_mean05_true <- matrix(,nrow = iter, ncol = length(sigma_k0))

k <- 0

for (j in 1:iter) {
  n <- 200
  
  # Visits
  visits  <-  c(0, 4, 8, 14, 20, 26)
  n_visit <-  length(visits)
  
  # Mean trajectory control
  
  muC <- c(7.92, 7.82, 7.8, 7.8, 7.78, 7.78)
  
  
  # Mean trajectory intervention
  muT <-  c(7.92, 7.55, 7.20, 7.1, 7.05, 7.05)
  
  
  # create Sigma
  
  # Spatial power correlation matrix - corresponds to AR(1) for equidistant visits
  rho      <- 0.8
  exponent <- abs(matrix(visits, nrow = n_visit, ncol = n_visit, byrow = TRUE) - visits) / (visits[2] - visits[1])
  corr_mat <- rho^exponent
  
  
  # Variance at baseline is pooled variance from each group
  sigma2_bl = 0.48 
  
  # Variance for treatment and control group at each visit
  sigma2_trt = c(sigma2_bl, 0.8, 1.1, 1.4, 1.23, 1.48)
  sigma2_ctl = c(sigma2_bl, 0.75, 0.80, 0.90, 1.06, 1.14)  
  # sigma2_ctl = c(sigma2_bl, 0.8, 1.1, 1.4, 1.23, 1.48)
  
  # Percentage of pts. that discontinued randomized trt
  rate_trtdiscont_trt <- (175 - c(175, 169, 164, 156, 151, 149)) / 175 
  rate_trtdiscont_ctl <- (178 - c(178, 172, 165, 158, 140, 133)) / 178 
  
  
  # covariance matrix for control group
  mat1 <- matrix(sigma2_ctl, nrow = n_visit, ncol = n_visit, byrow = TRUE)
  mat2 <-  matrix(sigma2_ctl, nrow = n_visit, ncol = n_visit, byrow = FALSE)
  Sigma_ctl <- sqrt(mat1 * mat2) * corr_mat
  
  
  # covariance matrix for treatment group
  mat1 <- matrix(sigma2_trt, nrow = n_visit, ncol = n_visit, byrow = TRUE)
  mat2 <-  matrix(sigma2_trt, nrow = n_visit, ncol = n_visit, byrow = FALSE)
  Sigma_trt <- sqrt(mat1 * mat2) * corr_mat
  
  
  # Set probability of discontinuation
  probDisc_C <- 0.02
  probDisc_T <- 0.03
  or_outcome <- 1.10 # +1 point increase => +10% odds of discontinuation
  
  # Set drop-out rate following discontinuation
  prob_dropout <- 0.5
  
  # Set simulation parameters of the control group
  parsC <- set_simul_pars(
    mu = muC,
    sigma = Sigma_ctl,
    n = n)
  
  # Set simulation parameters of the intervention group
  # Set simulation parameters of the control group
  parsT <- set_simul_pars(
    mu = muT,
    sigma = Sigma_trt,
    n = n)
  
  # Simulate data
  data <- simulate_data(
    pars_c = parsC,
    pars_t = parsT,
    post_ice1_traj = "MAR")
  
  data <- subset(data, select = c("id","visit","group","outcome_bl","outcome_noICE"))
  colnames(data) <- c("id","visit","group","y_bl","y")
  data$trt <- ifelse(data$group=="Control",0,1)
  data$groupx <- ifelse(data$group=="Control","Control","Treatment")
  
  # Beta parameter for logistic regression modeling of DAR trt discont
  beta_discont_dar <- data.frame(groupx = rep(c("Control", "Treatment"), each = 5),
                                 visitn = rep(1:5, times = 2),
                                 beta_main = rep(-13, times = 10),
                                 beta_prev_resp = c(1.42, 1.14, 1.33, 1.51, 1.46, 
                                                    1.42, 1.14, 1.47, 1.48, 1.40), 
                                 beta_current_resp = rep(0, times = 10),
                                 beta_baseline = c(0, 0.3, 0.1, 0.05, 0, 
                                                   0, 0.3, 0.1, 0.05, 0), 
                                 stringsAsFactors = FALSE)
  
  #_____________________________________________
  # simulate Discontinuation
  #______________________________________________
  is_dnar <- any(names(beta_discont_dar) == "beta_current_resp")
  
  data$visitn <- as.numeric(data$visit) -1
  data$baseline_var <- data$y_bl
  data$response_ontrt <- data$y
  
  data_discont <- full_join(data, beta_discont_dar, 
                            by = c("groupx", "visitn")) 
  
  data_discont <- data_discont %>% 
    group_by(id) %>% 
    mutate(prev_response = lag(response_ontrt, default = 0),
           logit_prob = beta_main + beta_prev_resp * prev_response + beta_baseline * baseline_var,
           logit_prob = if(is_dnar){logit_prob + beta_current_resp * response_ontrt} else {logit_prob},
           prob = 1/(1+exp(-logit_prob)),
           prob = if_else(is.na(prob) | visitn==1 , 0, prob),
           ontrt = 1-rbinom(n = n(), size = 1, prob = prob)) %>%  
    mutate(first_offtrt_visit = min(visitn[ontrt == 0], max(visitn)),
           ontrt = if_else(visitn > first_offtrt_visit, 0, ontrt),
           r = ifelse(max(first_offtrt_visit)==5 & min(ontrt)==1,0,first_offtrt_visit-1)) %>% 
    ungroup()
  
  data_discont <- data_discont %>% 
    group_by(id) %>% 
    mutate(Discontinue = if_else(any(ontrt == 0), 1, 0)) %>%
    ungroup()
  
  data_true <- data_discont %>% 
    mutate(
      y_mar = if_else(ontrt == 1 | (ontrt == 0 & groupx == "Control"), response_ontrt, as.numeric(NA))
    ) 
  
  
  # _______________________________________
  # 1. Conditional imputation (Σ uses trt variances before r_time, ctrl after)
  # _______________________________________
  
  # helper: build subject-specific mixed covariance given r_time
  build_sigma_mixed <- function(r_time, sigma2_trt, sigma2_ctl, corr_mat) {
    # visits are indexed 1..6 (0,4,8,14,20,26 weeks)
    # r_time is in 0..5 (your code uses visitn = visit - 1)
    v_idx <- 1:6
    # SDs per visit: trt up to and including r_time+1 (the last on-trt visit),
    # then control afterwards
    sd_vec <- sqrt(ifelse(v_idx <= (r_time+1), sigma2_trt[v_idx], sigma2_ctl[v_idx]))
    # covariance with shared correlations
    (sd_vec %o% sd_vec) * corr_mat
  }
  
  for (i in unique(data_true$id[data_true$groupx == "Treatment" & data_true$Discontinue == 1])) {
    person_data <- data_true[data_true$id == i, ]
    r_time <- max(person_data$r)
    
    # r_time is 0..5; if r_time==0 (no post-baseline obs) or r_time>=5 (no future), skip
    if (r_time < 5 && r_time != 0) {
      # observed indices are post-baseline up to r_time+1 (because visitn = 0..5 maps to 1..6)
      obs_idx <- 2:(r_time + 1)
      fut_idx <- (r_time + 2):6
      
      # build Σ that mixes trt and ctrl variances around r_time
      Sigma_mix <- build_sigma_mixed(
        r_time = r_time,
        sigma2_trt = sigma2_trt,
        sigma2_ctl = sigma2_ctl,
        corr_mat = corr_mat
      )
      
      
      
      # conditional blocks
      Sigma_11 <- as.matrix(Sigma_mix[obs_idx, obs_idx, drop = FALSE])
      Sigma_22 <- as.matrix(Sigma_mix[fut_idx, fut_idx, drop = FALSE])
      Sigma_21 <- as.matrix(Sigma_mix[fut_idx, obs_idx, drop = FALSE])
      Sigma_12 <- t(Sigma_21)
      
      beta <- Sigma_21 %*% solve(Sigma_11)
      Sigma_cond <- Sigma_22 - beta %*% Sigma_12
      
      # observed outcomes and means
      Y_obs <- as.numeric(person_data$y_mar[person_data$visitn %in% (obs_idx - 1)])
      mu_obs_trt     <- muT[obs_idx]
      mu_future_trt  <- muT[fut_idx]
      mu_future_ctrl <- muC[fut_idx]
      
      # mean after discontinuation: migrate from trt path to control path
      # (J2R-like: future mean = control, with conditioning on on-trt residuals)
      mu_cond <- as.vector(beta %*% (Y_obs - mu_obs_trt)) + mu_future_ctrl
      mu_cond <- mu_cond + k * (muT[r_time+1] - muC[r_time+1]) 
      
      # draw imputed future values
      imputed_vals <- MASS::mvrnorm(1, mu = mu_cond, Sigma = Sigma_cond)
      
      # write back
      for (p in seq_along(fut_idx)) {
        data_true$y_mar[data_true$id == i & data_true$visitn == (fut_idx[p] - 1)] <- imputed_vals[p]
      }
    }
  }
  
  data_true <- data_true[data_true$visitn==5,]
  lm <- summary(lm(y_mar ~ as.factor(group) , data = data_true))
  
  true_value[j] <- lm$coefficients[2, 1]
  
  print(j)
}
true_value_00_00 <- true_value
#Data sets
#res <- as.data.frame( cbind(true_value))


