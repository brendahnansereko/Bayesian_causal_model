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
iter <- 1
mar_mean <- j2r_mean <- mmrm_mean  <- cir_mean <- cr_mean  <-  0
j2r_se <- mmrm_se <- cir_se <- cr_se <- 0

j2r_constant_k0  <- cir_constant_k1 <- j2r_constant_k0_se  <- cir_constant_k1_se  <- 0


prop_check <- 0
causal_k_mean0  <- causal_k_mean05 <- causal_k_mean1 <- matrix(,nrow = iter, ncol = length(sigma_k0))
causal_k_mean0_se  <- causal_k_mean05_se <- causal_k_mean1_se <- matrix(,nrow = iter, ncol = length(sigma_k0))
causal_k_mean0_true <- causal_k_mean1_true <- causal_k_mean05_true <- matrix(,nrow = iter, ncol = length(sigma_k0))

for (j in 1:iter) {
  n <- 100
  
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
  sigma2_ctl = c(sigma2_bl, 0.8, 1.1, 1.4, 1.23, 1.48)
  
  
  # Percentage of pts. that discontinued randomized trt
  rate_trtdiscont_trt <- (175 - c(175, 169, 164, 156, 151, 149)) / 175 
  rate_trtdiscont_ctl <- (178 - c(178, 172, 165, 158, 140, 133)) / 178 
  
  
  # covariance matrix for treatment group
  mat1 <- matrix(sigma2_ctl, nrow = n_visit, ncol = n_visit, byrow = TRUE)
  mat2 <-  matrix(sigma2_ctl, nrow = n_visit, ncol = n_visit, byrow = FALSE)
  Sigma <- sqrt(mat1 * mat2) * corr_mat
  
  
  # Set probability of discontinuation
  probDisc_C <- 0.02
  probDisc_T <- 0.03
  or_outcome <- 1.10 # +1 point increase => +10% odds of discontinuation
  
  # Set drop-out rate following discontinuation
  prob_dropout <- 0.5
  
  # Set simulation parameters of the control group
  parsC <- set_simul_pars(
    mu = muC,
    sigma = Sigma,
    n = n)
  
  # Set simulation parameters of the intervention group
  # Set simulation parameters of the control group
  parsT <- set_simul_pars(
    mu = muT,
    sigma = Sigma,
    n = n)
  
  # Simulate data
  data <- simulate_data(
    pars_c = parsC,
    pars_t = parsT,
    post_ice1_traj = "MAR")
  
  data <- subset(data, select = c("id","visit","group","outcome_bl","outcome_noICE"))
  colnames(data) <- c("id","visit","group","y_bl","y")
  data$trt <- ifelse(data$group=="Control",0,1)
  
  data$groupx <- ifelse(data$group=="Control","Control","Intervention")
  
  # Beta parameter for logistic regression modeling of DAR trt discont
  beta_discont_dar <- data.frame(groupx = rep(c("Control", "Intervention"), each = 5),
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
           r = ifelse(max(first_offtrt_visit)==5 & min(ontrt)==1,0,first_offtrt_visit)) %>%
    ungroup()
  
  # Genarate discontinue patterns
  # data_discont$r <- if_else(data_discont$first_offtrt_visit==5 & 
  # data_discont$ontrt==1,0,data_discont$first_offtrt_visit)
  
  # generate missing values after discontinuation
  data_discont$y_mar <- if_else(data_discont$ontrt==1 ,data_discont$response_ontrt,NA)
  
  data_discont$visit <- as.numeric(data_discont$visit) -1 #generate appropriate factor levels
  
  data_discont <- data_discont[data_discont$visitn>0,]
  
  data <- data_discont
  #____________________________________________________
  #        MAR imputation analysis
  #____________________________________________________
  
  data$visit <- as.factor(data$visit)
  
  #expand the data by adding rows with NAs for the visits that patients were missing:
  expandedData <- expand_locf(
    data,
    id = levels(data$id),
    visit = levels(data$visit),
    vars = c("y_bl","group"),
    group = c("id"),
    order = c("id", "visit")
  )
  
  # Create another data frame that records, for each patient,
  # the first visit after the  (ICE) occurred and the imputation approach for missing data after
  # the ICE.
  ice_info <- expandedData %>%
    arrange(id, visit) %>%
    filter(is.na(y_mar)) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(id, visit) %>%
    mutate(strategy = "MAR")
  
  # Impute and analyse - The following code is for using Stan to fit the Bayesian model:
  bayesDraws <- draws(
    data = expandedData,
    data_ice = ice_info,
    vars = set_vars(
      outcome = "y_mar",
      visit = "visit",
      subjid = "id",
      group = "group",
      covariates = c("y_bl*visit", "group*visit")),
    method = method_bayes(control=control_bayes(warmup=200,thin=1),
                          n_samples = 5000,
                          same_cov = TRUE
                          #seed = 484381136
    ))
  
  # now thin the draws for the imputation methods
  thinnedDraws <- bayesDraws
  # take every 50th draw
  thinnedDraws$samples <- thinnedDraws$samples[seq(50, length(thinnedDraws$samples), by = 50)]
  # now modify method accordingly as if we had thinned in the first place
  thinnedDraws$method$n_samples <- length(thinnedDraws$samples)
  thinnedDraws$method$control$thin <- 50
  class(thinnedDraws$samples) <- c("sample_list", "list")
  
  #____________________________________________________
  #        RBI imputation analysis - JR
  #____________________________________________________
  
  #set.seed(837263)
  ice_info_jr <- ice_info %>%
    mutate(strategy = ifelse(strategy == "MAR", "JR", strategy))
  
  jrImps <- impute(thinnedDraws,
                   references = c("Intervention" = "Control", "Control" = "Control"),
                   update_strategy = ice_info_jr)
  ancovaVars <- set_vars(
    subjid = "id",
    outcome = "y_mar",
    visit = "visit",
    group = "group",
    covariates = c("y_bl"))
  jrMIFits <- analyse(
    imputations = jrImps,
    vars = ancovaVars
  )
  jrMIResults <- pool(jrMIFits)
  jrMIResults
  j2r_mean[j] <- jrMIResults$pars$trt_5$est
  j2r_se[j] <- jrMIResults$pars$trt_5$se
  
  #____________________________________________________
  #        RBI imputation analysis - CIR
  #____________________________________________________
  
  ice_info_cir <- ice_info_jr %>%
    mutate(strategy = ifelse(strategy == "JR", "CIR", strategy))
  
  #set.seed(837263)
  cirImps <- impute(thinnedDraws,
                    references = c("Intervention" = "Control", "Control" = "Control"),
                    update_strategy = ice_info_cir)
  
  cirMIFits <- analyse(
    imputations = cirImps,
    vars = ancovaVars
  )
  cirMIResults <- pool(cirMIFits)
  cirMIResults
  cir_mean[j] <- cirMIResults$pars$trt_5$est
  cir_se[j] <- cirMIResults$pars$trt_5$se

  
  #____________________________________________________
  #        MRMM
  #____________________________________________________
  mmrmFit <- mmrm(y_mar~y_bl*visit+group*visit+us(visit|id), data=data,
                  reml=TRUE, method = "Kenward-Roger")
  
  # Obtain estimates
  t_lastvisit <- numeric(length(mmrmFit$beta_est))
  
  #set elements corresponding to main effect of drug and visit 7 interaction to 1
  t_lastvisit[match("groupIntervention", names(mmrmFit$beta_est))] <- 1
  t_lastvisit[match("visit5:groupIntervention", names(mmrmFit$beta_est))] <- 1
  mmrm <- df_1d(mmrmFit, t_lastvisit)
  mmrm
  
  mmrm_mean[j] <- mmrm$est
  mmrm_se[j] <- mmrm$se
  
  #______________________________________
  #  Bayesian Approach Analysis by Liu
  #______________________________________
  
  # Obtain the Bayesian MMRM paramters from the rbmi package
  betas <- as.data.frame(do.call(rbind, do.call(rbind, bayesDraws[["samples"]])[,3]))
  
  # Obtain proportions for ICE at each visit
  a <- c(0,2,3,4,5)
  
  data_x <- subset(data,data$group=="Intervention")
  
  n_prop <- (table(factor(data_x$r, levels=a)))/5
  
  
  # Draws for proportions
  prop_sampled <- rBeta2009::rdirichlet(samples, (rep(1,5)) + n_prop)
  
  #J2R
  betas$j2r_constant_k0_values <- (betas$V2 + betas$V15)*prop_sampled [,1]
  j2r_constant_k0[j] <- mean(betas$j2r_constant_k0_values)
  j2r_constant_k0_se[j] <- sd(betas$j2r_constant_k0_values)
  
  #CIR
  betas$cir_constant_k1_values <- (betas$V2 + betas$V15)*prop_sampled [,1] + ((betas$V2)*prop_sampled [,2]) + 
    ((betas$V2 + betas$V12)*prop_sampled [,3]) + ((betas$V2 + betas$V13)*prop_sampled[,4]) +
    ((betas$V2 + betas$V14)*prop_sampled[,5])
  cir_constant_k1[j] <- mean(betas$cir_constant_k1_values)
  cir_constant_k1_se[j] <- sd(betas$cir_constant_k1_values)
  
  
  
  #__________________________________________
  #  Casual  Analysis
  #______________________________________
  
  for (i in 1:length(sigma_k0)) {
    #J2R
    k <- rnorm(samples,0,sigma_k0[i])
    betas$causal_k_mean0 <-  (betas$V2 + betas$V15)*prop_sampled [,1] + k*((betas$V2)*prop_sampled [,2]) + 
      k*((betas$V2 + betas$V12)*prop_sampled [,3]) + k*((betas$V2 + betas$V13)*prop_sampled[,4]) +
      k*((betas$V2 + betas$V14)*prop_sampled[,5])
    
    causal_k_mean0[j,i] <- mean(betas$causal_k_mean0)
    causal_k_mean0_se[j,i] <- sd(betas$causal_k_mean0)
    
    #CIR
    k_2 <- rnorm(samples,1,sigma_k0[i])
    betas$causal_k_mean1 <-  (betas$V2 + betas$V15)*prop_sampled [,1] + k_2*((betas$V2)*prop_sampled [,2]) + 
      k_2*((betas$V2 + betas$V12)*prop_sampled [,3]) + k_2*((betas$V2 + betas$V13)*prop_sampled[,4]) +
      k_2*((betas$V2 + betas$V14)*prop_sampled[,5])
    
    causal_k_mean1[j,i] <- mean(betas$causal_k_mean1)
    causal_k_mean1_se[j,i] <- sd(betas$causal_k_mean1)
    
    
    k_3 <- rnorm(samples,0.5,sigma_k0[i])
    betas$causal_k_mean05 <-  (betas$V2 + betas$V15)*prop_sampled [,1] + k_3*((betas$V2)*prop_sampled [,2]) + 
      k_3*((betas$V2 + betas$V12)*prop_sampled [,3]) + k_3*((betas$V2 + betas$V13)*prop_sampled[,4]) +
      k_3*((betas$V2 + betas$V14)*prop_sampled[,5])
    
    causal_k_mean05[j,i] <- mean(betas$causal_k_mean05)
    causal_k_mean05_se[j,i] <- sd(betas$causal_k_mean05)
    
    
    
    #TRue values
    prop_true <- c(0.53117, 0.16526, 0.17739, 0.09409, 0.03209) 
    
    
    
    # Draw for J2R true value
    k_true_1 <- rnorm(1,0,sigma_k0[i])
    causal_k_mean0_true[j,i] <-  ((muT[6] - muC[6])*prop_true [1]) + k_true_1*((muT[2] - muC[2])*prop_true [2]) + k_true_1*((muT[3] - muC[3])*prop_true [3])
    + k_true_1*((muT[4] - muC[4])*prop_true [4])  + k_true_1*((muT[5] - muC[5])*prop_true [5])
    
    
    # Draw for CIR true value
    k_true_2 <- rnorm(1,1,sigma_k0[i])
    causal_k_mean1_true[j,i] <- ((muT[6] - muC[6])*prop_true [1]) + k_true_2*((muT[2] - muC[2])*prop_true [2]) + k_true_2*((muT[3] - muC[3])*prop_true [3])
    + k_true_2*((muT[4] - muC[4])*prop_true [4])  + k_true_2*((muT[5] - muC[5])*prop_true [5])
    
    
    # Draw for CIR true value
    k_true_3 <- rnorm(1,1,sigma_k0[i])
    causal_k_mean05_true[j,i] <- ((muT[6] - muC[6])*prop_true [1]) + k_true_3*((muT[2] - muC[2])*prop_true [2]) + k_true_3*((muT[3] - muC[3])*prop_true [3])
    + k_true_3*((muT[4] - muC[4])*prop_true [4])  + k_true_3*((muT[5] - muC[5])*prop_true [5])
    
  }
  
}

#Data sets
res <- as.data.frame( cbind(j2r_constant_k0,cir_constant_k1,j2r_constant_k0_se,cir_constant_k1_se,
                            mmrm_mean,mmrm_se,cir_mean,cir_se,j2r_mean,j2r_se))

causal_k_mean0 <- as.data.frame(causal_k_mean0)
causal_k_mean05 <- as.data.frame(causal_k_mean05)
causal_k_mean1 <- as.data.frame(causal_k_mean1)

causal_k_mean0_se <- as.data.frame(causal_k_mean0_se)
causal_k_mean05_se <- as.data.frame(causal_k_mean05_se)
causal_k_mean1_se <- as.data.frame(causal_k_mean1_se)

causal_k_mean0_true <- as.data.frame(causal_k_mean0_true)
causal_k_mean1_true <- as.data.frame(causal_k_mean1_true)
causal_k_mean05_true <- as.data.frame(causal_k_mean05_true)

write.csv(res,"/home/lsh1901704/res_low")

write.csv(causal_k_mean0,"/home/lsh1901704/causal_k_mean0_low")
write.csv(causal_k_mean05,"/home/lsh1901704/causal_k_mean05_low")
write.csv(causal_k_mean1,"/home/lsh1901704/causal_k_mean1_low")

write.csv(causal_k_mean0_se,"/home/lsh1901704/causal_k_mean0_se_low")
write.csv(causal_k_mean05_se,"/home/lsh1901704/causal_k_mean05_se_low")
write.csv(causal_k_mean1_se,"/home/lsh1901704/causal_k_mean1_se_low")

write.csv(causal_k_mean0_true,"/home/lsh1901704/causal_k_mean0_true_low")
write.csv(causal_k_mean1_true,"/home/lsh1901704/causal_k_mean1_true_low")
write.csv(causal_k_mean05_true,"/home/lsh1901704/causal_k_mean05_true_low")