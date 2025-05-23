library(rbmi)
library(dplyr)
library(mmrm)
library(LaplacesDemon)
library(rstan)
library(purrr)
library(triangle)


#---------------------------Define parameters------------------------

causal_k_mean0 <- causal_k_mean05 <- causal_k_mean1 <- causal_k_mean0_beta  <- causal_k_mean05_beta <- causal_k_mean1_beta <- 0
causal_k_mean0_se  <- causal_k_mean05_se <- causal_k_mean1_se <- causal_k_mean0_beta_se  <- causal_k_mean05_beta_se <- causal_k_mean1_beta_se <- 0
#SD in the prior for k
sigma_k0 <- seq(.0, 3, .1)    

samples <- 10000

set.seed(1234)
#-------------------------load dataset--------------------------------
head(antidepressant_data)

antidepressant_data$THERAPY <- relevel(antidepressant_data$THERAPY, ref="PLACEBO")


expandedData <- expand_locf(
  antidepressant_data,
  PATIENT = levels(antidepressant_data$PATIENT),
  VISIT = levels(antidepressant_data$VISIT),
  vars = c("BASVAL", "THERAPY"),
  group = c("PATIENT"),
  order = c("PATIENT", "VISIT")
)


ice_info <- expandedData %>%
  arrange(PATIENT, VISIT) %>%
  filter(is.na(CHANGE)) %>%
  group_by(PATIENT) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(PATIENT, VISIT) %>%
  mutate(strategy = "MAR")

ice_info <- ice_info[-which(ice_info$PATIENT == 3618),]

# generate proportions of missingness
ice_info$r <- ice_info$VISIT
data_x <- left_join(expandedData,ice_info, by = "PATIENT") 
data_x$r <- ifelse(is.na(data_x$r),0,data_x$r)
prop <- (prop.table(table(data_x$r[data_x$THERAPY=="DRUG"])))
n_drop <- table(data_x$r[data_x$THERAPY=="DRUG"])/4

ice_info <- subset(ice_info, select = c("PATIENT", "VISIT","strategy"))

bayesDraws <- draws(
  data = expandedData,
  data_ice = ice_info,
  vars = set_vars(
    outcome = "CHANGE",
    visit = "VISIT",
    subjid = "PATIENT",
    group = "THERAPY",
    covariates = c("BASVAL*VISIT", "THERAPY*VISIT")),
  method = method_bayes(
    n_samples = samples,
    control =  control_bayes(
      warmup = 100,
      thin = 10),
    same_cov = FALSE))
  #  seed = 484381136))
#set.seed(837263)


marImputations <- impute(bayesDraws)
ancovaVars <- set_vars(
  subjid = "PATIENT",
  outcome = "CHANGE",
  visit = "VISIT",
  group = "THERAPY",
  covariates = c("BASVAL"))
marFits <- analyse(
  imputations=marImputations,
  vars = ancovaVars)
marResults <- pool(marFits)
marResults

#______________________________________
#  MMRM
#______________________________________
mmrmFit <- mmrm(CHANGE~BASVAL*VISIT+THERAPY*VISIT+us(VISIT|PATIENT),
                data=antidepressant_data,
                reml=TRUE, method = "Kenward-Roger")
summary(mmrmFit)

# mmrm treat effect
(mmrmFit$beta_est['THERAPYDRUG'] + mmrmFit$beta_est['VISIT7:THERAPYDRUG'])


#______________________________________
#  RBI imputation analysis - J2R
#______________________________________
#set.seed(837263)
ice_info_jr <- ice_info %>%
  mutate(strategy = ifelse(strategy == "MAR", "JR", strategy))

jrImps <- impute(bayesDraws,
                 references = c("DRUG" = "PLACEBO", "PLACEBO" = "PLACEBO"),
                 update_strategy = ice_info_jr)

jrMIFits <- analyse(
  imputations = jrImps,
  vars = ancovaVars
)
jrMIResults <- pool(jrMIFits)
jrMIResults


#______________________________________
#  RBI imputation analysis - CIR
#______________________________________
ice_info_cir <- ice_info %>%
  mutate(strategy = ifelse(strategy == "MAR", "CIR", strategy))

#set.seed(837263)
cirImps <- impute(bayesDraws,
                  references = c("DRUG" = "PLACEBO", "PLACEBO" = "PLACEBO"),
                  update_strategy = ice_info_cir)

cirMIFits <- analyse(
  imputations = cirImps,
  vars = ancovaVars
)
cirMIResults <- pool(cirMIFits)
cirMIResults


#_____________________________________________
#  Bayesian causal model Analysis - fixed k_0
#_____________________________________________
# Obtain bayesian samples
betas <- as.data.frame(do.call(rbind, do.call(rbind, bayesDraws[["samples"]])[,3]))

# Draws for proportions - we use a Dirichlet prior Dirichlet(12,1,1,2) 
prop_sampled <- rdirichlet(samples, c(12,1,1,2) + n_drop)


#J2R
betas$j2r <- (betas$V2 + betas$V12)*prop_sampled [,1]

#CIR
betas$cir <- (betas$V2 + betas$V12)*prop_sampled [,1] + (betas$V2)*prop_sampled [,2] + (betas$V2 + betas$V10)*prop_sampled [,3] + (betas$V2 + betas$V11)*prop_sampled[,4]
#MAR
betas$mar <- (betas$V2 + betas$V12)



#_____________________________________________________
#   Bayesian causal model Analysis - A prior on k_0
#_____________________________________________________


for (i in 1:length(sigma_k0)) {
  #J2R
  K <- rnorm(samples,0,1)
  k <- K*sigma_k0[i]
  betas$causal_k_mean0<- (betas$V2 + betas$V12)*prop_sampled [,1] + k*((betas$V2)*prop_sampled [,2]) +  k*((betas$V2 + betas$V10)*prop_sampled [,3]) + k*((betas$V2 + betas$V11)*prop_sampled[,4])
  
  causal_k_mean0[i] <- mean(betas$causal_k_mean0)
  causal_k_mean0_se[i] <- sd(betas$causal_k_mean0)
  
  #CIR
  k_2 <- 1 + K*sigma_k0[i]
  betas$causal_k_mean1 <-  ((betas$V2 + betas$V12)*prop_sampled [,1]) + k_2*((betas$V2)*prop_sampled [,2]) + k_2*((betas$V2 + betas$V10)*prop_sampled [,3]) + k_2*((betas$V2 + betas$V11)*prop_sampled[,4])
  
  causal_k_mean1[i] <- mean(betas$causal_k_mean1)
  causal_k_mean1_se[i] <- sd(betas$causal_k_mean1)
  
  
  k_3 <- 0.5 + K*sigma_k0[i]
  betas$causal_k_mean05 <-  ((betas$V2 + betas$V12)*prop_sampled [,1]) + k_3*((betas$V2)*prop_sampled [,2]) + k_3*((betas$V2 + betas$V10)*prop_sampled [,3]) + k_3*((betas$V2 + betas$V11)*prop_sampled[,4])
  causal_k_mean05[i] <- mean(betas$causal_k_mean05)
  causal_k_mean05_se[i] <- sd(betas$causal_k_mean05)
  
}


#----------------Triangular  distribution for k_0-------------------------------

# mean = 0
k <- rtriangle(samples, a=0,b=0.25,c=0.0)  

betas$causal_k_mean0_triangular <-  ((betas$V2 + betas$V12)*prop_sampled [,1]) + k*((betas$V2)*prop_sampled [,2]) + k*((betas$V2 + betas$V10)*prop_sampled [,3]) + k*((betas$V2 + betas$V11)*prop_sampled[,4])




causal_k_mean0_triangular = cbind(apply(as.data.frame(betas$causal_k_mean0_triangular),2,mean),apply(as.data.frame(betas$causal_k_mean0_triangular),2,sd),
                       apply(as.data.frame(betas$causal_k_mean0_triangular),2,quantile,c(0.025)),apply(as.data.frame(betas$causal_k_mean0_triangular),2,quantile,c(0.975)))

# mean = 0.5

k <- rtriangle(samples, a=0,b=1,c=0.5)  

betas$causal_k_mean05_triangular <-  ((betas$V2 + betas$V12)*prop_sampled [,1]) + k*((betas$V2)*prop_sampled [,2]) + k*((betas$V2 + betas$V10)*prop_sampled [,3]) + k*((betas$V2 + betas$V11)*prop_sampled[,4])



causal_k_mean05_triangular = cbind(apply(as.data.frame(betas$causal_k_mean05_triangular),2,mean),apply(as.data.frame(betas$causal_k_mean05_triangular),2,sd),
                    apply(as.data.frame(betas$causal_k_mean05_triangular),2,quantile,c(0.025)),apply(as.data.frame(betas$causal_k_mean05_triangular),2,quantile,c(0.975)))

#----------------Trancted normal  distribution for k_0-------------------------------

# mean = 0
k <- abs(rnorm(samples,mean=0.0,sd=0.4))  

betas$causal_k_mean0_tranctednorm <-  ((betas$V2 + betas$V12)*prop_sampled [,1]) + k*((betas$V2)*prop_sampled [,2]) + k*((betas$V2 + betas$V10)*prop_sampled [,3]) + k*((betas$V2 + betas$V11)*prop_sampled[,4])


causal_k_mean0_tranctednorm = cbind(apply(as.data.frame(betas$causal_k_mean0_tranctednorm),2,mean),apply(as.data.frame(betas$causal_k_mean0_tranctednorm),2,sd),
                    apply(as.data.frame(betas$causal_k_mean0_tranctednorm),2,quantile,c(0.025)),apply(as.data.frame(betas$causal_k_mean0_tranctednorm),2,quantile,c(0.975)))

# mean = 0.5

k <- abs(rnorm(samples,mean=0.5,sd=0.4))  

betas$causal_k_mean05_tranctednorm <-  ((betas$V2 + betas$V12)*prop_sampled [,1]) + k*((betas$V2)*prop_sampled [,2]) + k*((betas$V2 + betas$V10)*prop_sampled [,3]) + k*((betas$V2 + betas$V11)*prop_sampled[,4])



causal_k_mean05_tranctednorm = cbind(apply(as.data.frame(betas$causal_k_mean05_tranctednorm),2,mean),apply(as.data.frame(betas$causal_k_mean05_tranctednorm),2,sd),
                    apply(as.data.frame(betas$causal_k_mean05_tranctednorm),2,quantile,c(0.025)),apply(as.data.frame(betas$causal_k_mean05_tranctednorm),2,quantile,c(0.975)))

#Summarise results
res <- subset(betas, select = c(j2r,cir,mar))
sumRes = cbind(apply(res,2,mean),apply(res,2,sd),apply(res,2,quantile,c(0.025)),apply(res,2,quantile,c(0.975)))


res_norm <- as.data.frame(cbind(causal_k_mean0,causal_k_mean0_se,causal_k_mean05,causal_k_mean05_se,causal_k_mean1,causal_k_mean1_se))

res_beta <- as.data.frame(cbind(causal_k_mean0_beta,causal_k_mean0_beta_se,causal_k_mean1_beta,causal_k_mean1_beta_se))
res_beta_2 <- as.data.frame(cbind(causal_k_mean05_beta,causal_k_mean05_beta_se))

sumRes
