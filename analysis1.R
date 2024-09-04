# This code demonstrates the models for analysis one 
# It does not include all model checking or any graphical output

# Please note the data is randomly generated and serves only to demonstrate the code
# Any associations identified are likely spurious

# All code is my own, name not supplied to maintain blinding for marking

## Required packages - uncomment as necessary ##
#install.packages("tidyverse")
#install.packages("mice")
#install.packages("egg")
#install.packages("cowplot")
#install.packages("twopartm")
#install.packages("betareg")
#install.packages("marginaleffects")
#install.packages("brms")

# Loading essential packages
library(magrittr)

# Setting working directory to location of dataset
# setwd()

# Clearing the global environment and loading the dataset
rm(list = ls())

load("df1.Rda")

#--------------------#
# Defining functions #
#--------------------#

# Defining the function for Rubin's rules
rubins <- function(est,sde,data,k) {
  if(length(est)==length(sde)) {
    x <- length(est)
    mu <- mean(est[1:x])
    fun <- function(y){
      (y-mu)^2
    }
    z <- sapply(est, fun)
    var_between <- 1/(x-1)*sum(z)
    func <- function(w){
      w^2
    }
    a <- sapply(sde, func)
    var_within <- mean(a)
    var_total <- ((1+(1/x))*(var_between)) + var_within
    sigma <- sqrt(var_total)
    lambda <- (var_between + (var_between/x))/(var_total)
    nu_old <- (x-1)/lambda
    nu_obs <- ((nrow(data)-k) + 1)/((nrow(data)-k) + 3)*(nrow(data)-k)*(1-lambda)
    nu_adj <- (nu_old*nu_obs)/(nu_old+nu_obs)
    t <- qt(0.975,nu_obs)
    lci <- mu-(t*sigma)
    uci <- mu+(t*sigma)
    p1 <- pt((mu/sigma),nu_adj, lower.tail = FALSE)
    p2 <- pt((mu/sigma),nu_adj, lower.tail = TRUE)
    p <- ifelse(p1>p2, 2*p2, 2*p1)
    print(c("Estimate"=mu,"Std. Err"=sigma, "Lower" = lci, "Upper" = uci, "P-value"=p))
  } else {
    print("Length Mismatch")
  }
}

#----------------------------#
# Generating the Imputations #
#----------------------------#

# Preventing overwriting of original data
data <- df1 

# Removing the code variable from the prediction matrix
matrix <- mice::make.predictorMatrix(data)
matrix[,c("code")] <- 0 

set.seed(123456)
imp <- mice::mice(data, m=30, maxit=10, predictorMatrix = matrix,
             printFlag = FALSE)

# Extracting the imputations for ethnicity
imp_e <- imp$imp$rethnic %>%
  setNames(paste0("e",names(.))) %>%
  merge(data,.,by=0,all.x=T) %>%
  dplyr::select(-Row.names)

# Extracting the imputations for aetiology
imp_p <- imp$imp$prd_gp %>%
  setNames(paste0("p",names(.))) %>%
  merge(data,.,by=0,all.x=T) %>%
  dplyr::select(-Row.names)

#-------------------------------------------#
# Calculating propensity scores and weights #
#-------------------------------------------#

# Propensity Score formula
ps_formula <- transfusion ~ rage + rsex + rethnic + prd_gp + prev_tx + krt

# Calculating the propensity score and weight for each imputation
for (i in 1:30) {
  imp_e$rethnic <- dplyr::coalesce(imp_e[[paste0("e",i)]],imp_e[['rethnic']]) 
  imp_p$prd_gp <- dplyr::coalesce(imp_p[[paste0("p",i)]],imp_p[['prd_gp']])
  data <- imp_e %>% dplyr::select(code, rethnic) %>% 
    dplyr::full_join(data %>% dplyr::select(-rethnic), by = dplyr::join_by(code) )
  data <- imp_p %>% dplyr::select(code, prd_gp) %>% 
    dplyr::full_join(data %>% dplyr::select(-prd_gp), by = dplyr::join_by(code))
  m_ps <-glm(ps_formula, data=data ,family="binomial")
  data[[paste0("ps",i)]] <- predict(m_ps,type="response")
  data[[paste0("w",i)]] <- ifelse(data$transfusion==1, 
                                  1/data[[paste0("ps",i)]], 
                                  1/(1-data[[paste0("ps",i)]])) # ATE weights
}

#----------------#
# Two-Part Model #
#----------------#

# Running the 2 part model
imps <- data %>%
  dplyr::select(code, transfusion, starts_with("w"), crf)

# Creating vectors for results
p1i <- p1t <- p1d <- se1i <- se1t <- se1d <- NA
e2i <- e2t <- e2d <- se2i <- se2t <- se2d <- NA
ei <- sei <- et <- set <- ed <- sed <- NA

# Calculating each imputation estimate for the first and second part 
for (i in 1:30) {
  modela <- twopartm::tpm(formula_part1 = crf ~ transfusion, link_part1 = "logit", 
               family_part2 = gaussian(link = "identity"),
               data = imps, weights=imps[[paste0("w",i)]])
  moda <- twopartm::summary(modela)
  modelb <- twopartm::tpm(formula_part1 = crf ~ 
                          factor(transfusion, levels=c("Y","N")), link_part1 = "logit", 
                          family_part2 = gaussian(link = "identity"),
                          data = imps, weights=imps[[paste0("w",i)]])
  modb <- twopartm::summary(modelb)
  p1i[i] <- moda$Firstpart.model$coefficients["(Intercept)","Estimate"]
  p1t[i] <- modb$Firstpart.model$coefficients["(Intercept)","Estimate"]
  p1d[i] <- moda$Firstpart.model$coefficients["transfusionY","Estimate"]
  se1i[i] <- moda$Firstpart.model$coefficients["(Intercept)","Std. Error"]
  se1t[i] <- modb$Firstpart.model$coefficients["(Intercept)","Std. Error"]
  se1d[i] <- moda$Firstpart.model$coefficients["transfusionY","Std. Error"]
  e2i[i] <- moda$Secondpart.model$coefficients["(Intercept)","Estimate"]
  e2t[i] <- modb$Secondpart.model$coefficients["(Intercept)","Estimate"]
  e2d[i] <- moda$Secondpart.model$coefficients["transfusionY","Estimate"]
  se2i[i] <- moda$Secondpart.model$coefficients["(Intercept)","Std. Error"]
  se2t[i] <- modb$Secondpart.model$coefficients["(Intercept)","Std. Error"]
  se2d[i] <- moda$Secondpart.model$coefficients["transfusionY","Std. Error"]
  predn <- twopartm::predict(modela, newdata=data.frame(transfusion="N"), se.fit=TRUE)
  predy <- twopartm::predict(modela, newdata=data.frame(transfusion="Y"), se.fit=TRUE)
  ei[i] <- predn$fit
  et[i] <- predy$fit
  sei[i] <- predn$se.fit
  set[i] <- predy$se.fit
  ed[i] <- predy$fit - predn$fit
  sed[i] <- sqrt((predn$se.fit^2) + (predy$se.fit^2))
  remove(moda, modb)
}

# Applying Rubin's rules

# First part
p1i_r <- rubins(p1i,se1i,imps,1) # Log Odds(Y=0|X=0)
p1t_r <- rubins(p1t,se1t,imps,1) # Log Odds(Y=0|X=1)
p1d_r <- rubins(p1d,se1d,imps,1) # Log Odds ratio

# Second part
e2i_r <- rubins(e2i,se2i,imps,1) # E(Y|X=0, Y>0) 
e2t_r <- rubins(e2t,se2t,imps,1) # E(Y|X=1, Y>0)
e2d_r <- rubins(e2d,se2d,imps,1) # Difference

# Overall
ei_r <- rubins(ei, sei, imps, 1) # E(Y|X=0)
et_r <- rubins(et, set, imps,1) # E(Y|X=1)
ed_r <- rubins(ed, sed, imps, 1) # Difference

#-----------------#
# Beta Regression #
#-----------------#

rm(list=ls()[! ls() %in% c("df1", "imps", "rubins")])

# Scaling the data
imps2 <- imps %>%
  dplyr::mutate(crf = (crf*(nrow(imps)-1)+0.5)/nrow(imps)) 

# Creating vectors for results
e <- se <- e2 <- se2 <- NA
pi <- sei <- pt <- set <- AIC <-  NA
phii <- phisei <- phit <- phiset <- NA

# Calculating each imputation estimate for the beta regression
for (i in 1:30) {
  modela <- betareg::betareg(crf ~ transfusion | transfusion, 
                    data=imps2, weights=imps2[[paste0("w",i)]])
  modelb <- betareg::betareg(crf ~ factor(transfusion, levels=c("Y","N")) | 
                      factor(transfusion, levels=c("Y","N")), 
                    data=imps2, weights=imps2[[paste0("w",i)]])
  diff <- marginaleffects::avg_slopes(modela)
  mod1 <- summary(modela)
  mod2 <- summary(modelb)
  e[i] <- mod1$coefficients$mean["transfusionY","Estimate"]
  se[i] <- mod1$coefficients$mean["transfusionY","Std. Error"]
  e2[i] <- diff$estimate
  se2[i] <- diff$std.error
  pi[i] <- mod1$coefficients$mean["(Intercept)","Estimate"]
  sei[i] <- mod1$coefficients$mean["(Intercept)","Std. Error"]
  pt[i] <- mod2$coefficients$mean["(Intercept)","Estimate"]
  set[i] <- mod2$coefficients$mean["(Intercept)","Std. Error"]
  phii[i] <-mod1$coefficients$precision["(Intercept)","Estimate"]
  phisei[i] <-mod1$coefficients$precision["(Intercept)","Std. Error"]
  phit[i] <-mod1$coefficients$precision["transfusionY","Estimate"]
  phiset[i] <-mod1$coefficients$precision["transfusionY","Std. Error"]
  AIC[i] <- AIC(modela)
  remove(modela, modelb, mod1, mod2)
}


beta <- rubins(e, se, imps, 1) # Log Odds ratio

trfn <- rubins(pi, sei, imps, 1) # logit[E(Y|X=0)]
trfy <- rubins(pt, set, imps, 1) # logit[E(Y|X=1)]
betadiff <- rubins(e2, se2, imps,1) # Difference 

phin <- rubins(phii, phisei, imps, 1) # Log (Precision | X=0)
phiy <- rubins(phit, phiset, imps, 1) # Log (Precision | X=1)

# Beta parameters 

# Transfusion = N
a_n <- plogis(trfn[1]) * exp(phin[1])
b_n <- exp(phin[1])- a_n

# Transfusion = Y
a_y <- plogis(trfn[1]+beta[1]) * exp(phin[1]+phiy[1])
b_y <- exp(phin[1]+phiy[1]) - a_y

#-----------------#
# ZOIB Regression #
#-----------------#

rm(list=ls()[! ls() %in% c("df1", "rubins", "imps")])

# Running the model across the imputations and saving the results.

# Parameters
muN <- muN_se <- muY <- muY_se <- muC <- muC_se <- NA
zoiN <- zoiN_se <- zoiY <- zoiY_se <- zoiC <- zoiC_se <- NA
coiN <- coiN_se <- coiY <- coiY_se <- coiC <- coiC_se <- NA
phiN <- phiN_se <- phiY <- phiY_se <- NA

# Transformed parameters
p0N <- p0N_se <- p0Y <- p0Y_se <- p0C <- p0C_se <- NA
p1N <- p1N_se <- p1Y <- p1Y_se <- p1C <- p1C_se <- NA
totN <- totN_se <- totY <- totY_se <- totC <- totC_se <- NA

# Running the full model once
for (i in 1:1){
  model<- brms::brm(
    formula = brms::bf(
      reformulate("transfusion", response = sprintf("crf | weights(w%d)", i)),
      phi ~ transfusion,
      zoi ~ transfusion,
      coi ~ transfusion,
      family = brms::zero_one_inflated_beta()
    ),
    data = imps,
    silent = 2,
    refresh = 0,
    seed = 123456
  )
  
  draws <- brms::as_draws_df(model, pars = "b_") %>%
    as.data.frame() %>%
    dplyr::mutate(muN = plogis(b_Intercept)) %>%
    dplyr::mutate(muY = plogis(b_Intercept+b_transfusionY)) %>%
    dplyr::mutate(phiN = exp(b_phi_Intercept)) %>%
    dplyr::mutate(phiY = exp(b_phi_Intercept + b_phi_transfusionY)) %>%
    dplyr::mutate(zoiN = plogis(b_zoi_Intercept)) %>%
    dplyr::mutate(zoiY = plogis(b_zoi_Intercept+b_zoi_transfusionY)) %>%
    dplyr::mutate(zoiC = exp(b_zoi_transfusionY)) %>%
    dplyr::mutate(coiN = plogis(b_coi_Intercept)) %>%
    dplyr::mutate(coiY = plogis(b_coi_Intercept+b_coi_transfusionY)) %>%
    dplyr::mutate(coiC = exp(b_coi_transfusionY)) %>%
    dplyr::mutate(p0N = plogis(b_zoi_Intercept)*(1 - plogis(b_coi_Intercept))) %>%
    dplyr::mutate(p0Y = plogis(b_zoi_Intercept + b_zoi_transfusionY)*
             (1 - plogis(b_coi_Intercept + b_coi_transfusionY))) %>%
    dplyr::mutate(p1N = plogis(b_zoi_Intercept)*(plogis(b_coi_Intercept))) %>%
    dplyr::mutate(p1Y = plogis(b_zoi_Intercept + b_zoi_transfusionY)
           *(plogis(b_coi_Intercept + b_coi_transfusionY))) %>%
    dplyr::mutate(totN = plogis(b_zoi_Intercept)*plogis(b_coi_Intercept)
           + plogis(b_Intercept)*(1 - plogis(b_zoi_Intercept))) %>%
    dplyr::mutate(totY = plogis(b_zoi_Intercept + b_zoi_transfusionY)*
             plogis(b_coi_Intercept + b_coi_transfusionY) 
           + plogis(b_Intercept + b_transfusionY)*
             (1 - plogis(b_zoi_Intercept + b_zoi_transfusionY))) %>%
    dplyr::select(muN:totY) %>%
    brms::posterior_summary() %>% 
    as.data.frame() 
  
  muN[i] <- draws["muN","Estimate"]
  muN_se[i] <- draws["muN","Est.Error"]
  muY[i] <- draws["muY","Estimate"]
  muY_se[i] <- draws["muY","Est.Error"]
  phiN[i] <- draws["phiN","Estimate"]
  phiN_se[i] <- draws["phiN","Est.Error"]
  phiY[i] <- draws["phiY","Estimate"]
  phiY_se[i] <- draws["phiY","Est.Error"]
  zoiN[i] <- draws["zoiN","Estimate"]
  zoiN_se[i] <- draws["zoiN","Est.Error"]
  zoiY[i] <- draws["zoiY","Estimate"]
  zoiY_se[i] <- draws["zoiY","Est.Error"]
  zoiC[i] <- draws["zoiC","Estimate"]
  zoiC_se[i] <- draws["zoiC","Est.Error"]
  coiN[i] <- draws["coiN","Estimate"]
  coiN_se[i] <- draws["coiN","Est.Error"]
  coiY[i] <- draws["coiY","Estimate"]
  coiY_se[i] <- draws["coiY","Est.Error"]
  coiC[i] <- draws["coiC","Estimate"]
  coiC_se[i] <- draws["coiC","Est.Error"]
  p0N[i] <- draws["p0N","Estimate"]
  p0N_se[i] <- draws["p0N","Est.Error"]
  p0Y[i] <- draws["p0Y","Estimate"]
  p0Y_se[i] <- draws["p0Y","Est.Error"]
  p1N[i] <- draws["p1N","Estimate"]
  p1N_se[i] <- draws["p1N","Est.Error"]
  p1Y[i] <- draws["p1Y","Estimate"]
  p1Y_se[i] <- draws["p1Y","Est.Error"]
  totN[i] <- draws["totN","Estimate"]
  totN_se[i] <- draws["totN","Est.Error"]
  totY[i] <- draws["totY","Estimate"]
  totY_se[i] <- draws["totY","Est.Error"]
  
  hyp <- brms::hypothesis(model, 
                    c("mu" = "plogis(Intercept + transfusionY) = plogis(Intercept)",
                      "p0" = "plogis(zoi_Intercept + zoi_transfusionY)*(1 - plogis(coi_Intercept + coi_transfusionY)) = plogis(zoi_Intercept)*(1 - plogis(coi_Intercept))",
                      "p1" = "plogis(zoi_Intercept + zoi_transfusionY)*(plogis(coi_Intercept + coi_transfusionY)) = plogis(zoi_Intercept)*(plogis(coi_Intercept))",
                      "tot" = "plogis(zoi_Intercept + zoi_transfusionY)*plogis(coi_Intercept + coi_transfusionY) + plogis(Intercept + transfusionY)*(1 - plogis(zoi_Intercept + zoi_transfusionY)) = plogis(zoi_Intercept)*plogis(coi_Intercept) + plogis(Intercept)*(1 - plogis(zoi_Intercept))"
                    )
  )
  muC[i] <- hyp$hypothesis$Estimate[hyp$hypothesis$Hypothesis=="mu"]
  muC_se[i] <-hyp$hypothesis$Est.Error[hyp$hypothesis$Hypothesis=="mu"]
  p0C[i] <- hyp$hypothesis$Estimate[hyp$hypothesis$Hypothesis=="p0"]
  p0C_se[i] <-hyp$hypothesis$Est.Error[hyp$hypothesis$Hypothesis=="p0"]
  p1C[i] <- hyp$hypothesis$Estimate[hyp$hypothesis$Hypothesis=="p1"]
  p1C_se[i] <-hyp$hypothesis$Est.Error[hyp$hypothesis$Hypothesis=="p1"]
  totC[i] <- hyp$hypothesis$Estimate[hyp$hypothesis$Hypothesis=="tot"]
  totC_se[i] <-hyp$hypothesis$Est.Error[hyp$hypothesis$Hypothesis=="tot"]
  remove(draws,hyp)
}

# Updating for the other iterations
for (i in 2:30) {
  newform <- brms::bf(
    reformulate("transfusion", response = sprintf("crf | weights(w%d)", i)),
    phi ~ transfusion,
    zoi ~ transfusion,
    coi ~ transfusion,
    family = brms::zero_one_inflated_beta()
  )
  
  m_new <- update(model, newform, newdata = imps, 
                  silent = 2, refresh = 0, seed = 123456)
  
  draws <- brms::as_draws_df(m_new, pars = "b_") %>%
    as.data.frame() %>%
    dplyr::mutate(muN = plogis(b_Intercept)) %>%
    dplyr::mutate(muY = plogis(b_Intercept+b_transfusionY)) %>%
    dplyr::mutate(phiN = exp(b_phi_Intercept)) %>%
    dplyr::mutate(phiY = exp(b_phi_Intercept + b_phi_transfusionY)) %>%
    dplyr::mutate(zoiN = plogis(b_zoi_Intercept)) %>%
    dplyr::mutate(zoiY = plogis(b_zoi_Intercept+b_zoi_transfusionY)) %>%
    dplyr::mutate(zoiC = exp(b_zoi_transfusionY)) %>%
    dplyr::mutate(coiN = plogis(b_coi_Intercept)) %>%
    dplyr::mutate(coiY = plogis(b_coi_Intercept+b_coi_transfusionY)) %>%
    dplyr::mutate(coiC = exp(b_coi_transfusionY)) %>%
    dplyr::mutate(p0N = plogis(b_zoi_Intercept)*(1 - plogis(b_coi_Intercept))) %>%
    dplyr::mutate(p0Y = plogis(b_zoi_Intercept + b_zoi_transfusionY)*
             (1 - plogis(b_coi_Intercept + b_coi_transfusionY))) %>%
    dplyr::mutate(p1N = plogis(b_zoi_Intercept)*(plogis(b_coi_Intercept))) %>%
    dplyr::mutate(p1Y = plogis(b_zoi_Intercept + b_zoi_transfusionY)
           *(plogis(b_coi_Intercept + b_coi_transfusionY))) %>%
    dplyr::mutate(totN = plogis(b_zoi_Intercept)*plogis(b_coi_Intercept)
           + plogis(b_Intercept)*(1 - plogis(b_zoi_Intercept))) %>%
    dplyr::mutate(totY = plogis(b_zoi_Intercept + b_zoi_transfusionY)*
             plogis(b_coi_Intercept + b_coi_transfusionY) 
           + plogis(b_Intercept + b_transfusionY)*
             (1 - plogis(b_zoi_Intercept + b_zoi_transfusionY))) %>%
    dplyr::select(muN:totY) %>%
    brms::posterior_summary() %>% 
    as.data.frame() 
  
  muN[i] <- draws["muN","Estimate"]
  muN_se[i] <- draws["muN","Est.Error"]
  muY[i] <- draws["muY","Estimate"]
  muY_se[i] <- draws["muY","Est.Error"]
  phiN[i] <- draws["phiN","Estimate"]
  phiN_se[i] <- draws["phiN","Est.Error"]
  phiY[i] <- draws["phiY","Estimate"]
  phiY_se[i] <- draws["phiY","Est.Error"]
  zoiN[i] <- draws["zoiN","Estimate"]
  zoiN_se[i] <- draws["zoiN","Est.Error"]
  zoiY[i] <- draws["zoiY","Estimate"]
  zoiY_se[i] <- draws["zoiY","Est.Error"]
  zoiC[i] <- draws["zoiC","Estimate"]
  zoiC_se[i] <- draws["zoiC","Est.Error"]
  coiN[i] <- draws["coiN","Estimate"]
  coiN_se[i] <- draws["coiN","Est.Error"]
  coiY[i] <- draws["coiY","Estimate"]
  coiY_se[i] <- draws["coiY","Est.Error"]
  coiC[i] <- draws["coiC","Estimate"]
  coiC_se[i] <- draws["coiC","Est.Error"]
  p0N[i] <- draws["p0N","Estimate"]
  p0N_se[i] <- draws["p0N","Est.Error"]
  p0Y[i] <- draws["p0Y","Estimate"]
  p0Y_se[i] <- draws["p0Y","Est.Error"]
  p1N[i] <- draws["p1N","Estimate"]
  p1N_se[i] <- draws["p1N","Est.Error"]
  p1Y[i] <- draws["p1Y","Estimate"]
  p1Y_se[i] <- draws["p1Y","Est.Error"]
  totN[i] <- draws["totN","Estimate"]
  totN_se[i] <- draws["totN","Est.Error"]
  totY[i] <- draws["totY","Estimate"]
  totY_se[i] <- draws["totY","Est.Error"]
  
  hyp <- brms::hypothesis(m_new, 
                    c("mu" = "plogis(Intercept + transfusionY) = plogis(Intercept)",
                      "p0" = "plogis(zoi_Intercept + zoi_transfusionY)*(1 - plogis(coi_Intercept + coi_transfusionY)) = plogis(zoi_Intercept)*(1 - plogis(coi_Intercept))",
                      "p1" = "plogis(zoi_Intercept + zoi_transfusionY)*(plogis(coi_Intercept + coi_transfusionY)) = plogis(zoi_Intercept)*(plogis(coi_Intercept))",
                      "tot" = "plogis(zoi_Intercept + zoi_transfusionY)*plogis(coi_Intercept + coi_transfusionY) + plogis(Intercept + transfusionY)*(1 - plogis(zoi_Intercept + zoi_transfusionY)) = plogis(zoi_Intercept)*plogis(coi_Intercept) + plogis(Intercept)*(1 - plogis(zoi_Intercept))"
                    )
  )
  muC[i] <- hyp$hypothesis$Estimate[hyp$hypothesis$Hypothesis=="mu"]
  muC_se[i] <-hyp$hypothesis$Est.Error[hyp$hypothesis$Hypothesis=="mu"]
  p0C[i] <- hyp$hypothesis$Estimate[hyp$hypothesis$Hypothesis=="p0"]
  p0C_se[i] <-hyp$hypothesis$Est.Error[hyp$hypothesis$Hypothesis=="p0"]
  p1C[i] <- hyp$hypothesis$Estimate[hyp$hypothesis$Hypothesis=="p1"]
  p1C_se[i] <-hyp$hypothesis$Est.Error[hyp$hypothesis$Hypothesis=="p1"]
  totC[i] <- hyp$hypothesis$Estimate[hyp$hypothesis$Hypothesis=="tot"]
  totC_se[i] <-hyp$hypothesis$Est.Error[hyp$hypothesis$Hypothesis=="tot"]
  remove(draws,hyp)
  print(i)
}

# Alpha (P(Y=0/1))
zoiN_est <- rubins(zoiN, zoiN_se, imps, 1) # X = 0
zoiY_est <- rubins(zoiY, zoiY_se, imps, 1) # X = 1
zoiC_est <- rubins(zoiC, zoiC_se, imps, 1) # Odds ratio
palpha <- pnorm(1, mean=zoiC_est[1], sd = zoiC_est[2])
palpha <- ifelse(palpha>0.5, 2*(1-palpha), 2*palpha)
palpha <- 1-palpha # Prob(alpha<>1)

# Gamma (P(Y=1|Y=0/1))
coiN_est <- rubins(coiN, coiN_se, imps, 1) # X = 0
coiY_est <- rubins(coiY, coiN_se, imps, 1) # X = 1
coiC_est <- rubins(coiC, coiC_se, imps, 1) # Odds ratio
pgamma <- pnorm(1, mean=coiC_est[1], sd = coiC_est[2])
pgamma <- ifelse(pgamma>0.5, 2*(1-pgamma), 2*pgamma) 
pgamma <- 1-pgamma # Prob(gamma<>1)

# Mu (E(Y|Y<>0/1))
muN_est <- rubins(muN, muN_se, imps, 1) # X = 0
muY_est <- rubins(muY, muY_se, imps, 1) # X = 1
mu_diff <- rubins(muC, muC_se, imps, 1) # Difference
pmu <- pnorm(0, mean=mu_diff[1], sd = mu_diff[2]) 
pmu <- ifelse(pmu>0.5, 2*(1-pmu), 2*pmu) 
pmu <- 1-pmu # Prob(mu<>0)

# P(Y=0)
p0N_est <- rubins(p0N, p0N_se, imps, 1) # X = 0
p0Y_est <- rubins(p0Y, p0Y_se, imps, 1) # X = 1
p0_diff <- rubins(p0C, p0C_se, imps, 1) # Difference
p0 <- pnorm(0, mean=p0_diff[1], sd = p0_diff[2])
p0 <- ifelse(p0>0.5, 2*(1-p0), 2*p0) 
p0 <- 1-p0 # Prob(p0<>0)

# P(Y=1)
p1N_est <- rubins(p1N, p1N_se, imps, 1) # X = 0
p1Y_est <- rubins(p1Y, p1Y_se, imps, 1) # X = 1
p1_diff <- rubins(p1C, p1C_se, imps, 1) # Difference
p1 <- pnorm(0, mean=p1_diff[1], sd = p1_diff[2])
p1 <- ifelse(p1>0.5, 2*(1-p1), 2*p1) 
p1 <- 1-p1 # Prob(p1<>0)

# E(Y)
totN_est <- rubins(totN, totN_se, imps, 1) # X = 0
totY_est <- rubins(totY, totY_se, imps, 1) # X = 1
tot_diff <- rubins(totC, totC_se, imps, 1) # Difference
ptot <- pnorm(0, mean=tot_diff[1], sd = tot_diff[2])
ptot <- ifelse(ptot>0.5, 2*(1-ptot), 2*ptot) 
ptot <- 1-ptot # Prob(p1<>0)

# Precision
phiN_est <- rubins(phiN, phiN_se, imps, 1) # X = 0
phiY_est <- rubins(phiY, phiY_se, imps, 1) # X = 1

# Beta parameters 

# Transfusion = N
a_n <- muN_est[1] * phiN_est[1]
b_n <- phiN_est[1] - a_n

# Transfusion = Y
a_y <- muY_est[1] * phiY_est[1]
b_y <- phiY_est[1] - a_y

#---------------#
# Ordinal Model #
#---------------#

rm(list=ls()[! ls() %in% c("df1", "rubins")])

# Creating the ordinal categories

data <- df1 %>%
  dplyr::mutate(crf = ifelse(crf==0,0,
                      ifelse(crf<0.85,1,2))) %>%
  dplyr::mutate(crf = as.factor(crf))

# Rerunning the imputation model for the new variable 

# Removing the code variable from the prediction matrix
matrix <- mice::make.predictorMatrix(data)
matrix[,c("code")] <- 0 

set.seed(123456)
imp <- mice::mice(data, m=30, maxit=10, predictorMatrix = matrix,
                  printFlag = FALSE)

# Extracting the imputations for ethnicity
imp_e <- imp$imp$rethnic %>%
  setNames(paste0("e",names(.))) %>%
  merge(data,.,by=0,all.x=T) %>%
  dplyr::select(-Row.names)

# Extracting the imputations for aetiology
imp_p <- imp$imp$prd_gp %>%
  setNames(paste0("p",names(.))) %>%
  merge(data,.,by=0,all.x=T) %>%
  dplyr::select(-Row.names)

# Propensity Score formula
ps_formula <- transfusion ~ rage + rsex + rethnic + prd_gp + prev_tx + krt

# Calculating the propensity score and weight for each imputation
for (i in 1:30) {
  imp_e$rethnic <- dplyr::coalesce(imp_e[[paste0("e",i)]],imp_e[['rethnic']]) 
  imp_p$prd_gp <- dplyr::coalesce(imp_p[[paste0("p",i)]],imp_p[['prd_gp']])
  data <- imp_e %>% dplyr::select(code, rethnic) %>% 
    dplyr::full_join(data %>% dplyr::select(-rethnic), by = dplyr::join_by(code) )
  data <- imp_p %>% dplyr::select(code, prd_gp) %>% 
    dplyr::full_join(data %>% dplyr::select(-prd_gp), by = dplyr::join_by(code))
  m_ps <-glm(ps_formula, data=data ,family="binomial")
  data[[paste0("ps",i)]] <- predict(m_ps,type="response")
  data[[paste0("w",i)]] <- ifelse(data$transfusion==1, 
                                  1/data[[paste0("ps",i)]], 
                                  1/(1-data[[paste0("ps",i)]])) # ATE weights
}

# Running the ordinal model
imps <- data %>%
  dplyr::select(code, transfusion, starts_with("w"), crf)

# 3 categories
e <- se <- it01 <- se01 <- it12 <- se12 <- NA
it01_t <- se01_t <- it12_t <- se12_t <- NA
mn1e <- mn1se <- mn2e <- mn2se <- aic <- df <- NA

# Calculating each imputation estimate for ordinal and multinomial models
for (i in 1:30) {
  modpolr <- MASS::polr(crf ~ transfusion, data=imps, weights=imps[[paste0("w",i)]],
                Hess = TRUE)
  modpolrb <- MASS::polr(crf ~ factor(transfusion, levels=c("Y","N")), data=imps,
                  weights=imps[[paste0("w",i)]], Hess = TRUE)
  modmult <- nnet::multinom(crf ~ transfusion, data=imps, weights=imps[[paste0("w",i)]],
                      Hess = TRUE)
  mod <- summary(modpolr)
  mod2 <- summary(modpolrb)
  mod3 <- summary(modmult)
  aicp <- AIC(modpolr)
  aicm <- AIC(modmult)
  e[i] <- mod$coefficients["transfusionY","Value"]
  se[i] <- mod$coefficients["transfusionY","Std. Error"]
  it01[i] <-mod$coefficients["0|1","Value"]
  se01[i] <-mod$coefficients["0|1","Std. Error"]
  it12[i] <-mod$coefficients["1|2","Value"]
  se12[i] <-mod$coefficients["1|2","Std. Error"]
  it01_t[i] <-mod2$coefficients["0|1","Value"]
  se01_t[i] <-mod2$coefficients["0|1","Std. Error"]
  it12_t[i] <-mod2$coefficients["1|2","Value"]
  se12_t[i] <-mod2$coefficients["1|2","Std. Error"]
  mn1e[i] <- mod3$coefficients[1,"transfusionY"]
  mn1se[i] <- mod3$standard.errors[1,"transfusionY"]
  mn2e[i] <- mod3$coefficients[2,"transfusionY"]
  mn2se[i] <- mod3$standard.errors[2,"transfusionY"]
  aic[i] <- aicp-aicm
  df[i] <- modpolr$df.residual
  remove(modpolr, modpolrb, modmult, mod, mod2, mod3, aicp, aicm)
}

# Simultaneous log odds ratio
ord <- rubins(e, se, imps, 1)

int01 <- rubins(it01, se01, imps, 1) # Log Odds Cat<2 for X=0
int12 <- rubins(it12, se12, imps, 1) # Log Odds Cat<3 for X=0

int01_t <- rubins(it01_t, se01_t, imps, 1) # Log Odds Cat<2 for X=1
int12_t <- rubins(it12_t, se12_t, imps, 1) # Log Odds Cat<3 for X=0

summary(aic) # AIC Ordinal - Multinomial
