# This code demonstrates the models for analysis one 
# It does not include all model checking and only some graphical output

# Please note the data is randomly generated and serves only to demonstrate the code
# Any associations identified are likely spurious
# For simplicity each individual who receives a transfusion only receives one

# All code is my own, name not supplied to maintain blinding for marking

## Required packages - uncomment as necessary ##
#install.packages("tidyverse")
#install.packages("survival")
#install.packages("ggsurvfit")

# Loading essential packages
library(magrittr)

# Setting working directory to location of dataset
# setwd()

# Clearing the global environment
rm(list = ls())

#------------------------#
# Creating survival data #
#------------------------#

load("df2.Rda")

data <- df2 %>%
  tidyr::drop_na(rethnic, prd_gp)

# Converting the transfusion dates to days relative to activation
days <- data %>%
  dplyr::mutate(post = ifelse(days>(end-3), 0, 1)) %>%
  dplyr::filter(post==1) %>%
  dplyr::select(-post) %>%
  dplyr::mutate(pre = ifelse(days<0, 0, 1)) %>%
  dplyr::filter(pre==1) %>%
  dplyr::select(-pre) %>%
  dplyr::rename(trf=days) %>%
  dplyr::select(code, trf) %>%
  dplyr::full_join(data, ., by=dplyr::join_by(code)) %>%
  dplyr::select(-days)

# Adding an indicator for pre-eskd transfusion
days <- data %>%
  dplyr::mutate(pre = ifelse(days<0 & abs(days) > wait, 0, 1)) %>%
  dplyr::filter(pre==0) %>%
  dplyr::select(-pre) %>%
  dplyr::mutate(hist=1) %>%
  dplyr::select(code, hist) %>%
  dplyr::full_join(days, ., by=dplyr::join_by(code)) %>%
  dplyr::mutate(hist = ifelse(is.na(hist),0,1)) %>%
  dplyr::select(code:prev_tx, hist, trf)

# Adding an indicator for post-eskd, pre-activation transfusion

days <- data %>%
  dplyr::mutate(pre = ifelse(days<0 & abs(days) < wait, 0, 1)) %>%
  dplyr::filter(pre==0) %>%
  dplyr::select(-pre) %>%
  dplyr::mutate(dial=1) %>%
  dplyr::select(code, dial) %>%
  dplyr::full_join(days, ., by=dplyr::join_by(code)) %>%
  dplyr::mutate(dial = ifelse(is.na(dial),0,1)) %>%
  dplyr::select(code:prev_tx, hist, dial, trf)

days <- days %>%
  dplyr::select(code:end, hist, dial, trf) %>%
  dplyr::rename(day=trf) %>%
  dplyr::mutate(trf = ifelse(is.na(day), 0, 1))

#------------------------#
# Time-updated variables #
#------------------------#

# Creating the age variable
ages <- data %>%
  dplyr::select(code:end, rage) %>%
  dplyr::mutate(age_1 = floor(rage-end/365.25)+1) %>%
  dplyr::mutate(day_1 = 0)

for(i in 1:15) {
  x <- 365.25 * i
  ages[[paste0("age_",i+1)]] <- ifelse(ages$end>x, ages$rage-(i-1), NA)
  ages[[paste0("day_",i+1)]] <- ifelse(ages$end>x, ages$end-x, NA)
  remove(x)
}

ages <- ages %>%
  tidyr::pivot_longer(cols=c(age_1:day_16), 
               names_to =c(".value","time"), 
               names_sep ="_",
               values_drop_na = TRUE
  ) %>%
  dplyr::select(code, start, end, age, day)

# Converting to survival data
surv <- data %>%
  dplyr::left_join(days%>%dplyr::select(code, hist, dial), by = dplyr::join_by(code)) %>%
  dplyr::select(code:wait, rsex:prev_tx, hist, dial) %>%
  dplyr::mutate(status=1) %>% 
  survival::tmerge(., ., id=code, transplant=event(end,status)) %>%
  survival::tmerge(., ages, id=code, rage=tdc(day,age)) %>% # Adding the age variable
  survival::tmerge(., days, id=code, transfusion = tdc(day, trf)) %>% # Adding the transfusion variable 
  dplyr::mutate(transfusion = ifelse(is.na(transfusion), 0, 1)) %>%
  dplyr::select(code, rage, rsex, rethnic, prd_gp, prev_tx, wait, hist, dial, transfusion,
                tstart, tstop, transplant)

#-------------------------------#
# Estimand 1 - Any transfusions #
#-------------------------------#

rm(list=ls()[! ls() %in% c("surv")])

# Converting to binary transfusion variable (time-updated)
surv1 <- surv %>%
  dplyr::mutate(transfusion = ifelse(hist==0  & dial==0 & transfusion==0, 0,1)) %>%
  dplyr::mutate(transfusion = factor(transfusion , levels=c(0,1), labels=c("N","Y")) ) %>%
  dplyr::mutate(wait = wait/365.25) %>%
  dplyr::select(-c(hist,dial))

## Kaplan Meier Estimates for all transfusions
km_all <- 
  surv1 %>%
  dplyr::mutate(transfusion = relevel(transfusion, "Y")) %>%
  survival::survfit(survival::Surv(tstart/365.25,tstop/365.25,transplant) ~
            transfusion, data=.) 

km_all

# Running the cox model

model1 <- survival::coxph(survival::Surv(tstart,tstop,transplant)~transfusion 
                          + rage + rsex + rethnic + prd_gp + prev_tx + wait, 
                          data= surv1, ties="breslow")

summary(model1)

# An example Martingale Plot for age

age1 <- surv1 %>%
  dplyr::mutate(mtng = resid(model1, type="martingale")) %>%
  dplyr::select(rage:mtng) %>%
  ggplot2::ggplot(ggplot2::aes(x=rage, y=mtng)) + 
  ggplot2::geom_point(colour="grey",alpha=0.4) +
  ggplot2::geom_smooth(method = loess, span = 1, fill = "#ffcccc", color="#ff0000") +
  ggplot2::geom_hline(yintercept = 0, linetype="dotted") +
  ggplot2::xlab("Age (years)") +
  ggplot2::ylab("") +
  ggplot2::theme_classic()

age1

# Schoenfeld Residual Plots

sch.resid <- survival::cox.zph(model1, transform='identity')

sch.resid

sch.resid1 <- sch.resid$y %>%
  as.data.frame() %>%
  tibble::rownames_to_column("time") %>%
  dplyr::mutate(time = as.numeric(gsub("X", "", time))) %>%
  dplyr::rename("Transfusion" = "transfusion") %>%
  dplyr::rename("Age" = "rage") %>%
  dplyr::rename("Sex" = "rsex") %>%
  dplyr::rename("Ethnicity" = "rethnic") %>%
  dplyr::rename("Aetiology" = "prd_gp") %>%
  dplyr::rename("Previous Transplant" = "prev_tx") %>%
  dplyr::rename("Waiting Points" = "wait") 

list <- sch.resid1 %>% dplyr::select(-time) %>% names() 

sch_plot <- sch.resid1 %>%
  tidyr::pivot_longer(
    cols=dplyr::all_of(list)
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x=time, y=value)) +
  ggplot2::geom_point(colour="grey", alpha=0.2) +
  ggplot2::geom_smooth(method = loess, span = 1, fill = "#ccccff", color="#0000ff") + 
  ggplot2::geom_hline(yintercept = 0, linetype="dotted") + 
  ggplot2::ylab("") + 
  ggplot2::theme_classic() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()) +
  ggplot2::facet_wrap(~name, scales = "free") +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
        strip.background = ggplot2::element_rect(color="white"))

sch_plot

# Estimating the curves for individuals with median/modal values 

# Function to extract the modal value for categorical variables
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


fit <- ggsurvfit::survfit2(model1, newdata = data.frame(
  rage = median(surv1$rage),
  rsex = Mode(surv1$rsex), 
  rethnic = Mode(surv1$rethnic),
  prd_gp = Mode(surv1$prd_gp),
  prev_tx = Mode(surv1$prev_tx),
  wait = median(surv1$wait),
  transfusion=c("N","Y")))

fit

#-------------------------------------------#
# Estimand 2 - Post-activation transfusions #
#-------------------------------------------#

rm(list=ls()[! ls() %in% c("surv", "Mode")])

surv2 <- surv %>%
  dplyr::mutate(wait = wait/365.25) %>%
  dplyr::mutate(pre = ifelse(hist==0  & dial==0, 0,1)) %>%
  dplyr::mutate(pre = factor(pre , levels=c(0,1), labels=c("N","Y"))) %>%
  dplyr::mutate(transfusion = factor(transfusion , levels=c(0,1), labels=c("N","Y")))

# Kaplan Meier Estimates

# Transfusion prior to activation
km_pre <- 
  surv2 %>%
  dplyr::mutate(pre = relevel(pre, "Y")) %>%
  survival::survfit(survival::Surv(tstart/365.25,tstop/365.25,transplant) ~
                      pre, data=.) 

km_pre

# Transfusion after activation
km_post <- 
  surv2 %>%
  dplyr::mutate(transfusion = relevel(transfusion, "Y")) %>%
  survival::survfit(survival::Surv(tstart/365.25,tstop/365.25,transplant) ~
                      transfusion, data=.) 

km_post

# Running the cox model
model2 <- survival::coxph(survival::Surv(tstart,tstop,transplant)~transfusion + pre
                          + rage + rsex + rethnic + prd_gp + prev_tx + wait, 
                          data= surv2, ties="breslow")

summary(model2)


# Schoenfeld Residuals

sch.resid2 <- survival::cox.zph(model2, transform='identity')
sch.resid2

# Estimating the curves for individuals with median/modal values 
#NB the sample data has no individuals with transfusions in both times so this is an extrapolation

ggsurvfit::survfit2(model2, newdata = data.frame(
  rage = median(surv2$rage),
  rsex = Mode(surv2$rsex), 
  rethnic = Mode(surv2$rethnic),
  prd_gp = Mode(surv2$prd_gp),
  prev_tx = Mode(surv2$prev_tx),
  wait = median(surv2$wait),
  pre = c("N","Y", "N","Y"),
  transfusion=c("N","N","Y","Y")))
