library(survey)
library(glmnet)
library(tidyverse)
library(xtable)
library(future)
library(vcd)
library(doFuture)
registerDoFuture()
plan(multisession)
source("codes/adalasso-chen-phd.R")

totals <- readRDS("data/gus-woj-sek-boot-totals.rds") %>%
  filter(substr(kod2,1,1) != 6) %>%
  mutate(kod2 = ifelse(kod2 == 95, 96, kod2),
         kod2 = ifelse(kod2 == 92, 91, kod2),
         kod2 = ifelse(kod2 == 53, 54, kod2),
         #kod1 = as.character(kod1),
         kod2 = as.character(kod2),
         rok = as.character(rok)) 

head(totals)

final_data <- readRDS("data/bkl-final.rds") %>%
  filter(zawod1 != 6, rok!=2012, nace %in% unique(totals$sekcja)) %>%
  mutate(zawod2 = ifelse(zawod2 == 95, 96, zawod2),
         zawod2 = ifelse(zawod2 == 92, 91, zawod2),
         zawod2 = ifelse(zawod2 == 99, 96, zawod2),
         zawod2 = ifelse(zawod2 == 53, 54, zawod2),
         zeros = str_extract(zawod6, "0+$"),
         zeros = str_count(zeros, "0"),
         zeros = ifelse(is.na(zeros), 6, 6-zeros),
         rok = as.character(rok)) %>%
  rename(kod1 = zawod1,
         kod2 = zawod2,
         kod6 = zawod6,
         sekcja = nace) %>%
  filter(zeros > 1) 

tic <- Sys.time()
komps <- c("komp_techniczne", "komp_matematyczne", "komp_kulturalne", "komp_komputerowe", 
           "komp_kognitywne", "komp_kierownicze", "komp_interpersonalne", 
           "komp_indywidualne", "komp_fizyczne", "komp_dyspozycyjne", "komp_biurowe")

res_alasso_kod2 <- foreach(i = 1:500, .export = c("wynik"), .verbose = TRUE) %dopar% {
  
  set.seed(i)
  
  totals_boot <- totals %>% filter(b == i)
  
  final_data_kod2 <- final_data %>%  group_by(rok) %>%  
    sample_frac(1, replace = T) %>% ungroup() %>%
    add_count(rok, name = "m") %>%
    left_join( totals_boot %>% count(rok, wt = hat_wolne, name = "total")) %>%
    mutate(weight = total / m ) %>%
    as.data.frame()
  
  X <- Matrix::fac2sparse(final_data_kod2$kod2) %>% t() 
  
  wynik_est <- list()
  wynik_wt <- list()
  wynik_model <- list()
  for (k in 1:length(komps)) {
    
    m1 <- mycv.glmnet(x = X, y = as.matrix(final_data_kod2[,komps[k]]), intercept = T, alpha = 0)
    m1 <- mycv.glmnet(x = X, y = as.matrix(final_data_kod2[,komps[k]]), intercept = T, penalty.factor = 1/abs(m1$coef[-1]))
    X_tot <- Matrix::fac2sparse(totals_boot$kod2) %>% t() 
    lin_pred  <- as.numeric(cbind(1, X_tot) %*% m1$coef)
    t1 <- totals_boot
    t2 <- totals_boot
    t1$pred <- exp(lin_pred) / (1 + exp(lin_pred))
    t2$pred <- 0
    t1$flag <- "1"
    t2$flag <- "0"
    
    totals_cal <- as.data.frame(rbind(t1,t2))
    totals_cal$wt  <- ifelse(totals_cal$flag  == 1, totals_cal$pred*totals_cal$hat_wolne, totals_cal$hat_wolne)
    tab_totals <- xtabs(wt~ rok + flag, data = totals_cal)
    tab_totals[,1] <- tab_totals[,1] - tab_totals[,2]
    final_data_kod2$flag <- as.character(final_data_kod2[,komps[k]])
    lasso_kod2 <- svydesign(ids = ~1, weight = ~ weight, data = final_data_kod2)
    cal_lasso_kod2 <- calibrate(design = lasso_kod2, 
                                formula = list(~rok + flag), 
                                population = list(tab_totals))
    
    svyby(formula = as.formula(paste("~", komps[k])), 
          by = ~ rok, 
          FUN = svymean, 
          design = cal_lasso_kod2) %>%
      select(-se) -> wyn
    
    wynik_est[[k]] <- wyn
    wynik_wt[[k]] <- weights(cal_lasso_kod2)
    wynik_model[[k]] <- m1
  }
  
  
  wynik <- list(wynik_est = bind_cols(wynik_est), 
                wynik_wt = bind_cols(wynik_wt),
                wynik_model = wynik_model)
  
  
}

toc <- Sys.time() - tic
toc
