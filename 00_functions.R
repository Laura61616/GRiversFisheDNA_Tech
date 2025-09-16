#############################################################
# loading libraries
#############################################################
## reading and formatting data
library(readxl)
library(tidyr)
library(dplyr)
#library(plyr)
library(reshape2)
library(purrr)
## species standardize
#library(rfishbase)
library(taxize)
## plotting
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggpmisc)
library(ggspatial)
library(ggtext)
library(hrbrthemes)
## geo mapping
library(sf)
library(terra)
## biodiversity
library(vegan)
library(ape)
library(betapart)
## modelling
library(lmodel2)
library(lmerTest)
library(mgcv)
library(report)
library(partR2)
library(rstatix)
library(FSA)

## climate space
library(paran)
library(factoextra)
library(FactoMineR)
library(corrplot)
## 
library(decontam)
#############################################################
# loading functions
#############################################################
#########generating label from model
gam_fit_lab <- function(fit, fac, explained = NULL){
  p <- fit$p.table[fac,'Pr(>|t|)']
  r <- fit$r.sq
  if (is.null(explained)) {
    expl <- fit$dev.expl
  }else{
    expl <- explained
  }
  ind <- case_when(p < 0.001 ~ '***',
                   (p >= 0.001)&(p < 0.01) ~ '**', 
                   (p >= 0.01)&(p < 0.05) ~ '*', 
                   (p >= 0.05)&(p < 0.1) ~ '.', 
                   p >= 0.1 ~ '')
  if (is.nan(p)) {
    lab <- paste0('R2 = ', round(r,2), 
                  '\nexplained = ', round(expl*100, 2),
                  '%\np = NA')
    return(lab)
    
  }
  if (p < 0.001) {
    lab <- paste0('R2 = ', round(r,2), 
                  '\nexplained = ', round(expl*100, 2),
                  '%\np < 0.001', ' ',ind)
  }else{
    lab <- paste0('R2 = ', round(r,2), 
                  '\nexplained = ', round(expl*100, 2),
                  '%\np = ', round(p,3), ' ',ind)
  }
  return(lab)
}
lme_fit_lab <- function(fit, fac, explained = NULL){
  p <- fit$coefficients[fac,'Pr(>|t|)']
  r <- fit$AICtab
  if (is.null(explained)) {
    expl <- NULL
  }else{
    expl <- explained
  }
  ind <- case_when(p < 0.001 ~ '***',
                   (p >= 0.001)&(p < 0.01) ~ '**', 
                   (p >= 0.01)&(p < 0.05) ~ '*', 
                   (p >= 0.05)&(p < 0.1) ~ '.', 
                   p >= 0.1 ~ '')
  if (p < 0.001) {
    lab <- paste0('AIC = ', round(r,2), 
                  '\nexplained = ', round(expl*100, 2),
                  '%\np < 0.001', ' ',ind)
  }else{
    lab <- paste0('AIC = ', round(r,2), 
                  '\nexplained = ', round(expl*100, 2),
                  '%\np = ', round(p,3), ' ',ind)
  }
  return(lab)
}
lmod2_fit_lab <- function(fit, method = 'OLS'){
  tmp <- subset(fit$regression.results, Method == method)
  p <- tmp$`P-perm (1-tailed)`
  slope <- tmp$Slope
  r <- fit$rsquare
  ind <- case_when(p < 0.001 ~ '***',
                   (p >= 0.001)&(p < 0.01) ~ '**', 
                   (p >= 0.01)&(p < 0.05) ~ '*', 
                   (p >= 0.05)&(p < 0.1) ~ '.', 
                   p >= 0.1 ~ '')
  if (p < 0.001) {
    lab <- paste0('R2 = ', round(r,2), 
                  '\nslope = ', round(slope, 2),
                  '\np < 0.001', ' ',ind)
  }else{
    lab <- paste0('R2 = ', round(r,2), 
                  '\nslope = ', round(slope, 2),
                  '\np = ', round(p,3), ' ',ind)
  }
  return(lab)
}
summary_lme <- function(rep, r, useFactors, other_para = 'Continent'){
  est <- as.data.frame(r$Ests)
  rownames(est) <- est$term
  est <- est[c(useFactors, grep(other_para, rownames(est), value = T)),]
  
  explained <- as.data.frame(r$R2)
  rownames(explained) <- explained$term
  explained <- explained[c(useFactors, other_para),]
  
  rep <- rep[!is.na(rep$Parameter),]
  rownames(rep) <- gsub(' ','_',rep$Parameter)
  cond_r2 <- na.omit(rep[rep$Parameter=='R2 (conditional)','Fit'])
  if (length(cond_r2) == 0) {
    cond_r2 <- NA
  }
  
  paras <- c(useFactors, grep(other_para,rownames(rep),value = T))
  paras[grep('Continent',paras)[1]] <- 'Continent'
  full_res <- data.frame(
    'AIC' = rep[rep$Parameter=='AIC','Fit'],
    'cond_r2' = cond_r2,
    'marg_r2' = rep[rep$Parameter=='R2 (marginal)','Fit'],
    rep[c(useFactors, grep(other_para,rownames(rep),value = T)),],
    'estimate' = est[, 'estimate'],
    'expained' = explained[paras, 'estimate']
  )
  return(full_res)
}
