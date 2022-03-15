####################################################################################################################
## Table S1: in each region for the summer season: regression of NIRv on NO2  and regression of NIRV on albedo
####################################################################################################################


rm(list=ls())
library(fixest)
require(tidyverse)
library(here)
library(sjPlot)
here()


df_maize<-read.csv(here("data","final_allvars_dfs_filteredModisPeak","samples100kUSA_maize_allvars_allCountries_filteredpoints.csv"))



saveR2 <- function(df){
  

  N <- 20 
  DF <- data.frame(CNTRY_NAME=rep(NA, N), crop=rep(NA, N), year_gws=rep(NA, N), greenness_var=rep(NA, N), 
                   Adjusted.R2_NIRv_vs_albedo=rep(NA, N), Within.R2_NIRv_vs_albedo=rep(NA, N),
                   Adjusted.R2_NIRv_vs_NO2=rep(NA, N), Within.R2_NIRv_vs_NO2=rep(NA, N),
                   stringsAsFactors=FALSE) 
  
  country_vars <- as.character(df$CNTRY_NAME) %>% unique()
  crop_vars <- as.character(df$crop) %>% unique()
  this.vi <- 'NIRv_max'
  
  i<-1
  for (this.country in country_vars){
 
    print(this.country)
    
    dat <- df %>% filter(CNTRY_NAME==this.country)
    allyears<-paste(dat$year_gws %>% unique(), collapse = '-')
    
    
    if(nrow(dat)>1){
      

        linform_fe = as.formula('NIRv_max ~ OMLER440nm| cell_geomid + year_gws')
        linmod_fe = feols(linform_fe,data=dat)

        ar2_NIRv_vs_albedo <- r2(linmod_fe, type = "ar2",full_names=TRUE)
        wr2_NIRv_vs_albedo <- r2(linmod_fe, type = "wr2",full_names=TRUE)

        linform_fe2 = as.formula('NIRv_max ~ NO2| cell_geomid + year_gws')
        linmod_fe2 = feols(linform_fe2,data=dat)

        ar2_NIRv_vs_NO2 <- r2(linmod_fe2, type = "ar2",full_names=TRUE)
        wr2_NIRv_vs_NO2 <- r2(linmod_fe2, type = "wr2",full_names=TRUE)

        
        DF[i, ] <- list(this.country, crop_vars, allyears, this.vi,
                        ar2_NIRv_vs_albedo, wr2_NIRv_vs_albedo,
                        ar2_NIRv_vs_NO2,    wr2_NIRv_vs_NO2
                        )
        i<-i+1
    }else{
      
      print(paste0(this.country," ",this.year," ####   EMPTY DF   ####"))
    }
  }
  
  DF <- DF %>% drop_na(CNTRY_NAME)
  return(DF)
}


DF_withR2<-saveR2(df_maize)