################################################################################################################
## Table S2. Characteristics of regions used in the study                                                      
## Data were gathered using Google Earth Engine. GEE scripts for exporting the data are available upon request
################################################################################################################


rm(list=ls())
library(fixest)
require(tidyverse)
library(here)
library("data.table")
library(sjPlot)
here()


df_maize<-read.csv(here("data","final_allvars_dfs_filteredModisPeak","samples100kUSA_maize_allvars_allCountries_filteredpoints.csv"))
df_wheat<-read.csv(here("data","final_allvars_dfs_filteredModisPeak","samples100kUSA_wheat_allvars_allCountries_filteredpoints.csv"))
df_rice <-read.csv(here("data","final_allvars_dfs_filteredModisPeak","samples100kUSA_ricet_allvars_allCountries_filteredpoints.csv"))



# merge all files 
df_raw <- bind_rows(list(df_maize,df_wheat,df_rice))


# info for Table S2 in paper: Number of points per Region/Season
df_raw %>%group_by(CNTRY_NAME,crop)%>% summarise(n = n())
table_w_years_used<-df_raw %>%group_by(CNTRY_NAME,crop, year_gws)%>% summarise(n = n())




