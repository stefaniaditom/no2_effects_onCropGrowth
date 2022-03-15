#########################################################################################
### Plot NO2 coeff vs No2averages by Region/ years - (fig S4)
## scatter anomalies to explain Fixed Effects - (fig 3 (c))
#########################################################################################

rm(list=ls())
library(fixest)
require(tidyverse)
library(here)
library("data.table")
library(ggridges) 
library(sjPlot)
here()


df_maize<-read.csv(here("data","final_allvars_dfs_filteredModisPeak","samples100kUSA_maize_allvars_allCountries_filteredpoints.csv"))
df_wheat<-read.csv(here("data","final_allvars_dfs_filteredModisPeak","samples100kUSA_wheat_allvars_allCountries_filteredpoints.csv"))
df_rice <-read.csv(here("data","final_allvars_dfs_filteredModisPeak","samples100kUSA_ricet_allvars_allCountries_filteredpoints.csv"))




process.overall <- function(df,weath_control){
  
  if(weath_control==TRUE){
    control_form = " ~ pr_early + pr_early_squared + pr_late + pr_late_squared + vpd + NO2 | cell_geomid"
  }else if(weath_control==FALSE){
    control_form = " ~ NO2 | cell_geomid"
  }
  
  this.regime <- "overall"
  o3_control <-FALSE
  
  #Create empty DF
  N <- 100  # total number of rows to preallocate
  DF <- data.frame(CNTRY_NAME=rep(NA, N), crop=rep(NA, N), year_gws=rep(NA, N), greenness_var=rep(NA, N), regime=rep(NA, N),
                   NO2coef=rep(NA, N), NO2pval=rep(NA, N), no2stdErr=rep(NA, N), no2CI_2.5 = rep(NA, N), no2CI_97.5 = rep(NA, N),
                   NO2avg=rep(NA, N), NO2p05=rep(NA, N),Npoints=rep(NA, N),
                   greennVarAvg = rep(NA, N),
                   weath_control=rep(NA, N), o3_control=rep(NA, N),
                   stringsAsFactors=FALSE)
  
  
  country_vars <- as.character(df$CNTRY_NAME) %>% unique()
  year_vars <- df$year_gws %>% unique()
  crop_vars <- as.character(df$crop) %>% unique()
  vi_vars <- c('NDVI_max','NIRv_max')
  
  i<-1
  for (this.country in country_vars){
    
    for (this.year in year_vars){
      
      # this.year<-2019
      dat <- df %>% filter(CNTRY_NAME==this.country) %>% filter(year_gws==this.year)
      # print(dim(dat))
      
      
      if(nrow(dat)>1){
        
        for (this.vi in vi_vars){
          print(paste0(this.country," ",this.year," ",this.vi))
          
          linform_fe = as.formula(paste0(this.vi,control_form))
          linmod_fe = feols(linform_fe,data=dat)

          
          no2cf = unname(linmod_fe$coefficients['NO2']) 
          no2stdErr = unname(summary(linmod_fe, se = "cluster")$se['NO2'])
          no2pval = (summary(linmod_fe)$coeftable)['NO2','Pr(>|t|)']
          ci = confint(linmod_fe, 'NO2', level=0.95) # adding 95% Confidence Interval
          no2CI_2.5 = unname(unlist(ci[1]))
          no2CI_97.5 = unname(unlist(ci[2]))
          
          # save NO2 average and 5th percentile
          NO2avg <- dat$NO2 %>% mean(na.rm = TRUE)
          NO2p05 <- quantile(dat$NO2, 0.05)
          greennVarAvg<- dat[[this.vi]] %>% mean(na.rm = TRUE)
          
          Npoints <- nrow(dat)
          
          DF[i, ] <- list(this.country,crop_vars , this.year, this.vi, this.regime,
                          no2cf, no2pval, no2stdErr, no2CI_2.5,no2CI_97.5,
                          NO2avg, NO2p05, Npoints, greennVarAvg,
                          weath_control, o3_control
          )
          
          i<-i+1
          
        }
        
      }else{
        print(paste0(this.country," ",this.year," ####   EMPTY DF   ####"))
      }
      
    }
    
  }
  DF <- DF %>% drop_na(CNTRY_NAME)
  DF$signif <- ifelse(DF$NO2pval < 0.05, '< 0.05', '>= 0.05')
  DF<- DF %>% group_by(CNTRY_NAME,year_gws,greenness_var) %>% mutate(NpointsPerc = Npoints/sum(Npoints))
  
  
  return(DF)
}




outSplit.maize_overall<- process.overall(df_maize, TRUE)
outSplit.wheat_overall <- process.overall(df_wheat, TRUE)
outSplit.rice_overall <- process.overall(df_rice, TRUE)

outSplit.overall<-rbind(outSplit.maize_overall,outSplit.wheat_overall,outSplit.rice_overall)
outSplit.overall<- outSplit.overall %>% mutate(country_year = paste0(CNTRY_NAME,year_gws))




# Make zeros print as "0" always
prettyZero <- function(l){
  max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1
  lnew = formatC(l, replace.zero = T, zero.print = "0",
                 digits = max.decimals, format = "f", preserve.width=T)
  return(lnew)
}



plot_NO2coeffs_v2<-function(Dataframe,vi,shapes,fname){
  
  # # Dataframe<-outAll.overall
  Dataframe$crop[Dataframe$crop=='maize'] <- 'summer'
  Dataframe$crop[Dataframe$crop=='wheat'] <- 'winter'
  Dataframe$crop[Dataframe$crop=='rice'] <- 'rice'
  
  Dataframe<-Dataframe %>%filter(greenness_var==vi)
  Dataframe<-Dataframe %>%filter(NpointsPerc>0.01)
  
  ggplot(data = Dataframe, mapping = aes(x = NO2coef, y = c(country_year)))+
    geom_point(aes(colour=crop, size=NpointsPerc, shape=signif),alpha=1,stroke = 2)+
    geom_errorbar(aes(xmin=NO2coef-no2stdErr, xmax=NO2coef+no2stdErr), width=.0002)+
    ylab(NULL)+
    xlab(expression(NO[2]~coefficient))+
    geom_vline(xintercept=0,color="grey")+
    scale_x_continuous(labels = prettyZero)+
    guides(size=FALSE,colour = guide_legend(order = 1), 
           shape = guide_legend(order = 2))+
    labs(shape="p-value",colour="season")+
    scale_shape_manual(values=c(16,1))+ 
    scale_colour_manual(values=c("red","goldenrod2", "royalblue1")) 

  
  # ggsave(
  #   fname,
  #   plot = last_plot(),
  #   device = "png",
  #   path = here("figs"),
  #   scale = 1,
  #   width = 18,
  #   height = 7,
  #   units = "cm",
  #   dpi = 200,
  #   limitsize = TRUE
  # )
  
}

plot_NO2coeffs_v2(outSplit.overall,'NIRv_max',c(16),"fig.png")



###################################################################################
## Figure S4: PLOT NO2 coeff vs No2averages by Region/ years
###################################################################################


# Plot coeff vs no2 averages by region
plot_coeffvsAvg <- function(df,fname){
  
  df$year_gws <- as.factor(df$year_gws)
  df$crop[df$crop=='maize'] <- 'summer'
  df$crop[df$crop=='wheat'] <- 'winter'
  df$crop[df$crop=='rice'] <- 'rice'
  df$crop = factor(df$crop, levels=c('winter','summer','rice'))
  
  df$greenness_var[df$greenness_var=='NDVI_max'] <- 'NDVI'
  df$greenness_var[df$greenness_var=='NIRv_max'] <- 'NIRv'
  
  ggplot(data = df, mapping = aes(x = NO2avg, y = NO2coef))+
    geom_point(aes(color = CNTRY_NAME, shape = year_gws),size=3)+
    geom_errorbar(aes(ymin=no2CI_2.5, ymax=no2CI_97.5), width=.0002)+ 
    facet_grid(greenness_var ~ crop) +
    expand_limits(y=0)+
    geom_hline(yintercept=0,color="grey")+
    ylab(expression(NO[2]~coefficient))+
    xlab(expression(NO[2]~average))+
    labs(shape="",color = "")+
    # guides(fill=guide_legend(title=""))+
    theme_bw()
    
  ggsave(
    fname,
    plot = last_plot(),
    device = "png",
    path = here("figs"),
    scale = 1,
    width = 20,
    height = 9,
    units = "cm",
    dpi = 200,
    limitsize = TRUE
  )
    
}

plot_coeffvsAvg(outSplit.overall,'figS4_no2coeffVSno2Average.png')
#############################################################################




#######################################################################
## scatter anomalies to explain Fixed Effects - Figure 3 (c)
## CHINA - 2020 - wheat
# only plot for China winter season- to use to explain LL FE in paper
#######################################################################


scatter_Anomalies <- function(df, fname){
  

  df <- df %>%
    group_by_at(c('crop','CNTRY_NAME','year_gws','cell_geomid')) %>%
    mutate_at(.vars = vars(NIRv_max, NO2,pr_early,pr_early_squared,pr_late,pr_late_squared,vpd,AOD_055,AOD_047), .funs = funs('dm' = . - mean(.)))
  
  
  control_form = "NIRv_max ~ pr_early + pr_early_squared + pr_late + pr_late_squared + vpd + NO2 | cell_geomid + year_gws"
  linform_fe = as.formula(control_form)
  linmod_fe = feols(linform_fe,data=df)
  
  
  control_form_dm = "NIRv_max_dm ~ pr_early_dm + pr_early_squared_dm + pr_late_dm + pr_late_squared_dm + vpd_dm + NO2_dm "
  linform_fe = as.formula(control_form_dm)
  linmod_fe = feols(linform_fe,data=df)

  
  ggplot(df, aes(x=NO2_dm, y=NIRv_max_dm)) +
    geom_point(alpha = 0.05)+
    geom_smooth(method=lm)+
    ylab('NIRv de-meaned')+
    xlab(expression(NO[2]~de-meaned))+
    theme_bw()
  
  ggsave(
    fname,
    plot = last_plot(),
    device = "png",
    path = here("figs"),
    scale = 1,
    width = 12,
    height = 10,
    units = "cm",
    dpi = 600,
    limitsize = TRUE
  )
}


scatter_Anomalies(df_wheat%>% filter(CNTRY_NAME=='China') %>% filter(year_gws==2020),"fig3_scatterAnomalies_China_2020winter.png")





