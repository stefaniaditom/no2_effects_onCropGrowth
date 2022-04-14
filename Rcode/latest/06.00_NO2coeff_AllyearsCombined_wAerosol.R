################################################################################################################################
### Plot NO2 coeff all Years Combined WITH WEATHER  overall adding aerosol as control (fig S2)
#############################################################################################################################


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


###################################################################################
## Figure S2
###################################################################################


process.overall <- function(df, weath_control){
  
  if(weath_control==TRUE){
    control_form = " ~ pr_early + pr_early_squared + pr_late + pr_late_squared + vpd + NO2 + AOD_055| cell_geomid + year_gws"
  }else if(weath_control==FALSE){
    control_form = " ~ NO2 | cell_geomid + year_gws"
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
  crop_vars <- as.character(df$crop) %>% unique()
  vi_vars <- c('NDVI_max','NIRv_max')
  
  i<-1
  for (this.country in country_vars){

    dat <- df %>% filter(CNTRY_NAME==this.country) 
    allyears<-paste(dat$year_gws %>% unique(), collapse = '-')

    
    
    if(nrow(dat)>1){
      
      for (this.vi in vi_vars){
        print(paste0(this.country," ",this.vi))
        
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
        
        #save number of points used in the regression after all the filtering
        Npoints <- nrow(dat)
        
        DF[i, ] <- list(this.country, crop_vars, allyears, this.vi, this.regime,
                        no2cf, no2pval, no2stdErr,no2CI_2.5,no2CI_97.5,
                        NO2avg, NO2p05, Npoints, greennVarAvg,
                        weath_control, o3_control

        )
        
        i<-i+1
        
      }
      
    }else{
      print(paste0(this.country," ",this.year," ####   EMPTY DF   ####"))

    }
    

    
  }
  DF <- DF %>% drop_na(CNTRY_NAME)
  DF$signif <- ifelse(DF$NO2pval < 0.05, '< 0.05', '>= 0.05')
  
  DF<- DF %>% group_by(CNTRY_NAME,greenness_var) %>% mutate(NpointsPerc = Npoints/sum(Npoints))
  
  
  return(DF)
}



outAll.maize_overall<- process.overall(df_maize, TRUE)
outAll.wheat_overall <- process.overall(df_wheat, TRUE)
outAll.rice_overall <- process.overall(df_rice, TRUE)

outAll.overall<-rbind(outAll.maize_overall, outAll.wheat_overall)



# Make zeros print as "0" always
prettyZero <- function(l){
  max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1
  lnew = formatC(l, replace.zero = T, zero.print = "0",
                 digits = max.decimals, format = "f", preserve.width=T)
  return(lnew)
}

crop_names <- c(
  "maize" = "summer crop",
  "wheat" = "winter crop",
  "rice" = "rice summer crop"
)


plot_NO2coeffs_v2<-function(Dataframe,vi,shapes,fname){
  
  Dataframe$crop[Dataframe$crop=='maize'] <- 'summer'
  Dataframe$crop[Dataframe$crop=='wheat'] <- 'winter'
  Dataframe$crop[Dataframe$crop=='rice'] <- 'rice'

  Dataframe<-Dataframe %>%filter(greenness_var==vi)
  Dataframe<-Dataframe %>%filter(NpointsPerc>0.01)
  
  Dataframe<- Dataframe%>%mutate(CNTRY_NAME=factor(CNTRY_NAME, levels = c("Western Europe", "United States", "South America","India","China")))

  ggplot(data=Dataframe, mapping=aes(y= NO2coef, x = CNTRY_NAME, ymin=no2CI_2.5, ymax=no2CI_97.5, colour=crop))+
   geom_point(aes(size=NpointsPerc, shape=signif),alpha=1,stroke=2, position=position_dodge(width=0.3))+
    geom_errorbar(width=.002, position=position_dodge(0.3))+ # using 95% CI
    coord_flip()+
    xlab(NULL)+
    ylab(expression(NO[2]~coefficient))+
    geom_hline(yintercept=0,color="grey")+
    scale_y_continuous(labels = prettyZero)+
    guides(size=FALSE,
           colour = guide_legend(order = 1,reverse=TRUE),
           shape = guide_legend(order = 2)
    )+# get rid of Npoints legend
    labs(shape="p-value",colour="season")+
    scale_shape_manual(values=c(16,1))+ # shapes need to be from21-25 to have a fill and color property
    scale_colour_manual(values=c("#00B680", "#007CA6"))+ #green/blue
    theme_bw()

  ggsave(
    fname,
    plot = last_plot(),
    device = "png",
    path = here("figs"),
    scale = 1,
    width = 18,
    height = 6,
    units = "cm",
    dpi = 200,
    limitsize = TRUE
  )
  
}

plot_NO2coeffs_v2(outAll.overall,'NIRv_max',c(16),"figS2_No2coeffAllyears_wAerosolControl.png")



