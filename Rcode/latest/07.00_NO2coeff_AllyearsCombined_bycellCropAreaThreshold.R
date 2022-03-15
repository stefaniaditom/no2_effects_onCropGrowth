#################################################################################################
### Figure S3. Restricting analysis to areas with a high proportion of cropland.
#################################################################################################


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





process.overall <- function(df, weath_control){
  
  if(weath_control==TRUE){
    control_form = " ~ pr_early + pr_early_squared + pr_late + pr_late_squared + vpd + NO2 | cell_geomid + year_gws"
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
  DF <- DF %>% drop_na(CNTRY_NAME)
  DF$signif <- ifelse(DF$NO2pval < 0.05, '< 0.05', '>= 0.05')

  DF<- DF %>% group_by(CNTRY_NAME,greenness_var) %>% mutate(NpointsPerc = Npoints/sum(Npoints)) 
  
  return(DF)
}




df_maize$cell_cropArea %>%max() 
df_rice$cell_cropArea %>%max() 


ths = seq(0, 2900000000, by=500000000)
datalist_overall = list()



for (i in 1:length(ths)){

  th <-ths[i]
  print(th)
  df_maize_sub<-df_maize %>% filter(cell_cropArea>=th)
  df_wheat_sub<-df_wheat %>% filter(cell_cropArea>=th)

  outAll.maize_overall<- process.overall(df_maize_sub, TRUE)
  outAll.maize_overall$th <- th
  outAll.wheat_overall <- process.overall(df_wheat_sub, TRUE)
  outAll.wheat_overall$th <- th

  # only maize and wheat
  outAll.overall_temp<-rbind(outAll.maize_overall, outAll.wheat_overall)
  outAll.overall_temp$i <- i
  datalist_overall[[i]] <- outAll.overall_temp

}

outAll.overall = do.call(rbind, datalist_overall)





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


plot_NO2coeffs_vs_th<-function(Dataframe,vi,shapes,fname){
  
  Dataframe$crop[Dataframe$crop=='maize'] <- 'summer'
  Dataframe$crop[Dataframe$crop=='wheat'] <- 'winter'
  Dataframe$crop[Dataframe$crop=='rice'] <- 'rice'
  
  Dataframe<-Dataframe %>%filter(greenness_var==vi)
  Dataframe<-Dataframe %>%filter(NpointsPerc>0.01)
  Dataframe<- Dataframe%>%mutate(CNTRY_NAME=factor(CNTRY_NAME, levels = c("Western Europe", "United States", "South America","India","China")))

  
  ggplot(data = Dataframe, mapping = aes(x = NO2coef, y = th))+
    geom_point(aes(colour=crop, size=NpointsPerc, shape=signif),alpha=1,stroke = 2)+
    geom_errorbar(aes(xmin=no2CI_2.5, xmax=no2CI_97.5), width=.0002)+ # using 95% CI
    facet_grid(vars(CNTRY_NAME))+
    ylab('crop area threshold (m2/gridcell)')+
    xlab(expression(NO[2]~coefficient))+
    geom_vline(xintercept=0,color="grey")+

    coord_cartesian(ylim = c(-300000000, 2700000000))+
    scale_x_continuous(labels = prettyZero)+

    guides(size=FALSE,
           colour = guide_legend(order = 1,reverse=TRUE),
           shape = guide_legend(order = 2)
           )+
    labs(shape="p-value",colour="season")+
    scale_shape_manual(values=c(16,1))+
    scale_colour_manual(values=c("goldenrod2", "royalblue1"))+ 
    theme_bw()
  
  ggsave(
    fname,
    plot = last_plot(),
    device = "png",
    path = here("figs"),
    scale = 1,
    width = 12,
    height = 15,
    units = "cm",
    dpi = 600,
    limitsize = TRUE
  )
}

plot_NO2coeffs_vs_th(outAll.overall,'NIRv_max',c(16),"S3_no2coeff_vs_cropAreaTh.png")








plot_Npoints_vs_th<-function(Dataframe,vi,shapes,fname){
  
  Dataframe$crop[Dataframe$crop=='maize'] <- 'summer'
  Dataframe$crop[Dataframe$crop=='wheat'] <- 'winter'
  Dataframe$crop[Dataframe$crop=='rice'] <- 'rice'
  
  Dataframe<-Dataframe %>%filter(greenness_var==vi)
  Dataframe<-Dataframe %>%filter(NpointsPerc>0.01)
  
  Dataframe<- Dataframe%>%mutate(CNTRY_NAME=factor(CNTRY_NAME, levels = c("Western Europe", "United States", "South America","India","China")))
  
  ggplot(data = Dataframe, mapping = aes(x = Npoints, y = th))+
    geom_point(aes(colour=crop, size=NpointsPerc, shape=signif),alpha=1,stroke = 2)+
    facet_grid(vars(CNTRY_NAME))+
    ylab('crop area threshold (m2/gridcell)')+
    xlab('Number of points')+
    geom_vline(xintercept=0,color="grey")+
    coord_cartesian(ylim = c(-300000000, 2700000000))+
    guides(size=FALSE,
           colour = guide_legend(order = 1,reverse=TRUE),
           shape = guide_legend(order = 2)
    )+
    labs(shape="p-value",colour="season")+
    scale_shape_manual(values=c(16,1))+ 
    scale_colour_manual(values=c("goldenrod2", "royalblue1"))+ 
    theme_bw()
  
  ggsave(
    fname,
    plot = last_plot(),
    device = "png",
    path = here("figs"),
    scale = 1,
    width = 12,
    height = 15,
    units = "cm",
    dpi = 600,
    limitsize = TRUE
  )
}

plot_Npoints_vs_th(outAll.overall,'NIRv_max',c(16),"S3_npoints_vs_cropAreaTh.png")




