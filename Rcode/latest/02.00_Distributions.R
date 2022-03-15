#########################################################################################
### NO2 Distribution (fig 2 and fig S7)
### NO2 Distribution - Anomalies - for 0.5 and 1 deg FE  (fig S1)
### Ratio distribution (fig 5)
### NH3 Distribution (fig S5)
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


###################################################################################
## Figure 2 an S7:  NO2 Distribution
###################################################################################


plot_NO2Dist <- function(df,fname){
  
  df<- df %>% mutate(country_year = paste0(CNTRY_NAME,year_gws))
  df$year_gws<- as.factor(df$year_gws)
  # country_names = as.character(df$CNTRY_NAME) %>% unique()
  # years = df$year_gws %>% unique()
  crop_type = df$crop %>% unique()
  if(crop_type=='maize'){
    title='summer crop'
  }else if(crop_type=='wheat'){
    title='winter crop'
  }else if(crop_type=='rice'){
    title='Rice summer crop'
  }
  
  ggplot(df, aes(x=NO2, y=fct_rev(as_factor(country_year)), fill=year_gws)) +
    geom_density_ridges(scale=3)+
    facet_grid(vars(CNTRY_NAME), scales = "free",switch="y")+
    coord_cartesian(xlim = c(0, 100), ylim = c(1,5.5))+
    ylab(NULL)+
    xlab('NO2 Tropospheric column (1e-6 mol/m2)')+
    # xlab(expression(NO[2]~Tropospheric~column~(1e^-6~mol/m^2)))+
    ggtitle(title)+
    # theme(axis.title.x = element_text(angle=0, hjust = 0.5))+
    # theme(plot.title = element_text(hjust = 0.5))+
    # scale_fill_manual(values = c("gray30","gray50","gray70"))+
    scale_fill_manual(values = c("#F8766D","#00BA38","#619CFF"))+
    geom_text(aes(label = CNTRY_NAME), x = 80, y = Inf, vjust=1.5) +
    guides(fill=guide_legend(title=""))+
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'gray90')
    )
  
  ggsave(
    fname,
    plot = last_plot(),
    device = "png",
    path = here("figs"),
    scale = 1,
    width = 12,
    height = 10,#maize,wheat
    # height = 4.8, #rice
    units = "cm",
    dpi = 600,
    limitsize = TRUE
  )
}

plot_NO2Dist(df_maize,"figS7_b_summer.png")
plot_NO2Dist(df_wheat,"fig2_b_winter.png")
# plot_NO2Dist(df_rice,"fig2_rice.png")



###################################################################################
## Figure S1:  NO2 Distribution - Anomalies
###################################################################################

plot_NO2AnomaliesDist <- function(df, LL_FE, fname){
  
  # df: maize or wheat
  # LL_FE: if using 1 or 1/2 degree FE
  
  df <- df %>%
    # group_by(crop, CNTRY_NAME, year_gws, LL_FE) %>%
    group_by_at(c('crop','CNTRY_NAME','year_gws',LL_FE)) %>%
    mutate_at(.vars = vars(NDVI_max,NIRv_max, NO2), .funs = funs('dm' = . - mean(.)))
  
  df<- df %>% mutate(country_year = paste0(CNTRY_NAME,year_gws))
  df$year_gws<- as.factor(df$year_gws)
  # country_names = as.character(df$CNTRY_NAME) %>% unique()
  # years = df$year_gws %>% unique()
  crop_type = df$crop %>% unique()
  if(crop_type=='maize'){
    title='summer crop'
  }else if(crop_type=='wheat'){
    title='winter crop'
  }
  if(LL_FE=='LL_grid'){
    title2='1/2 degree FE'
  }else if(LL_FE=='LL1d_grid'){
    title2='1 degree FE'
  }
  
  ggplot(df, aes(x=NO2_dm, y=fct_rev(as_factor(country_year)), fill=year_gws)) +
    geom_density_ridges(scale=3)+
    facet_grid(vars(CNTRY_NAME), scales = "free",switch="y")+
    coord_cartesian(xlim = c(-20, 20), ylim = c(1,5.5))+
    ylab(NULL)+
    xlab('NO2 Tropospheric column (1e-6 mol/m2)')+
    # xlab(bquote('NO2 Tropospheric column (1'~e^-6 ~ 'mol/' ~ m^2 ~')' ))+
    ggtitle(paste0(title,' ',title2))+
    scale_fill_manual(values = c("#F8766D","#00BA38","#619CFF"))+
    # theme(axis.title.x = element_text(angle=0, hjust = 0.5))+
    # theme(plot.title = element_text(hjust = 0.5))+
    # guides(fill=guide_legend(title=""))+
    # theme(strip.text.y.left = element_text(angle = 0),
    #       axis.title.y=element_blank(),
    #       axis.text.y=element_blank(),
    #       axis.ticks.y=element_blank(),
    #       panel.background = element_rect(fill = 'gray96', colour = 'gray90'))
    
    geom_text(aes(label = CNTRY_NAME), x = 15, y = Inf, vjust=1.5) +
    
    guides(fill=guide_legend(title=""))+
    theme(
      # strip.text.y.left = element_text(angle = 0),
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'gray90')
      # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray95")
      # panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "red")
    )
  
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


plot_NO2AnomaliesDist(df_maize,'LL_grid',"figS1_NO2anomalies_0.5dFEsummer.png")
plot_NO2AnomaliesDist(df_wheat,'LL_grid',"figS1_NO2anomalies_0.5dFEwinter.png")

plot_NO2AnomaliesDist(df_maize,'LL1d_grid',"figS1_NO2anomalies_1dFEsummer.png")
plot_NO2AnomaliesDist(df_wheat,'LL1d_grid',"figS1_NO2anomalies_1dFEwinter.png")



###################################################################################
## Figure 5:  RATIO Distribution
###################################################################################

# By year - only one year in each plot
plot_RatioDist_byYear <- function(df,fname){
  
  
  df<-df %>%drop_na(monthlyRationHCHO_NO2)
  df$year_gws<- as.factor(df$year_gws)
  crop_type = df$crop %>% unique()
  if(crop_type=='maize'){
    title='summer crop'
  }else if(crop_type=='wheat'){
    title='winter crop'
  }else if(crop_type=='rice'){
    title='rice crop'
  }
  
  # df$high_ratio <- as.factor((df[,'monthlyRationHCHO_NO2']>2)*1)
  df$high_ratio <- (df[,'monthlyRationHCHO_NO2']>2)*1
  df$high_ratio[df$high_ratio==0] <- 'low'
  df$high_ratio[df$high_ratio==1] <- 'high'
  df$high_ratio <- as.factor(df$high_ratio)
  df$high_ratio <- factor(df$high_ratio, levels = c("low", "high"))
  
  # HISTOGRAM ##
  ggplot() +
    geom_histogram(data = df, aes(x=monthlyRationHCHO_NO2, fill=high_ratio),bins = 200)+
    # facet_grid(vars(CNTRY_NAME,year_gws), scales = "free")+
    facet_grid(vars(CNTRY_NAME), scales = "free")+
    # facet_grid(vars(CNTRY_NAME))+
    # coord_cartesian(xlim = c(-1, 15))+
    scale_x_continuous(limits = c(0, 10))+
    geom_vline(xintercept=2, colour="gray")+
    ylab('Number of sample points')+
    # xlab('HCHO/NO2 ratio')+
    xlab(expression(HCHO/NO[2]))+
    ggtitle(title)+
    labs(fill = "")+#remove legend title+
  
  
    geom_text(data = df, mapping = aes(label = CNTRY_NAME), x = 8, y = Inf, vjust=1.5) +
    # 
    # guides(fill=guide_legend(title=""))
    theme(
      # strip.text.y.left = element_text(angle = 0),
      strip.background = element_blank(),
      strip.text = element_blank(),
      # axis.title.y=element_blank(),
      # axis.text.y=element_blank(),
      # axis.ticks.y=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'gray90')
    )
  ggsave(
    fname,
    plot = last_plot(),
    device = "png",
    path = here("figs"),
    scale = 1,
    width = 12,
    height = 10,
    # height = 4.8, #rice
    units = "cm",
    dpi = 600,
    limitsize = TRUE
  )
  
}
plot_RatioDist_byYear(df_maize %>% filter(year_gws==2020),"fig5_ratiohist_summer2020.png")
plot_RatioDist_byYear(df_wheat %>% filter(year_gws==2020),"fig5_ratiohist_winter2020.png")
# plot_RatioDist_byYear(df_rice %>% filter(year_gws==2020),"fig4_rice2020.png")

# plot_RatioDist_byYear(df_maize %>% filter(year_gws==2019),"fig5_summer2019.png")
# plot_RatioDist_byYear(df_wheat %>% filter(year_gws==2019),"fig5_winter2019.png")
# plot_RatioDist_byYear(df_rice %>% filter(year_gws==2019),"fig4_rice2019.png")







###################################################################################
## Figure S5:  NH3 Distribution
###################################################################################



plot_NH3Dist_v2 <- function(df,fname){
  
  df<-df %>%drop_na(NH3)
  df<- df %>% mutate(country_year = paste0(CNTRY_NAME,year_gws))
  df$year_gws<- as.factor(df$year_gws)
  # country_names = as.character(df$CNTRY_NAME) %>% unique()
  # years = df$year_gws %>% unique()
  crop_type = df$crop %>% unique()
  if(crop_type=='maize'){
    title='summer crop'
  }else if(crop_type=='wheat'){
    title='winter crop'
  }else if(crop_type=='rice'){
    title='rice crop'
  }
  
  ggplot(df, aes(x=NH3, y=fct_rev(as_factor(country_year)), fill=year_gws)) +
    geom_density_ridges(scale=3)+
    facet_grid(vars(CNTRY_NAME), scales = "free",switch="y")+
    # coord_cartesian(xlim = c(-1, 40))+
    coord_cartesian(xlim = c(-1, 40), ylim = c(1,4.5))+
    ylab(NULL)+
    # xlab('NO2 Tropospheric column (1e-6 mol/m2)')+
    # xlab(expression(NO[2]~Tropospheric~column~(1e^-6~mol/m^2)))+
    ggtitle(title)+
    # scale_fill_manual(values = c("gray30","gray50","gray70"))+
    scale_fill_manual(values = c("#F8766D","#00BA38","#619CFF"))+
    # theme(axis.title.x = element_text(angle=0, hjust = 0.5))+
    # theme(plot.title = element_text(hjust = 0.5))+
    
    geom_text(aes(label = CNTRY_NAME), x = 35, y = Inf, vjust=1.5) +
    guides(fill=guide_legend(title=""))+
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'gray90')
    )
  
  ggsave(
    fname,
    plot = last_plot(),
    device = "png",
    path = here("figs"),
    scale = 1,
    width = 12,
    height = 10,#maize,wheat
    # height = 4.8, #rice
    units = "cm",
    dpi = 600,
    limitsize = TRUE
  )
}
plot_NH3Dist_v2(df_maize,"figS3_summer.png")
plot_NH3Dist_v2(df_wheat,"figS3_winter.png")
plot_NH3Dist_v2(df_rice,"figS3_rice.png")


############################################################################################

