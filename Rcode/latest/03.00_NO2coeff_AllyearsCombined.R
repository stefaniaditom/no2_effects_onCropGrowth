#########################################################################################
### Plot NO2 coeff all Years Combined WITH WEATHER  overall and by Regime (fig 4,5)
### Plot Yield percentage change by Region (fig 6)
#########################################################################################


rm(list=ls())
library(fixest)
require(tidyverse)
library(here)
library("data.table")
library(ggridges)
library(sjPlot)
# library("ggstance")
here()


df_maize<-read.csv(here("data","final_allvars_dfs_filteredModisPeak","samples100kUSA_maize_allvars_allCountries_filteredpoints.csv"))
df_wheat<-read.csv(here("data","final_allvars_dfs_filteredModisPeak","samples100kUSA_wheat_allvars_allCountries_filteredpoints.csv"))
df_rice <-read.csv(here("data","final_allvars_dfs_filteredModisPeak","samples100kUSA_ricet_allvars_allCountries_filteredpoints.csv"))


###################################################################################
## Figure 4,5: Plot NO2 coeff for overall and by splitting by regime ALL YEARS COMBINED
###################################################################################



process.overall <- function(df, weath_control){
  
  if(weath_control==TRUE){
    control_form = " ~ pr_early + pr_early_squared + pr_late + pr_late_squared + vpd + NO2 | cell_geomid + year_gws"
  }else if(weath_control==FALSE){
    control_form = " ~ NO2 | cell_geomid + year_gws"
  }
  
  this.regime <- "overall"
  o3_control <-FALSE
  
  #Create empty DF
  N <- 100  # total number of rows to preallocate--possibly an overestimate
  DF <- data.frame(CNTRY_NAME=rep(NA, N), crop=rep(NA, N), year_gws=rep(NA, N), greenness_var=rep(NA, N), regime=rep(NA, N),
                   NO2coef=rep(NA, N), NO2pval=rep(NA, N), no2stdErr=rep(NA, N), no2CI_2.5 = rep(NA, N), no2CI_97.5 = rep(NA, N),
                   NO2avg=rep(NA, N), NO2p05=rep(NA, N),Npoints=rep(NA, N),
                   greennVarAvg = rep(NA, N),
                   # O3avg=rep(NA, N), ratioavg=rep(NA, N),O3tropavg=rep(NA, N), #NO2sd=rep(NA, N), O3sd=rep(NA, N), ratiosd=rep(NA, N), O3tropsd=rep(NA, N),
                   # monthlyRatioAvg=rep(NA, N), , 
                   weath_control=rep(NA, N), o3_control=rep(NA, N),
                   # ,inputdir=rep(NA, N),
                   stringsAsFactors=FALSE) 
  
  
  country_vars <- as.character(df$CNTRY_NAME) %>% unique()
  # year_vars <- df$year_gws %>% unique()
  crop_vars <- as.character(df$crop) %>% unique()
  vi_vars <- c('NDVI_max','NIRv_max')
  
  i<-1
  for (this.country in country_vars){
    # print(this.country)
    
    # for (this.year in year_vars){
    
    # print(this.year)
    dat <- df %>% filter(CNTRY_NAME==this.country) #%>% filter(year_gws==this.year)
    allyears<-paste(dat$year_gws %>% unique(), collapse = '-')
    # print(dim(dat))
    
    
    if(nrow(dat)>1){
      
      for (this.vi in vi_vars){
        # print(this.vi)
        # print(paste0(this.country," ",this.year," ",this.vi))
        print(paste0(this.country," ",this.vi))
        
        linform_fe = as.formula(paste0(this.vi,control_form))
        # print(linform_fe)
        linmod_fe = feols(linform_fe,data=dat)
        # print(linmod_fe)
        
        no2cf = unname(linmod_fe$coefficients['NO2']) ##(linmod_fe$coeftable)['NO2','Estimate']
        no2stdErr = unname(summary(linmod_fe, se = "cluster")$se['NO2'])
        no2pval = (summary(linmod_fe)$coeftable)['NO2','Pr(>|t|)']
        # print(no2cf)
        # print(no2stdErr)
        # print(no2pval)
        ci = confint(linmod_fe, 'NO2', level=0.95) # adding 95% Confidence Interval
        no2CI_2.5 = unname(unlist(ci[1]))
        no2CI_97.5 = unname(unlist(ci[2]))
        # print(no2CI_2.5)
        # print(typeof(no2CI_2.5))
 
        
        
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
                        # O3avg, ratioavg, O3tropavg,# NO2sd, O3sd, ratiosd, O3tropsd,monthlyRatioAvg,
        )
        
        i<-i+1
        
      }
      
    }else{
      print(paste0(this.country," ",this.year," ####   EMPTY DF   ####"))
      # print('#################   EMPTY DF   ##################################')
    }
    
    # }
    
  }
  DF <- DF %>% drop_na(CNTRY_NAME)
  DF$signif <- ifelse(DF$NO2pval < 0.05, '< 0.05', '>= 0.05')
  # DF<- DF %>% group_by(CNTRY_NAME,year_gws,greenness_var) %>% mutate(NpointsPerc = Npoints/sum(Npoints)) #to check if works for multiple regimes
  DF<- DF %>% group_by(CNTRY_NAME,greenness_var) %>% mutate(NpointsPerc = Npoints/sum(Npoints)) #to check if works for multiple regimes
  
  
  return(DF)
}





process.byregime <- function(df, weath_control){

  if(weath_control==TRUE){
    control_form = " ~ pr_early + pr_early_squared + pr_late + pr_late_squared + vpd + NO2 | cell_geomid + year_gws"
  }else if(weath_control==FALSE){
    control_form = " ~ NO2 | cell_geomid + year_gws"
  }

  o3_control <-FALSE

  #Create empty DF
  N <- 100  # total number of rows to preallocate--possibly an overestimate
  DF <- data.frame(CNTRY_NAME=rep(NA, N), crop=rep(NA, N), year_gws=rep(NA, N), greenness_var=rep(NA, N), regime=rep(NA, N),
                   NO2coef=rep(NA, N), NO2pval=rep(NA, N), no2stdErr=rep(NA, N), no2CI_2.5 = rep(NA, N), no2CI_97.5 = rep(NA, N),
                   NO2avg=rep(NA, N), NO2p05=rep(NA, N),Npoints=rep(NA, N),
                   greennVarAvg = rep(NA, N),
                   # O3avg=rep(NA, N), ratioavg=rep(NA, N),O3tropavg=rep(NA, N), #NO2sd=rep(NA, N), O3sd=rep(NA, N), ratiosd=rep(NA, N), O3tropsd=rep(NA, N),
                   # monthlyRatioAvg=rep(NA, N), ,
                   weath_control=rep(NA, N), o3_control=rep(NA, N),
                   # ,inputdir=rep(NA, N),
                   stringsAsFactors=FALSE)


  country_vars <- as.character(df$CNTRY_NAME) %>% unique()
  # year_vars <- df$year_gws %>% unique()
  crop_vars <- as.character(df$crop) %>% unique()
  vi_vars <- c('NDVI_max','NIRv_max')
  regime_vars <- c("low","high")


  i<-1
  for (this.country in country_vars){

    # for (this.year in year_vars){

    for (this.regime in regime_vars){

      dat <- df %>% filter(CNTRY_NAME==this.country) %>% filter(regime==this.regime)
      allyears<-paste(dat$year_gws %>% unique(), collapse = '-')
      # print(dim(dat))

      if(nrow(dat)>1){

        for (this.vi in vi_vars){
          # print(this.vi)
          print(paste0(this.country," ",this.regime," ",this.vi))

          linform_fe = as.formula(paste0(this.vi,control_form))
          # print(linform_fe)
          
          tryCatch( # to make it work even if feols fails beacuse of collinaraity- will skip that and keep on next coruntry/regime
                    # this usually happens when regime has few points
            expr = {
              
              linmod_fe = feols(linform_fe,data=dat)
              # print(linmod_fe)
              
              no2cf = unname(linmod_fe$coefficients['NO2']) ##(linmod_fe$coeftable)['NO2','Estimate']
              no2stdErr = unname(summary(linmod_fe, se = "cluster")$se['NO2'])
              no2pval = (summary(linmod_fe)$coeftable)['NO2','Pr(>|t|)']
              # print(no2cf)
              # print(no2stdErr)
              # print(no2pval)
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
                              # O3avg, ratioavg, O3tropavg,# NO2sd, O3sd, ratiosd, O3tropsd,monthlyRatioAvg,
              )
              i<-i+1
              
            },
            
            error = function(e){ 
              # (Optional)
              # Do this if an error is caught...
              print('Error Here')
            },
            warning = function(w){
              # (Optional)
              # Do this if an warning is caught...
            },
            finally = {
              # (Optional)
              # Do this at the end before quitting the tryCatch structure...
            }
          )
            
        }

      }else{
        print(paste0(this.country," ",this.regime," ####   EMPTY DF   ####"))
        # print('#################   EMPTY DF   ##################################')
      }
    }#close regime
    # }#close year
  }#close country

  DF <- DF %>% drop_na(CNTRY_NAME)
  DF$signif <- ifelse(DF$NO2pval < 0.05, '< 0.05', '>= 0.05')
  # DF<- DF %>% group_by(CNTRY_NAME,year_gws,greenness_var) %>% mutate(NpointsPerc = Npoints/sum(Npoints)) #to check if works for multiple regimes
  DF<- DF %>% group_by(CNTRY_NAME,greenness_var) %>% mutate(NpointsPerc = Npoints/sum(Npoints)) #to check if works for multiple regimes

  return(DF)
}






outAll.maize_overall<- process.overall(df_maize, TRUE)
outAll.maize_multi <- process.byregime(df_maize, TRUE)

outAll.wheat_overall <- process.overall(df_wheat, TRUE)
outAll.wheat_multi <- process.byregime(df_wheat, TRUE)

outAll.rice_overall <- process.overall(df_rice, TRUE)
outAll.rice_multi <- process.byregime(df_rice, TRUE)


###################################################################################
## Fig 6:  TABLE Yield percentage change                                 ##
###################################################################################

NIRv0 <- 0.07

yield_perc.maize <- outAll.maize_overall %>% filter(greenness_var=='NIRv_max')
yield_perc.maize$NO2Delta <- yield_perc.maize$NO2p05 - yield_perc.maize$NO2avg
yield_perc.maize$NIRvGain <- yield_perc.maize$NO2coef * yield_perc.maize$NO2Delta
yield_perc.maize$percYieldChange <-  yield_perc.maize$NIRvGain / (yield_perc.maize$greennVarAvg - NIRv0) * 100
yield_perc.maize$percYieldChange_stderror <- yield_perc.maize$no2stdErr * abs(yield_perc.maize$NO2Delta / (yield_perc.maize$greennVarAvg - NIRv0) * 100)
yield_perc.maize$percYieldChange_no2CI_2.5 <-  abs(yield_perc.maize$NO2coef - yield_perc.maize$no2CI_2.5) * abs(yield_perc.maize$NO2Delta / (yield_perc.maize$greennVarAvg - NIRv0) * 100)
yield_perc.maize$percYieldChange_no2CI_97.5 <- abs(yield_perc.maize$NO2coef - yield_perc.maize$no2CI_97.5) * abs(yield_perc.maize$NO2Delta / (yield_perc.maize$greennVarAvg - NIRv0) * 100)



yield_perc.wheat <- outAll.wheat_overall %>% filter(greenness_var=='NIRv_max')
yield_perc.wheat$NO2Delta <- yield_perc.wheat$NO2p05 - yield_perc.wheat$NO2avg
yield_perc.wheat$NIRvGain <- yield_perc.wheat$NO2coef * yield_perc.wheat$NO2Delta
yield_perc.wheat$percYieldChange <- yield_perc.wheat$NIRvGain / (yield_perc.wheat$greennVarAvg - NIRv0) * 100
yield_perc.wheat$percYieldChange_stderror <- yield_perc.wheat$no2stdErr * abs(yield_perc.wheat$NO2Delta / (yield_perc.wheat$greennVarAvg - NIRv0) * 100)
yield_perc.wheat$percYieldChange_no2CI_2.5 <- abs(yield_perc.wheat$NO2coef - yield_perc.wheat$no2CI_2.5) * abs(yield_perc.wheat$NO2Delta / (yield_perc.wheat$greennVarAvg - NIRv0) * 100)
yield_perc.wheat$percYieldChange_no2CI_97.5 <-abs(yield_perc.wheat$NO2coef - yield_perc.wheat$no2CI_97.5) * abs(yield_perc.wheat$NO2Delta / (yield_perc.wheat$greennVarAvg - NIRv0) * 100)


yield_perc<-rbind(yield_perc.maize,yield_perc.wheat)

yield_perc_small <- yield_perc[,c("CNTRY_NAME","crop", "year_gws","NO2coef","no2stdErr","NO2avg","NO2p05","NO2Delta","NIRvGain","greennVarAvg", "percYieldChange","percYieldChange_stderror","percYieldChange_no2CI_2.5","percYieldChange_no2CI_97.5")]

yield_perc_small<- yield_perc_small %>% 
  rename(
    Region = CNTRY_NAME,
    NO2stdError = no2stdErr,
    NIRvavg = greennVarAvg
  )



bar_yieldchange <-function(df,fname){
  # df<-df%>% mutate(country_crop = paste0(Region,'_',crop))
  df$crop[df$crop=='maize'] <- 'summer'
  df$crop[df$crop=='wheat'] <- 'winter'
  df<- df%>%mutate(Region=factor(Region, levels = c("Western Europe", "United States", "South America","India","China")))
  # ggplot(data = df, mapping = aes(x=country_crop,y=percYieldChange, fill=crop))+
  ggplot(data = df, mapping = aes(x=Region,y=percYieldChange, fill=crop))+
    geom_bar(stat="identity", position=position_dodge())+
    # scale_fill_manual(values=c("goldenrod2", "royalblue1"))+
    # scale_fill_manual(values=c("goldenrod2", "royalblue1"))+
    # scale_fill_manual(values=c("#DC3220", "#005AB5"))+
    # scale_fill_manual(values=c("#ffc20a", "#0c7bdc"))+
    scale_fill_manual(values=c("#00B680", "#007CA6"))+ #green/blue
    # coord_cartesian(ylim = c(-2, 40))+
    ylab("Yield Change %")+
    xlab(NULL)+
    labs(fill=NULL)+
    guides(fill = guide_legend(reverse=TRUE))+ # reverse order legend: first winter then summer
    # geom_errorbar(aes(ymin=percYieldChange-percYieldChange_stderror, ymax=percYieldChange+percYieldChange_stderror), width=.2,
    #               position=position_dodge(.9))+
    geom_errorbar(aes(ymin=percYieldChange-percYieldChange_no2CI_2.5, ymax=percYieldChange+percYieldChange_no2CI_97.5), width=.2,
                  position=position_dodge(.9))+
    
    coord_flip()+
    theme_bw()
  
  
  # ggsave(
  #   fname,
  #   plot = last_plot(),
  #   device = "png",
  #   path = here("figs"),
  #   scale = 1,
  #   width = 18,
  #   height = 6,
  #   units = "cm",
  #   dpi = 200,
  #   limitsize = TRUE
  # )
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
bar_yieldchange(yield_perc_small,'fig6_barplotYieldChange.png')



# ###################################################################################################################################









# only maize and wheat
outAll.overall<-rbind(outAll.maize_overall, outAll.wheat_overall)
outAll.multi<-rbind(outAll.maize_multi,outAll.wheat_multi)


# # test adding RICE!
# outAll.overall<-rbind(outAll.maize_overall,outAll.wheat_overall,outAll.rice_overall)
# outAll.multi<-rbind(outAll.maize_multi,outAll.wheat_multi,outAll.rice_multi)



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



# One Panel for Fig 4- NIRv and Summer and Winter color coded
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

plot_NO2coeffs_v2(outAll.overall,'NIRv_max',c(16),"fig4_No2coeffAllyears.png")






# Fig 5
plot_NO2coeffs_v2<-function(Dataframe,vi,shapes,fname){

  Dataframe$crop[Dataframe$crop=='maize'] <- 'summer'
  Dataframe$crop[Dataframe$crop=='wheat'] <- 'winter'
  Dataframe$crop[Dataframe$crop=='rice'] <- 'rice'
  
  Dataframe<-Dataframe %>%filter(greenness_var==vi)
  Dataframe<-Dataframe %>%filter(NpointsPerc>0.01)

  Dataframe<- Dataframe%>%mutate(CNTRY_NAME=factor(CNTRY_NAME, levels = c("Western Europe", "United States", "South America","India","China")))
  
  Dataframe$crop = factor(Dataframe$crop, levels=c('winter','summer','rice'))
  Dataframe$regime = factor(Dataframe$regime, levels=c('low','high'))
  
  
  # ggplot(data = Dataframe, mapping = aes(x = NO2coef, y = CNTRY_NAME))+
  #   geom_point(aes(colour=regime, size=NpointsPerc, shape=signif),alpha=1,stroke=2)+
  #   geom_errorbar(aes(xmin=no2CI_2.5, xmax=no2CI_97.5), width=.0002)+ # using 95% CI
  #   facet_wrap(facets = vars(crop))+#, scales = "free_x")
  #   ylab(NULL)+
  #   xlab(expression(NO[2]~coefficient))+
  #   geom_vline(xintercept=0,color="grey")+
  #   coord_cartesian(xlim = c(-0.005, 0.001))+
  #   scale_x_continuous(labels = prettyZero)+
  #   guides(size=FALSE,# get rid of Npoints legend
  #          colour = guide_legend(order = 1), 
  #          shape = guide_legend(order = 2)
  #   )+
  #   labs(shape="p-value",colour="regime")+
  #   scale_shape_manual(values=c(16,1))+ # shapes need to be from21-25 to have a fill and color property
  #   theme_bw()
  
  
  
  
  ggplot(data = Dataframe, mapping = aes(y = NO2coef, x = CNTRY_NAME,colour=regime, ymin=no2CI_2.5, ymax=no2CI_97.5))+
    geom_point(aes(size=NpointsPerc, shape=signif),alpha=1,stroke=2, position=position_dodge(width=0.5))+
    geom_errorbar(width=.0002, position=position_dodge(width=0.5))+ # using 95% CI
    facet_wrap(facets = vars(crop))+#, scales = "free_x")
    xlab(NULL)+
    ylab(expression(NO[2]~coefficient))+
    geom_hline(yintercept=0,color="grey")+
    coord_cartesian(ylim = c(-0.005, 0.001))+
    scale_y_continuous(labels = prettyZero)+
    guides(size=FALSE,# get rid of Npoints legend
           colour = guide_legend(order = 1), 
           shape = guide_legend(order = 2)
    )+
    labs(shape="p-value",colour="regime")+
    scale_shape_manual(values=c(16,1))+ # shapes need to be from21-25 to have a fill and color property
    theme_bw()+coord_flip()

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
library('svglite')
plot_NO2coeffs_v2(outAll.multi,'NIRv_max',c(15,17),"fig5_No2coeffAllyears_byregime.png")

