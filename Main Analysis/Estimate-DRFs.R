# load("~/work/causal-wildfire-smoke-2-duplicate/P2-Causal-Analysis-NB-totalcum2-1.14.RData")
load("~/work/causal-wildfire-smoke-2-duplicate/P2-Causal-Analysis-NB-totalcum2-2.7.RData") #ILI prevalence in PS model
# load("~/work/causal-wildfire-smoke-2-duplicate/P2-Causal-Analysis-NB-totalcum2-hugegraph-1.28.RData") #graph estimated with HUGE
# load("~/work/causal-wildfire-smoke-2-duplicate/P2-Causal-Analysis-NB-totalcum2-gfilter1-2.13.RData")  #0.33 gfilter cutoff
# load("~/work/causal-wildfire-smoke-2-duplicate/P2-Causal-Analysis-NB-totalcum2-gfilter3-2.14.RData")  #0.99 gfilter cutoff
# load("~/work/P2-Causal-Analysis-NB-totalcum12-8.12.RData")

set.seed(2)

library(dplyr) 
library(tidyverse)
library(lubridate)
library(stringr)
library(zoo)
library(ggplot2)
library(urbnmapr)
library(devtools)
library(readxl)
library(spdep)
library(sp)
library(huge)
library(INLA)
library(HMMpa)
library(invgamma)
library(brinla)
library(reshape2)
library(patchwork)
library(jsonlite)
library(geosphere)
library(urbnmapr)
library(RAQSAPI)
library(con2aqi)
library(pscl)
library(reshape2)
#library(corrplot)
#library(superheat)
#library(shapes)
library(scales)
library(splines)
library(fmesher)

#Defining dose grid: 

# dose_grid = c(0,10,50,200,500,1000,2000,3000)
# dose_grid = c(0,10,20,50,100,200,500,1000,2000,3000)
# dose_grid = c(0,0.5,1,10,20,50,100,500,1000,2500)
# dose_grid = c(0,0.1,0.5,1,5,10,50,100,500,1000)

cum2_quantiles = quantile(inla_df$cum_monthly_total_pm25_2,probs = c(seq(0.1,0.9,by=0.1),0.95))
print(cum2_quantiles)
dose_grid = round(cum2_quantiles,3)

#Setup: 

#inla_df$cum_monthly_total_pm25_2_binned = inla.group(inla_df$cum_monthly_total_pm25_2,n=25,method = "quantile")
inla_df$cum_monthly_total_pm25_2_binned = inla.group(inla_df$cum_monthly_total_pm25_2,n=10,method = "cut")
inla_df$Cluster = as.numeric(inla_df$Cluster)
inla_df$Group = as.numeric(inla_df$Group)
inla_df$sw_binned = inla.group(inla_df$sw_prod,n=10,method = "quantile")
inla_df$pm25_dev = inla_df$cum_monthly_total_pm25_2_binned
#inla_df$pm25_dev = inla.group(inla_df$cum_monthly_total_pm25_2,n=10,method = "cut")

#Need to make predictions in batches - otherwise they are too unstable
bins = list(c(1:15),c(16:30),c(31:45),c(46:58))

DRF_curves = function(dataset,dose_grid){
  cf_preds_list = list()
  
  for(i in 1:length(dose_grid)){
    cat("Calculating DRF for dose level:", dose_grid[i], "\n")
    
    # Initialize combined predictions dataframe for this dose
    combined_cf_preds = NULL
    
    for (j in 1:4){
      cat("  Processing batch", j, "of 4\n")
      
      #Get counties for this batch
      batch_counties = counties[bins[[j]]]
      
      #Create duplicate dataset for just these counties with dose grid z
      duplicate_df = dataset[dataset$NAME %in% batch_counties,]
      duplicate_df$cum_monthly_total_pm25_2 = dose_grid[i]
      duplicate_df$deaths = NA
      
      # Combine original dataset with counterfactual data for this batch
      inla_cf_df = rbind(dataset, duplicate_df)
      cf_idx = c((nrow(dataset)+1):nrow(inla_cf_df))
      
      # kgr_formula_hierarchical3 = deaths ~ -1 +
      #   Intercept_1 + Intercept_2 + Intercept_3 + Intercept_4 +
      #   Intercept_5 + Intercept_6 + Intercept_7 + Intercept_8 +
      #   Intercept_9 + Intercept_10 + Intercept_11 + Intercept_12 +
      #   Intercept_13 + Intercept_14 + Intercept_15 + Intercept_16 +
      #   Intercept_17 + Intercept_18 + Intercept_19 + Intercept_20 +
      #   Intercept_21 + Intercept_22 + Intercept_23 + Intercept_24 +
      #   Intercept_25 + Intercept_26 + Intercept_27 + Intercept_28 +
      #   Intercept_29 + Intercept_30 + Intercept_31 + Intercept_32 +
      #   Intercept_33 + Intercept_34 + Intercept_35 + Intercept_36 +
      #   Intercept_37 + Intercept_38 + Intercept_39 + Intercept_40 +
      #   Intercept_41 + Intercept_42 + Intercept_43 + Intercept_44 +
      #   Intercept_45 + Intercept_46 + Intercept_47 + Intercept_48 +
      #   Intercept_49 + Intercept_50 + Intercept_51 + Intercept_52 +
      #   Intercept_53 + Intercept_54 + Intercept_55 + Intercept_56 + Intercept_57 + Intercept_58 +
      #   factor(month) +
      #   aqi_resids +
      #   cum_total_precip_2 +
      #   avg_avg_air_temp_2 +
      #   avg_avg_dew_point_2 +
      #   max_avg_wind_speed_2 +
      #   f(sw_binned,model = "rw2") +
      #   # Global PM2.5 dose-response curve
      #   ns(cum_monthly_total_pm25_2,df=3) +
      #   # Cluster-specific deviations from global curve (shrunken toward zero)
      #   f(Z_cluster, model = spde, group = Cluster, control.group = list(model = "exchangeable")) +
      #   f(ID,model = "generic0",Cmatrix = solve(H^2), group = time, control.group = list(model = "ar",order = 6))
      
      kgr_formula_hierarchical = deaths ~ offset(log(Total_Pop/100000)) +
        f(NAME,model = "iid",hyper = list(prec = list(prior = "loggamma", param = c(10, 5)))) +
        factor(month) + 
        sw_prod2 + I(sw_prod2^2) + I(sw_prod2^3) + 
        # Spline for treatment * SDI group 
        Group*ns(cum_monthly_total_pm25_2,df=3) +
        f(ID,model = "generic0",Cmatrix = solve(H^2), group = time, control.group = list(model = "ar",order = 6))
      
      kgr_formula_hierarchical2 <- deaths ~ offset(log(Total_Pop/100000)) +
        factor(month) + 
        sw_prod2 + I(sw_prod2^2) + I(sw_prod2^3) + #use sw_prod2 for lognormal model PSs
        f(NAME, model = "iid",constr=TRUE) + 
        ns(cum_monthly_total_pm25_2,knots = c(10,50,100,200)) + 
        f(pm25_dev, model = "rw2", replicate = Group, scale.model = TRUE,
          hyper = list(prec = list(prior = "pc.prec", param = c(0.05, 0.5)))) +
        f(ID, model = "generic0", Cmatrix = solve(H^2), 
          group = time, control.group = list(model = "ar", order = 6))
      
      kgr_formula_hierarchical3 <- deaths ~ offset(log(Total_Pop/100000)) +
        factor(month) + 
        f(time,model = "iid",hyper = 
            list(prec = list(prior = "pc.prec",param = c(0.1,0.1))),constr=TRUE) + 
        sw_prod2 + I(sw_prod2^2) + I(sw_prod2^3) +
        f(NAME, model = "iid",constr=TRUE) + 
        ns(cum_monthly_total_pm25_2,df=3) + 
        f(pm25_dev, model = "rw2", replicate = Group, scale.model = TRUE,
          hyper = list(prec = list(prior = "pc.prec", param = c(0.05, 0.5)))) +
        f(ID, model = "generic0", Cmatrix = solve(H^2), 
          group = time, control.group = list(model = "ar", order = 6))
      
      #Fit model with error handling for inla()
      kgr_model <- tryCatch({
        inla(formula = kgr_formula_hierarchical2, family = "nbinomial", data = inla_cf_df, num.threads = 20,
             control.compute = list(dic=TRUE,waic=TRUE,return.marginals.predictor=TRUE),
             #control.fixed = prior_fixed,
             control.family = list(hyper = list(theta = list(prior = "pc.mgamma", param = 1))),
             control.inla = list(strategy = "simplified.laplace"),
             control.predictor = list(compute = TRUE, link = 1))
      }, error = function(e) {
        cat("Error fitting INLA model at dose", dose_grid[i], "batch", j, ":", conditionMessage(e), "\n")
        return(NULL)
      })
      
      if (is.null(kgr_model)) {
        cat("Skipping batch", j, "due to model fitting error\n")
        next
      }
      
      #Extract predictions for counterfactual data
      preds_model <- kgr_model$summary.fitted.values
      preds_model <- cbind(inla_cf_df$NAME, inla_cf_df$time, inla_cf_df$date, 
                           inla_cf_df$cum_monthly_total_pm25_2, inla_cf_df$deaths, preds_model)
      colnames(preds_model) <- c("county", "time", "date", "exposure", "deaths", "mean", 
                                 "sd", "0.025quant", "0.5quant", "0.975quant", "mode")
      
      #Get predictions for this batch's counterfactual data
      batch_cf_preds = preds_model[cf_idx,]
      
      #Append to combined predictions
      if (is.null(combined_cf_preds)) {
        combined_cf_preds = batch_cf_preds
      } else {
        combined_cf_preds = rbind(combined_cf_preds, batch_cf_preds)
      }
    }
    
    #Store combined predictions for this dose level
    cf_preds_list[[i]] = combined_cf_preds
  }
  
  return(cf_preds_list)
}

full_data_DRFs = DRF_curves(inla_df,dose_grid) #takes about 8 hrs to run on 384 cores

#Save results
saveRDS(full_data_DRFs,"full_data_drfs_NB_totalcum2_nocov_newspline_grouprw2_polyswprod_iliconfounder_ar6xgfilter_2.17.rds")
