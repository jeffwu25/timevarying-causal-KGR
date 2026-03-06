#load("~/work/P2-Causal-Analysis-NB-total-v2-10.29.RData") #cum6 workspace
load("~/work/causal-wildfire-smoke-2-duplicate/P2-Causal-Analysis-NB-totalcum2-2.7.RData") #cum2 workplace

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
library(corrplot)
library(superheat)
library(shapes)
library(scales)
library(splines)
library(fmesher)
library(utils)
library(tictoc)
library(parallel)
library(foreach)
library(doParallel)
library(future)
library(future.apply)
library(RhpcBLASctl)
library(pryr)

#Need to set threads bc INLA uses a lot of threads, can lead to issues in parallel computing 
Sys.setenv(
  OPENBLAS_NUM_THREADS = "10",
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

#Setup: 
inla_df$cum_monthly_total_pm25_2_binned = inla.group(inla_df$cum_monthly_total_pm25_2,n=25,method = "quantile")
inla_df$Cluster = as.numeric(inla_df$Cluster)
inla_df$Group = as.numeric(inla_df$Group)
inla_df$sw_binned = inla.group(inla_df$sw_prod,n=10,method = "quantile")
# inla_df$pm25_dev = inla_df$cum_monthly_total_pm25_2_binned
inla_df$pm25_dev = inla.group(inla_df$cum_monthly_total_pm25_2,n=10,method = "cut")


#Need to make predictions in batches - otherwise they are too unstable
bins = list(c(1:15),c(16:30),c(31:45),c(46:58))

#Estimating DRF function (same)
DRF_curves = function(dataset,dose_grid,sample_id){
  cf_preds_list = list()
  
  for(i in 1:length(dose_grid)){
    cat("Calculating DRF for dose level:", dose_grid[i], "\n")
    
    # Initialize combined predictions dataframe for this dose
    combined_cf_preds = NULL
    
    for (j in 1:4){
      cat("  Processing batch", j, "of 4\n")
      
      # Get counties for this batch
      batch_counties = counties[bins[[j]]]
      
      # Create duplicate dataset for just these counties
      duplicate_df = dataset[dataset$NAME %in% batch_counties,]
      duplicate_df$cum_monthly_total_pm25_2 = dose_grid[i]
      duplicate_df$deaths = NA
      
      # Combine original dataset with counterfactual data for this batch
      inla_cf_df = rbind(dataset, duplicate_df)
      cf_idx = c((nrow(dataset)+1):nrow(inla_cf_df))
      
      # kgr_formula_hierarchical = deaths ~ offset(log(Total_Pop/100000)) +
      #   f(NAME,model = "iid",hyper = list(prec = list(prior = "loggamma", param = c(10, 5))))+
      #   factor(month) + 
      #   aqi_resids +
      #   cum_total_precip_6 +
      #   avg_avg_air_temp_6 +
      #   avg_avg_dew_point_6 +
      #   max_avg_wind_speed_6 +
      #   f(sw_binned,model = "rw2") +
      #   # Spline for treatment * SDI group 
      #   Group*ns(cum_monthly_total_pm25_6_binned,df=3) +
      #   f(ID,model = "generic0",Cmatrix = solve(H^2), group = time, control.group = list(model = "ar",order = 6))
      
      kgr_formula_hierarchical2 = deaths ~ offset(log(Total_Pop/100000)) +
        f(NAME,model = "iid",hyper = list(prec = list(prior = "loggamma", param = c(10, 5)))) +
        factor(month) + 
        sw_prod + I(sw_prod^2) + I(sw_prod^3) + 
        # Spline for treatment * SDI group 
        Group*ns(cum_monthly_total_pm25_2_binned,df=3) +
        f(ID,model = "generic0",Cmatrix = solve(H^2), group = time, control.group = list(model = "ar",order = 6))
      
      kgr_formula_hierarchical3 <- deaths ~ offset(log(Total_Pop/100000)) +
        factor(month) + 
        sw_prod2 + I(sw_prod2^2) + I(sw_prod2^3) + #use sw_prod2 for lognormal model PSs
        f(NAME, model = "iid",constr=TRUE) + 
        ns(cum_monthly_total_pm25_2,df=3) + 
        f(pm25_dev, model = "rw2", replicate = Group, scale.model = TRUE,
          hyper = list(prec = list(prior = "pc.prec", param = c(0.05, 0.5)))) +
        f(ID, model = "generic0", Cmatrix = solve(H^2), 
          group = time, control.group = list(model = "ar", order = 6))
      
      # Error handling for inla()
      kgr_model <- tryCatch({
        inla(formula = kgr_formula_hierarchical3, family = "nbinomial", data = inla_cf_df, num.threads = 1,#verbose = TRUE,
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
      
      # Extract predictions for counterfactual data
      preds_model <- kgr_model$summary.fitted.values
      preds_model <- cbind(sample_id,inla_cf_df$NAME, inla_cf_df$time, inla_cf_df$date, 
                           inla_cf_df$cum_monthly_total_pm25_2, inla_cf_df$deaths, preds_model)
      colnames(preds_model) <- c("sample_id","county", "time", "date", "exposure", "deaths", "mean", 
                                 "sd", "0.025quant", "0.5quant", "0.975quant", "mode")
      
      # Get predictions for this batch's counterfactual data
      batch_cf_preds = preds_model[cf_idx,]
      
      # Append to combined predictions
      if (is.null(combined_cf_preds)) {
        combined_cf_preds = batch_cf_preds
      } else {
        combined_cf_preds = rbind(combined_cf_preds, batch_cf_preds)
      }
    }
    
    # Store combined predictions for this dose level
    cf_preds_list[[i]] = combined_cf_preds
  }
  
  return(cf_preds_list)
}

#Wild bootstrap implementation function: 
wild_bootstrap_DRF_curves = function(sample_id){
  #Generate random multipliers
  n_obs = nrow(fitted_values)
  #multipliers = rnorm(n_obs,0,1)
  multipliers = sample(c(-1,1),n_obs,replace = TRUE) #draw from Rademacher distribution 
  
  #Create wild bootstrap pseudo-responses
  wild_resids = multipliers*resids
  bootstrap_deaths = round(fitted_values$mean + wild_resids)
  bootstrap_deaths = pmax(0,bootstrap_deaths) #any counties with all 0s??
  
  #Create bootstrap dataset
  bootstrap_df = inla_df
  bootstrap_df$deaths = bootstrap_deaths
  
  #Estimate DRF from bootstrap dataset
  sample_DRF = DRF_curves(bootstrap_df,dose_grid,sample_id)
  saveRDS(sample_DRF, paste0("cum2_bootstrap_sample_", sample_id, ".rds"))
  
  gc()
  return(sample_DRF)
}

# Set up parallel processing (each core will handle one sample DRF curve)
plan(multisession, workers = 25) #set to number of samples in a batch
blas_set_num_threads(1)
set.seed(2)

#Import fitted values from final outcome model 
fitted_values = final_model_results$fitted_values
resids = (fitted_values$deaths - round(fitted_values$mean))

#Define dose grid

# dose_grid = c(0,10,20,50,100,200,500,1000,2000,3000)
# dose_grid = c(0,10,20,50,100,200,400,600,800,1000) #cum6
# dose_grid = c(0,0.1,0.5,1,5,10,50,100,500,1000) #cum2

cum2_quantiles = quantile(inla_df$cum_monthly_total_pm25_2,probs = c(seq(0.1,0.9,by=0.1),0.95))
#cum3_quantiles = quantile(inla_df$cum_monthly_total_pm25_3,probs = c(seq(0.1,0.9,by=0.1),0.95))
#cum6_quantiles = quantile(inla_df$cum_monthly_total_pm25_6,probs = c(seq(0.1,0.9,by=0.1),0.95))

dose_grid = round(cum2_quantiles,3)


#We had to run in batches of 25 samples at a time given computational resources 
tic()

batch_indices = c(75:100)
batch_results = future_lapply(batch_indices, wild_bootstrap_DRF_curves, future.seed = 2L) #25 samples takes about 3 days to run on 384 cores

# saveRDS(batch_results, "cum2_bootstrap_results_samples_1_to_25.rds")
# saveRDS(batch_results, "cum2_bootstrap_results_samples_26_to_50.rds")
# saveRDS(batch_results, "cum2_bootstrap_results_samples_51_to_75.rds")
saveRDS(batch_results, "cum2_bootstrap_results_samples_76_to_100.rds")

toc()
