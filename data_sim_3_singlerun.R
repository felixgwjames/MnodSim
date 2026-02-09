# Stochastic recruitment

#set.seed(12345)

#rm(list = ls())

setwd("/home/karina/Simulation/MqSim/")

## Save files after
  #Int
#write.csv(live_size_df, row.names = FALSE, file=paste0("ShinyAppRuns/SIZE_int_MR_", MR_death_impact, "_restoration_", intercept_togg, "_numb_",intercept_indiv, "_level_", intercept_MR_mean ,".csv"))
#write.csv(MR_df, row.names = FALSE,file=paste0("ShinyAppRuns/MR_int_MR_", MR_death_impact, "_restoration_", intercept_togg, "_numb_",intercept_indiv, "_level_", intercept_MR_mean ,".csv"))

  #Base
#write.csv(live_size_df, file="Intervention/SIZE_base_0.4_bad_Intro1000.csv")
#write.csv(MR_df, file="Intervention/MR_base_0.4_bad_Intro1000.csv")

  #MR runs
#write.csv(live_size_df, file="SupplPlot_SupplementaryRuns/SIZE_int_MR_2.csv")
#write.csv(MR_df, file="SupplPlot_SupplementaryRuns/MR_int_MR_2.csv")
## Load in libraries
library(tidyverse)

## Load in parameters
#source("MainPlot_MR_tests/Rerun_Mult_iter/MR_ensrun_configurations_Low.txt")
source("Functions/mortality_functions_MRintro_hill.R")
source("Functions/mortality_functions_hill.R")
source("Functions/recruitment_functions_3.R")
source("Functions/disturbance_functions.R")

## Initiate random population
init_MR <- rnorm(n=population_size, mean = MR_mean, sd = MR_sd); init_MR[init_MR<0]=0; init_MR[init_MR>1]=1
init_age <- as.integer(runif(n=population_size, min=1, max=100))
indiv_ID <- seq(from=1, to=population_size)

pop_df <- list(
  indiv_ID = indiv_ID,
  time = rep(1, length(indiv_ID)),
  MR = init_MR,
  mortality = rep(0, length(indiv_ID)),
  age = init_age
)


## Initialise objects
indiv_pop=NULL 
MR_pop=NULL 
curr_pop=NULL 
new_recruit_pop=NULL 
indiv_count_end=0
death_df=NULL 
age_df=NULL
MR_df=NULL 
live_size_df=NULL
pop_timepoints=NULL
dist_event=FALSE
time_point=1
intercept_indiv_original = intercept_indiv

# Workflow: 
## 1) Population goes through disturbance
## 2A) Population recruits + intervention score (no mortality)
## 2B) Population mortality calculated
## 3) New population generated: Population initial -> Recruitment + intervention



for (time_point in 1:time_max){
  
  ##### Apply disturbance results
  if (dist_imp){
    disturbance_event <- disturbance_event_chance (dist_togg = dist_imp, disturbance_age_struct = disturbance_age_struct_type, dist_impact_val = dist_impact, dist_age_impact_val = dist_age_impact)
    disturbance_event_res <- disturbance_event[1]; age_impact = disturbance_event[2]; MR_death_impact = disturbance_event[3]; recruitment_const = disturbance_event[4]
  } else {
    disturbance_event_res=0
  }
  
  #### Initiating population
  if(time_point==1){
    curr_pop_start <- pop_df
  } else {
    curr_pop_start <- lapply(curr_pop_end, function(x) x[curr_pop_end$time==time_point])
  }
  
  if((length(curr_pop_start$indiv_ID)==0)) {
    stop("All dead at time ", time_point, "\n")
  } else{
    
    #### Count number of individuals
    indiv_count_start=length(curr_pop_start$indiv_ID) + indiv_count_end
    indiv_alive_count=length(curr_pop_start$indiv_ID)
    
    # Intercept inputs
    if (intercept_togg & time_point > intercept_timepoint & intercept_reducMort){intercept_pop_indiv_ID=intercept_pop$indiv_ID; int_togg=TRUE} else {int_togg=FALSE; intercept_pop_indiv_ID=NULL}
    
    if (time_point == intercept_timepoint+1 & intercept_togg){cat ("intercept_pop_indiv_ID = ", length(intercept_pop_indiv_ID), " ; int_togg = ",int_togg,"\n" )}
    
    #### Recruitment on intial start pop
    # If MR has not activated then regardless MR impact on rec is 0 ðŸ˜º
    if (time_point>=MR_timepoint & MR_lateintro & MR_imp){
      MR_recruit_impact_tp = MR_recruit_impact
    } else if (time_point<MR_timepoint & MR_lateintro & MR_imp) {
      MR_recruit_impact_tp = 0
    } else if (!MR_lateintro | !MR_imp) {
      MR_recruit_impact_tp = 0
    }
    
    curr_pop_recruited <- recruit_rate(pop=curr_pop_start, recruitment_age=recruitment_age, population_min_size=population_minimum_size, population_max_size=population_carrying_capacity, density_recruit_togg=density_recruit_toggle, recruitment_size_mean=recruitment_mean, recruitment_size_sd=recruitment_sd, recruitment_constant=recruitment_const, MR_togg=MR_rec_toggle, MR_recruit_impact_val=MR_recruit_impact_tp, MR_rec_adjusted=MR_rec_adj, age_togg=age_rec_toggle, age_recruit_impact_val=age_recruit_impact_value, rec_age_shiftch=rec_age_shift, MR_parents=MR_inherit_par_num)
    
    recruited_indivs = length(curr_pop_recruited$indiv_ID) - indiv_alive_count
    indiv_count_end = recruited_indivs + indiv_count_end
    
    ### !!Intervention on initial start pop!! âœŠ
    
    if (time_point == intercept_timepoint & intercept_togg){
      intercept_indiv = intercept_indiv_original
      if (intercept_indiv <= 0) {stop("too few intercept indivs")}
      
      cat("Time at:", time_point,"\n",
          "Individuals alive:", length(curr_pop_recruited$indiv_ID), "\n",
          "Mean MR of live individuals:", mean(curr_pop_recruited$MR), "\n")
      
      int_MR <- rnorm(n=intercept_indiv, mean = intercept_MR_mean, sd = intercept_MR_sd); int_MR[int_MR<0]=0; int_MR[int_MR>1]=1
      
      intercept_pop <- list(
        indiv_ID = seq(from = indiv_count_end + 1,
                       to   = indiv_count_end + intercept_indiv),
        time = rep(time_point, intercept_indiv),
        MR=as.numeric(int_MR),
        #mortality = rep(0, intercept_indiv),
        age = rep(2, intercept_indiv))
      
      curr_pop_og <- curr_pop_recruited
      curr_pop_int <- list( # Merging onto previous population recruited
        indiv_ID=c(curr_pop_recruited$indiv_ID, intercept_pop$indiv_ID), 
        age=c(curr_pop_recruited$age, intercept_pop$age), 
        MR=c(curr_pop_recruited$MR, intercept_pop$MR), 
        time=c(curr_pop_recruited$time, intercept_pop$time))
      
      curr_pop_recruited <- curr_pop_int 
      
      indiv_count_end=length(curr_pop_recruited$indiv_ID) + indiv_count_end # Adding onto indiv_ID count
      
      cat(  "## Intervention ##\n",
            "Individuals alive:", length(curr_pop_recruited$indiv_ID), "\n",
            "Mean MR of live individuals after:", mean(curr_pop_recruited$MR), "\n",
            "#################\n")
      
      # Plots
      cat(str(curr_pop_recruited))
      #print(ggplot() + geom_point(aes(x=curr_pop_recruited$age, y=curr_pop_recruited$MR)) + theme_bw() + labs(title=paste("MR by age at", time_point)))
    } else {intercept_indiv=0}
    
    #### Population mortality on initial population
    
    # Different ifelse for each scenario
    
    if (time_point>=MR_timepoint & MR_lateintro & MR_imp){
      indiv_death <- mortality_death_rate_MRlate(pop=curr_pop_start, population_capacity=population_carrying_capacity, comp_togg=comp_imp, comp_impact_val=comp_impact, MR_death_impact_val=MR_death_impact, MR_age_impact_val=MR_age_impact, age_impact_val=age_impact, mortality_age_shiftch=mortality_age_shift, MR_intro=MR_lateintro, MR_intro_timepoint=MR_timepoint, int_togg=int_togg, intercept_pop_indiv_ID=intercept_pop_indiv_ID)
      if(time_point==MR_timepoint){cat("Using mortality_death_rate_MRlate \n")}
      
    } else if (time_point<MR_timepoint & MR_lateintro & MR_imp){
      MR_death_impact_beforeintro = 0
      indiv_death <- mortality_death_rate(pop=curr_pop_start, population_capacity=population_carrying_capacity, comp_togg=comp_imp, comp_impact_val=comp_impact, MR_togg=MR_imp, MR_death_impact_val=MR_death_impact_beforeintro, MR_age_impact_val=MR_age_impact, age_impact_val=age_impact, mortality_age_shiftch=mortality_age_shift, int_togg=int_togg, intercept_pop_indiv_ID=intercept_pop_indiv_ID)
      if(time_point==1){cat("Using mortality_death_rate & MR before imp \n")}
      
    } else if (!MR_lateintro | !MR_imp) {
      indiv_death <- mortality_death_rate(pop=curr_pop_start, population_capacity=population_carrying_capacity, comp_togg=comp_imp, comp_impact_val=comp_impact, MR_togg=MR_imp, MR_death_impact_val=MR_death_impact, MR_age_impact_val=MR_age_impact, age_impact_val=age_impact, mortality_age_shiftch=mortality_age_shift, int_togg=int_togg, intercept_pop_indiv_ID=intercept_pop_indiv_ID)
      if(time_point==1){cat("Using mortality_death_rate \n")}
      
    }
    
    #### Tracking death counts
    death_df_curr <- data.frame(Dead_ID=curr_pop_start$indiv_ID[as.logical(indiv_death)], age=curr_pop_start$age[as.logical(indiv_death)], MR=curr_pop_start$MR[as.logical(indiv_death)], time=curr_pop_start$time[as.logical(indiv_death)])
    death_df <- rbind(death_df, death_df_curr)
    
    indiv_death = c(indiv_death, rep(0, recruited_indivs + intercept_indiv)) 
    
    #### Final time point pop
    curr_pop_end <- list(indiv_ID=curr_pop_recruited$indiv_ID[!as.logical(indiv_death)], age=curr_pop_recruited$age[!as.logical(indiv_death)]+1, MR=curr_pop_recruited$MR[!as.logical(indiv_death)], time=curr_pop_recruited$time[!as.logical(indiv_death)]+1)
    
    # Generate summary data if not dead
    if (length(curr_pop_recruited$indiv_ID)>sum(indiv_death)){ 
      MR_summ <- data.frame(time=time_point, MR_mean_summ=mean(curr_pop_end$MR, na.rm=TRUE), MR_sd_summ=sd(curr_pop_end$MR))
      MR_df <- rbind(MR_df, MR_summ)
      
      age_summ <- data.frame(time=time_point, age_mean_summ=mean(curr_pop_end$age), age_sd_summ=sd(curr_pop_end$age),pop_size=length(curr_pop_end$indiv_ID))
      age_df <- rbind(age_df, age_summ) 
    }
    
    live_size <- data.frame(time=time_point, sum_size=length(curr_pop_end$indiv_ID))
    live_size_df <- rbind(live_size_df, live_size) 
    
    # Save populations at user designated timepoints
    if (!is.null(timepoint_pop_grab) && (time_point %in% timepoint_pop_grab)){
      i <- which(time_point == timepoint_pop_grab)
      pop_timepoints[[i]] <- curr_pop_end
    }
    
    # Return to base
    if (as.logical(disturbance_event_res)){
      age_impact = disturbance_event[5]
      MR_death_impact = disturbance_event[6]
      recruitment_const = disturbance_event[7]
    }
    
    # Verbose 
    if(time_point%%output_timept == 0){ 
      
      cat("Time at:", time_point,"\n",
          "Individuals alive:", length(curr_pop_end$indiv_ID), "\n",
          "Mean MR of live individuals:", mean(curr_pop_end$MR), "\n",
          "#################\n")
      
      # Plots
      #print(ggplot() + geom_point(data=data.frame(curr_pop_end), aes(x=age, y=MR)) + theme_bw() + labs(title=paste("MR by age at", time_point)))
      
    }
  }
}
# 
# # Post-run plots
# 
# mean_MR_time_death <- death_df %>% group_by(time) %>% summarise(mean_MR=mean(MR), sd_MR=sd(MR))
# 
# plot_livesize <- ggplot() + 
#   geom_point(data=live_size_df, aes(x=time, y=sum_size)) + 
#   stat_smooth(data=live_size_df, aes(x=time, y = sum_size), linewidth = 0.85, linetype="dashed", colour="grey40", span=10) +
#   geom_hline(yintercept=population_carrying_capacity, linewidth = 0.75, linetype="dashed", colour="chocolate") +
#   geom_hline(yintercept=population_minimum_size, linewidth = 0.75, linetype="dashed", colour="chocolate") +
#   geom_vline(xintercept=MR_timepoint, linewidth = 0.75, linetype="dashed", colour="chocolate") +
#   ggforce::facet_zoom(xlim=c(1100,1200)) +
#   labs(title="Live population size") 
# if(intercept_togg){plot_livesize <- plot_livesize + geom_vline(xintercept=intercept_timepoint, linewidth = 0.75, linetype="dashed", colour="forestgreen")}
#   
# plot_liveage  <- ggplot() + 
#   geom_point(data=age_df, aes(x=time, y=age_mean_summ)) +
#   geom_errorbar(data=age_df, aes(x=time, ymax = age_mean_summ + age_sd_summ, ymin = age_mean_summ - age_sd_summ)) + 
#   labs(title="Live age") 
# 
# plot_deadMR   <- ggplot(mean_MR_time_death, aes(x=time, y=mean_MR)) + geom_point() + labs(title="Death MR")
# plot_liveMR   <- ggplot() +
#   geom_point(data=MR_df, aes(x=time, y = MR_mean_summ)) +
#   geom_errorbar(data=MR_df, aes(x=time, ymax = MR_mean_summ + MR_sd_summ, ymin = MR_mean_summ - MR_sd_summ)) + 
#   stat_smooth(data=MR_df, aes(x=time, y = MR_mean_summ), linewidth = 0.75, linetype="dashed", colour="grey40", span=10) +
#   geom_vline(xintercept=MR_timepoint, linewidth = 0.75, linetype="dashed", colour="chocolate") +
#   ggforce::facet_zoom(xlim=c(1100,1200)) +
#   labs(title="Live MR")
# if(intercept_togg){plot_liveMR <- plot_liveMR + geom_vline(xintercept=intercept_timepoint, linewidth = 0.75, linetype="dashed", colour="forestgreen")}
# 
# library(patchwork)
# plot_deadMR / plot_liveMR + plot_layout(heights = c(1,3))
# plot_livesize / plot_liveage + plot_layout(heights = c(3,1))

