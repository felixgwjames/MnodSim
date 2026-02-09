## Function - Disturbance function impacts

## Inputs
# pop - Population
# dist_togg - Turn on and off disturbance presence
# disturbance_age_struct - Impact value flat or varied by age
# dist_prob - Base chance of a disturbance event
# dist_impact_val - Level of impact 
# dist_age_impact_val - Impact of age on disturbance (inflection point of exponential function)

##################################################################################
disturbance_event_chance <- function(dist_togg, disturbance_age_struct, dist_prob, dist_impact_val, dist_age_impact_val, pop=curr_pop, age_imp=age_impact, MR_death_imp=MR_death_impact){
  
  
  dist_event=FALSE
  
  if (dist_togg){
    
    # Does a disturbance occur?
    disturbance_size_chance <- rbeta(n=1, 5, 1) # exponential chance of size from 0 -1 impact
    disturbance_size <- 1-disturbance_size_chance
    dist_event <- rbinom(1,1,disturbance_size_chance) # lower disturbances have a higher chance
    
    
    if (as.logical(dist_event) & disturbance_age_struct=="complex"){
      #cat("Disturbance size of", round(disturbance_size,2),"at:", time_point,"!🔥🔥🔥 \n")
      
      ages <- pop$age
      disturbance_age_imp  <- disturbance_size * dist_impact_val * (1/(1+exp((ages-dist_age_impact_val)/10)))
      
      age_impact_new <- age_imp * (1+disturbance_age_imp)
      MR_death_impact_new <- MR_death_imp * (1+disturbance_age_imp)
      recruitment_const_new <- recruitment_const * (1+disturbance_age_imp)
    } 
    
    if (as.logical(dist_event) & disturbance_age_struct=="flat"){
      age_impact_new <- age_imp * (1+dist_impact_val)
      MR_death_impact_new <- MR_death_imp * (1+dist_impact_val)
      recruitment_const_new <- recruitment_const * (1+dist_impact_val)
    }
    
    else {
      age_impact_new=age_imp
      MR_death_impact_new=MR_death_imp
      recruitment_const_new=recruitment_const
    }
    
    results <- as.vector(c(dist_event, age_impact_new, MR_death_impact_new, recruitment_const_new, age_imp, MR_death_imp, recruitment_const))
    
    return(results)
    
  } else {
    return()
  }
}
