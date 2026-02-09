## Function - Mortality death rate _ Hill Function

######################################################################################

Hill_eqn <- function(age, steepness, mortality_age_shiftch, flip) {
  K <- mortality_age_shiftch / (((steepness - 1) / (steepness + 1))^(1 / steepness))
  y <- if(flip){0.75*K^steepness / (K^steepness + age^steepness)} else {age^steepness / (K^steepness + age^steepness)}
  return(y)
}

young_mortality <- function(age_x, age_impact_val, mortality_age_shiftch){
  # Mortality chance for individuals <age_impact_val
  death_perc <- Hill_eqn(age=age_x, mortality_age_shiftch=2, steepness = 2, flip=TRUE) * age_impact_val # Decreasing death with age
  return(death_perc)
}


mature_mortality <- function(age_x, age_impact_val, mortality_age_shiftch){
  # Mortality chance for individuals >=age_impact_val
  death_perc <- Hill_eqn(age=age_x, mortality_age_shiftch=mortality_age_shiftch, steepness = 15, flip=FALSE) * age_impact_val# Rising death with age
  return(death_perc)
}

######################################################################################
mortality_death_rate  <- function(pop, population_capacity, population_min_size, comp_togg, comp_impact_val, MR_togg, MR_death_impact_val, MR_age_impact_val, age_impact_val, mortality_age_shiftch, intercept_pop_indiv_ID, int_togg){
  
  require(scales)
  # Age death
  ages <- pop$age
  
  if(int_togg){
    #cat("int_togg", int_togg, "\n")
    match_idx <- pop$indiv_ID %in% intercept_pop_indiv_ID
    #cat("!!!", unique(match_idx), "\n")
    ages[match_idx] <- ifelse(ages[match_idx] < mortality_age_shiftch/4, ages[match_idx]*4, ages[match_idx])
  }
  
  MR <- pop$MR
  age_mortality_chance <- numeric(length(ages))
  
  for (i in seq_along(ages)) {
    x <- ages[i]
    if (x <= mortality_age_shiftch/2) {
      age_mortality_chance[i] <- young_mortality(x, age_impact_val, mortality_age_shiftch)
    } else {
      age_mortality_chance[i] <- mature_mortality(x, age_impact_val, mortality_age_shiftch)
    }
  }
  
  age_mortality_chance[age_mortality_chance > 1] <- 1 # remove anomaly high and low vals
  age_mortality_chance[age_mortality_chance < 0] <- 0 # remove anomaly high and low vals
  
  ## Plot age
  #ggplot() + geom_point(aes(x=ages, y=age_mortality_chance))
  
  # MR chance by death
  #MR <- rescale(pop$MR, c(0,1))
  
  if(!(MR_age_impact_val==0)){
    MR_chance <- (1 / (1 + (ages / MR_age_impact_val))) * MR * MR_death_impact_val
  } else {
    MR_chance <- MR * MR_death_impact_val
  }
  
  
  MR_chance[MR_chance > 1] <- 1; MR_chance[MR_chance < 0] <- 0 # remove anomaly high and low vals
  
  ## Plot MR
  # ggplot() + geom_point(aes(x=MR, y=MR_chance,colour=ages))
  
  # If both toggle is off 
  if (!MR_togg & !comp_togg)  { # If MR toggle is on but not competition
    final_mortality_chance_norm <- age_mortality_chance
  }
  
  # If MR toggle is on but not competition
  if (MR_togg & !comp_togg)  { # If MR toggle is on but not competition
    final_mortality_chance_norm <- MR_chance + age_mortality_chance - (MR_chance*age_mortality_chance) # conditional probability
    #ggplot() + geom_point(aes(x=MR_chance, y=age_mortality_chance, colour=final_mortality_chance_norm))
  }
  
  # If compeititon toggle is on but not MR
  if (!MR_togg & comp_togg) {
    pop_size = length(pop$indiv_ID)
    
    if(pop_size > population_capacity){
      comp_chance = (pop_size - population_capacity)/pop_size  * (1+comp_impact_val)  # Scale competition impact by how much over carrying capacity of population size
      final_mortality_chance_norm <- rescale(age_mortality_chance, from=c(0,1), to=c(comp_chance,1))
    } else {
      comp_chance=0
      final_mortality_chance_norm <- age_mortality_chance 
    }
  }
  
  # If both MR and compeititon toggle is on
  if (comp_togg & MR_togg) {
    pop_size = length(pop$indiv_ID)
    
    if(pop_size > population_capacity){
      comp_chance = (pop_size - population_capacity)/population_capacity  * (1+comp_impact_val)  # Scale competition impact by how much over carrying capacity of population size
      final_mortality_chance_norm <- rescale((MR_chance + age_mortality_chance - (MR_chance*age_mortality_chance)), from=c(0,1), to=c(comp_chance,1))
    } else {
      comp_chance=0
      final_mortality_chance_norm <- (MR_chance + age_mortality_chance - (MR_chance*age_mortality_chance)) 
    }
    
  }
  
  # Mortality
  mortality_base <- rbinom(n = length(final_mortality_chance_norm), size = 1, prob = final_mortality_chance_norm)
  return(mortality_base)
}

