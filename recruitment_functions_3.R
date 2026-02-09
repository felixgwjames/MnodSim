## Function - Recruitment - density dependent recruitment

## Inputs
# pop - Population
# recruitment_age - age to maturity
# recruitment_constant - chance of recruitment (1 = recruits, 0 = none)
# population_min_size - Increased fecundity due to disturbance

# MR_togg - Turn on and off the Myrtle rust affect
# MR_recruit_impact_val - impact of Myrtle rust susceptibility on recruitment (multiplier)

# Note: MR is recruited based on parent pheno using beta distribution. To visualise parent MR pheno of 1-4: for (i in 0:4){hist(rbeta(n=4000, shape1=i+1, shape2 = 3)) }

##############################################
recruit_rate <- function(pop, population_min_size, population_max_size, recruitment_age, recruitment_size_mean, density_recruit_togg, recruitment_size_sd, recruitment_constant, age_togg, age_recruit_impact_val, MR_togg, MR_recruit_impact_val, MR_rec_adjusted, rec_age_shiftch, MR_parents){
  
  # Current pop_size
  if (length(pop$indiv_ID) < population_min_size){
    recruitment_constant <- recruitment_constant*10
  } else if (length(pop$indiv_ID) >= population_min_size & length(pop$indiv_ID) < population_max_size & density_recruit_togg) {
    recruitment_adjust = 1-length(pop$indiv_ID)/population_max_size
  } else {recruitment_adjust=1}
  
  # How many individuals are fecund and recruit?
  ages <- pop$age
  fecund_indivs <- ages >= recruitment_age
  recruitment_indivs <- lapply(pop, function(x) x[fecund_indivs])
  
  #recruitment_indivs_MR <- rescale(recruitment_indivs$MR, to = c(0,1))
  recruitment_indivs_MR <- recruitment_indivs$MR
  recruitment_indivs_ages <- recruitment_indivs$age
  
  ### Reducing recruitment chances
  # Default impacts
  MR_impact  <- 1
  age_impact <- 1
  
  MR_impact <- 1 - recruitment_indivs_MR ^ MR_recruit_impact_val
  age_scaled <- recruitment_indivs_ages / rec_age_shiftch; age_scaled[age_scaled>1]=1
  age_impact <- age_scaled ^ age_recruit_impact_val
  
  # Both MR and age impacts
  if (MR_togg & age_togg) {
    total_impact <- MR_impact + age_impact - (MR_impact * age_impact) 
  } else if (MR_togg & !age_togg) {
    total_impact <- MR_impact
  } else if (age_togg & !MR_togg) {
    total_impact <- age_impact
  } else {
    total_impact <- 1
  }
  
  # Total impact on individuals recruitment
  total_impact <- total_impact*recruitment_adjust
  
  #tmp <- data.frame(recruitment_indivs_ages=recruitment_indivs_ages, recruitment_indivs_MR=recruitment_indivs_MR, MR_impact=MR_impact, age_impact=age_impact, total_impact=total_impact)
  #ggplot(data=tmp, aes(x=recruitment_indivs_MR, y=total_impact, colour=recruitment_indivs_ages)) + geom_point() +theme_bw()
  #ggplot(data=tmp, aes(x=recruitment_indivs_ages, y=total_impact, colour=recruitment_indivs_ages)) + geom_point() +theme_bw()
  #ggplot(data=tmp, aes(x=recruitment_indivs_MR, y=MR_impact, colour=recruitment_indivs_ages)) + geom_point() +theme_bw()
  #ggplot(data=tmp, aes(x=recruitment_indivs_ages, y=age_impact, colour=recruitment_indivs_ages)) + geom_point() +theme_bw()
  
  # Recruitment draw
  indiv_recruitment <- rbinom(n = length(recruitment_indivs$indiv_ID), size = 1, prob = recruitment_constant * total_impact)
  
  # For recruited individuals, what are their MR statuses
  
  if (sum(indiv_recruitment)>0){
    recruitment_indiv_MR <- recruitment_indivs_MR[as.logical(indiv_recruitment)]
    
    new_recruit <- as.integer(rnorm(n=sum(indiv_recruitment), mean = recruitment_size_mean, sd = recruitment_size_sd)); new_recruit[new_recruit<1]=1
    
    new_recruit_MR=NULL
    for (i in 1:length(recruitment_indiv_MR)) { # For each new recruit, use parent phenotype to generate MR, dependent on MR
      
      if(MR_parents == 2){
        parent_1_MR = recruitment_indiv_MR[i]
        parent_2_MR <- sample(recruitment_indiv_MR, size=1)
        AGV = mean(c(parent_1_MR, parent_2_MR)) # Additive genetic variance 
      } else {
        parent_1_MR = recruitment_indiv_MR[i]
        AGV = parent_1_MR # Additive genetic variance 
      }
      
      MR_rec_PDF <- rbeta(n=new_recruit[i], shape1=(AGV+1)^2, shape2=(2-AGV)^2) # distribution of recruited individual's MR between 0 to 1
      new_recruit_MR_new <- MR_rec_PDF+(AGV+MR_rec_adjusted-mean(MR_rec_PDF)) # recalibrate to make the mean the MR of parent pheno
      new_recruit_MR <- append(new_recruit_MR, new_recruit_MR_new)
    }
    
    new_recruit_MR[new_recruit_MR<0]=0; new_recruit_MR[new_recruit_MR>=1]=1
    
    new_recruit_pop <- list(indiv_ID=seq(from=indiv_count_start+1, to=indiv_count_start+sum(new_recruit)), 
                            time=rep(time_point, sum(new_recruit)), 
                            MR=new_recruit_MR, 
                            mortality=rep(0, sum(new_recruit)), 
                            age=rep(1, sum(new_recruit)))
    
    curr_pop <- list(
      indiv_ID=c(pop$indiv_ID, new_recruit_pop$indiv_ID), 
      age=c(pop$age, new_recruit_pop$age), 
      MR=c(pop$MR, new_recruit_pop$MR), 
      time=c(pop$time, new_recruit_pop$time), 
      mortality=c(pop$mortality, new_recruit_pop$mortality))
    
    return(curr_pop)
    
  } else {
    curr_pop <- list(
      indiv_ID=c(pop$indiv_ID, NULL), 
      age=c(pop$age, NULL), 
      MR=c(pop$MR, NULL), 
      time=c(pop$time, NULL), 
      mortality=c(pop$mortality, NULL))
  }
  return(curr_pop)
  
}
