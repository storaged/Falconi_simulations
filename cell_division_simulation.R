source("library.R")
source("Queue.R")


### CONSTANT VALUES
GENOME_SIZE <- 3.055 * 10^9
EXOME_SIZE  <- 40 * 10^6  

run_simulation <- function(par_list, 
                           #NO_GENERATIONS, NO_GENOME_MUT_PER_GEN = 400, 
                           #ONCOGENES_PERCENT = 0.01, 
                           debug = F, info = T){
  
  source("Queue.R")
  source("library.R")
  
  generate_child_cell <- function(cell, current_generation, number_of_exom_mutations,
                                  cell_ID){
    return(list(curr_gen = current_generation, 
                remaining_evolution_events = 2,
                next_generation_evolutionary_events = 0,
                mut_moments = ifelse(is.na(cell$mut_moments), 
                                     cell$curr_gen, 
                                     c(cell$mut_moments, cell$curr_gen)), 
                total_mut = cell$total_mut + number_of_exom_mutations,
                parent_ID = cell$parent_ID,
                cell_ID = cell_ID))
  }
  
  NO_GENOME_MUT_PER_GEN <- par_list$NO_GENOME_MUT_PER_GEN
  ONCO_MUTATION_PROB <- par_list$ONCOGENES_PERCENT * par_list$EXOME_SIZE / par_list$GENOME_SIZE
  NO_GENERATIONS <- par_list$NO_GENERATIONS

  current_generation <- 1
  
  history <- Queue$new()
  next_generation <- Queue$new()
  
  #number_of_exom_mutations <- rbinom(1, NO_GENOME_MUT_PER_GEN, ONCO_MUTATION_PROB)
  #mut_moments <- c()
  #if(number_of_exom_mutations) { mut_moments <- 1 }
  history$push(list(curr_gen = 1, 
                    remaining_evolution_events = 2,
                    next_generation_evolutionary_events = 0,
                    mut_moments = NA, #ifelse(number_of_exom_mutations, 1, 0), 
                    total_mut = 0, #number_of_exom_mutations,
                    parent_ID = 0,
                    cell_ID = 0))
  history$data %>% length
  
  iter <- 1
  cell_ID <- 1
  while (current_generation < NO_GENERATIONS) {
    if (info) cat(paste0("Generation number: ", current_generation))
    if (debug) print(paste0("generations: peek$curr_gen=", history$peek()$curr_gen, 
                 " : current_generation=",current_generation, 
                 " : NO_GENERATIONS=", NO_GENERATIONS))
    if (debug) print(paste0("remaining evololution events = ", history$peek()$remaining_evolution_events))
    while (!is.null(history$peek())){
    #print(history$peek()$curr_gen)
    #print(current_generation)
      current_cell <- history$poll()
      
      if (debug) print(paste0("(it:",iter,")\n",
                   ", current_cell: [$curr_gen=",current_cell$curr_gen,
                   ", $next_generation_evolutionary_events=", current_cell$next_generation_evolutionary_events, 
                   ", $remaining_evolution_events=", current_cell$remaining_evolution_events,"]"))
      while (current_cell$remaining_evolution_events > 0){
        
        number_of_exom_mutations <- rbinom(1, NO_GENOME_MUT_PER_GEN, ONCO_MUTATION_PROB)
        
        if (number_of_exom_mutations) { 
            cellA <- generate_child_cell(current_cell, current_generation, 
                                         number_of_exom_mutations, cell_ID)
            next_generation$push(cellA)
            cell_ID <- cell_ID + 1
        }
        else {
            current_cell$next_generation_evolutionary_events <- current_cell$next_generation_evolutionary_events + 1
        }
        
        current_cell$remaining_evolution_events <- current_cell$remaining_evolution_events - 1
      }
      current_cell$remaining_evolution_events <- 2 * current_cell$next_generation_evolutionary_events
      current_cell$next_generation_evolutionary_events <- 0  
      next_generation$push(current_cell)
    }
    current_generation <- current_generation + 1
    history <- next_generation
    next_generation <- Queue$new()
    iter <- iter + 1
    if (debug) if (iter > 20) { print("Breaking."); break }
  }
  results <- do.call(rbind, history$data)
  results$NO_GENERATIONS <- par_list$NO_GENERATIONS
  results$ONCOGENES_PERCENT <- par_list$ONCOGENES_PERCENT
  results$NO_GENOME_MUT_PER_GEN <- par_list$NO_GENOME_MUT_PER_GEN
  results$RUN_NUMBER <- par_list$iterations
  results
}

cl <- makeCluster(20)

iterations <- 1:100
oncogene_percent <- seq(0.002, 0.02, by = 0.005)
number_of_generations <- 10
mutations_per_cell <- seq(100, 400, by = 100)

parameters_grid <- expand.grid(NO_GENERATIONS = number_of_generations, 
            ONCOGENES_PERCENT = oncogene_percent, 
            NO_GENOME_MUT_PER_GEN = mutations_per_cell,
            EXOME_SIZE = EXOME_SIZE,
            GENOME_SIZE = GENOME_SIZE,
            iterations = iterations)

parameters_grid$ID <- 1:nrow(parameters_grid)

parametrs_list <- lapply(split(parameters_grid,parameters_grid$ID), as.list)


results_parallel <- parLapply(cl, parametrs_list, run_simulation)

# 35 is the number of divisions of a newborn



sampling_number_of_mutations <- lapply(c(100,150,200,250,300), function(mut_per_division){
  bootstrapped <- sapply(1:100, function(i){
    
    if (i %% 50 == 0 ) print(paste0("I: ", i, "Mut.per.div: ", mut_per_division))
    run_simulation(10, mut_per_division, debug = F, info = F) -> simulation_result
    res_tab <- sapply(simulation_result$data, function(cell) 
      c(cell$curr_gen, cell$total_mut, 
        cell$remaining_evolution_events, 
        cell$parent_ID, cell$cell_ID)
    ) 
    
    rownames(res_tab) <- c("Current Generation", "Total mutation", 
                           "Remaining evolution events", 
                           "parent ID", "cell_ID")
  
    c(summary(res_tab["Current Generation", res_tab["Total mutation",] != 0]), mut_per_division=mut_per_division)
  })
  
  data.frame(t(bootstrapped))
  
}) %>% do.call(rbind, .)

ggplot(sampling_res, aes(x =  mut_per_division, y = `Min.`, group = mut_per_division)) + 
  geom_boxplot() + 
  theme_minimal() + 
  xlab("Expected number of mutation events in a cell") + 
  ylab("Expected number of divisions until the first mutation event")

cat(paste0("Divisions until the first mutation: median = ", median(bootstrapped["Min.",]),", 95% Confidence interval: [", quantile(bootstrapped["Min.",], 0.025), 
           ", ", quantile(bootstrapped["Min.",], 0.975), "]"))



