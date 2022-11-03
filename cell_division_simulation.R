generate_child_cell <- function(cell, NO_EXOME_MUT_PER_GEN, EXOM_MUTATION_PROB){
  number_of_exom_mutations <- rbinom(1, NO_EXOME_MUT_PER_GEN, EXOM_MUTATION_PROB)
  
  if (number_of_exom_mutations) { 
    return(list(curr_gen = cell$curr_gen + 1, 
                remaining_evolution_events = 0,
                next_generation_evolutionary_events = 1,
                mut_moments = c(cell$mut_moments, cell$curr_gen), 
                total_mut = cell$total_mut + number_of_exom_mutations))
  } else {
    return(NULL)
  }
}

run_simulation <- function(NO_GENERATIONS){

  GENOME_SIZE <- 6.27 * 10^9
  EXOME_SIZE  <- 10 * 10^6
  
  NO_GENOME_MUT_PER_GEN <- 400
  NO_EXOME_MUT_PER_GEN  <- 4
  
  EXOM_MUTATION_PROB <- 0.0 #NO_EXOME_MUT_PER_GEN / EXOME_SIZE
  
  NO_GENERATIONS <- NO_GENERATIONS

  current_generation <- 1
  
  history <- Queue$new()
  
  number_of_exom_mutations <- rbinom(1, NO_EXOME_MUT_PER_GEN, EXOM_MUTATION_PROB)
  mut_moments <- c()
  if(number_of_exom_mutations) { mut_moments <- 1 }
  history$push(list(curr_gen = 1, 
                    remaining_evolution_events = 1,
                    next_generation_evolutionary_events = 0,
                    mut_moments = ifelse(number_of_exom_mutations, 1, 0), 
                    total_mut = number_of_exom_mutations))
  history$data %>% length
  iter <- 1
  while (current_generation < NO_GENERATIONS) {
    print(paste0(current_generation, ":",NO_GENERATIONS))
    #print(history$peek()$curr_gen)
    #print(current_generation)
    if( history$peek()$curr_gen == current_generation & cell$remaining_evolution_events == 0) {
       print(paste0("IF CLAUSE (it:",iter,"), CELL: cell$curr_gen=",cell$curr_gen,", current_generation=",current_generation, ", cell$rem_events=", cell$remaining_evolution_events,"]"))
       current_generation <- current_generation + 1
    } else {
      print(paste0("CELL: [",cell$curr_gen,", ",cell$remaining_evolution_events,",","]"))
      cell <- history$poll()
      #print(cell)
      unchanged <- 0
      for(i in 1:2) {
        cellA <- generate_child_cell(cell, NO_EXOME_MUT_PER_GEN, EXOM_MUTATION_PROB )
        if(!is.null(cellA)) {
          history$push(cellA)
        }
        else {
          unchanged <- unchanged + 1
        }
      }
      if (cell$remaining_evolution_events > 1) {
        cell$next_generation_evolutionary_events <- cell$next_generation_evolutionary_events + unchanged
        cell$remaining_evolution_events <- cell$remaining_evolution_events - 1
        history$push(cell)
      } else {
        cell$curr_gen <- cell$curr_gen + 1
        cell$remaining_evolution_events <- cell$next_generation_evolutionary_events
        cell$next_generation_evolutionary_events <- 0
        history$push(cell)
      }
    }
    #iter <- iter + 1
    #if(iter > 20) break
  }

  history
}

run_simulation(10) -> simulation_result
sapply(simulation_result$data, function(cell) 
  c(cell$curr_gen, cell$total_mut, cell$remaining_evolution_events)
) 
