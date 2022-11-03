source("library.R")
source("Queue.R")


### CONSTANT VALUES
GENOME_SIZE <- 6.27 * 10^9
EXOME_SIZE  <- 10 * 10^6

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
    if( history$peek()$curr_gen == current_generation & history$peek()$remaining_evolution_events == 0) {
       print(paste0("IF CLAUSE (it:",iter,")",
                    ", CELL: cell$curr_gen=",history$peek()$curr_gen,
                    ", current_generation=",history$peek()$current_generation, 
                    ", cell$rem_events=", history$peek()$remaining_evolution_events,"]"))
       current_generation <- current_generation + 1
    } else {
      
      cell <- history$poll()
      print(paste0("CELL: [",cell$curr_gen,", ",cell$remaining_evolution_events,",","]"))
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
    iter <- iter + 1
    if (iter > 20) break
  }

  history
}

run_simulation(10) -> simulation_result
sapply(simulation_result$data, function(cell) 
  c(cell$curr_gen, cell$total_mut, cell$remaining_evolution_events)
) 

get_error_var <- function(size = 1000*2^4, sd = .015) {rnorm(size, 0, sd)}
did_happened <- function(size = 1000*2^4) {sample(c(0,1), size, replace = T)}
first <- 0.25 + get_error_var()
second <- first * 1/(did_happened()+1) + get_error_var()
third <- second * 1/(did_happened()+1) + get_error_var()

hist(third, breaks = 100)

my_error_var <- 0.025
first <-  0.25 + get_error_var(sd = my_error_var)
second <- 0.125 + get_error_var(sd = my_error_var)
third <-  0.0625 + get_error_var(sd = my_error_var)
fourth <-  0.03125 + get_error_var(sd = my_error_var)
divisions_data <- rbind(
  data.frame(division = "first", percent = first),
  data.frame(division = "second", percent =  second),
  data.frame(division = "third", percent = third),
  data.frame(division = "fourth", percent = fourth)
)

stack <- ggplot(divisions_data, aes(x = percent)) + 
  geom_density(position = "stack", alpha = .25) + theme_minimal()
dodge <- ggplot(divisions_data, aes(x = percent, fill = division)) + 
  geom_density(position = "dodge", alpha = .25) + 
  geom_vline(xintercept = c(0.25,0.125, 0.0625, 0.03125), linetype = "dotted") + 
  theme_minimal() +
  theme(legend.position = "bottom") 

grid.arrange(stack, dodge, nrow=2, heights=c(4, 5))

