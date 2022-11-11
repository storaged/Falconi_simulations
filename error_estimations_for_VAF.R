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