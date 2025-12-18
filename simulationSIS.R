library(igraph)

# Parameters
n <- 1000             # Network size
p0 <- 0.05            # Initial infected fraction
beta <- 0.05          # Probability of infection
gamma <- 0.40         # Probability of recovery
time_steps <- 100     # Duration of simulation

# Network Generation
networks <- list(
  ER = erdos.renyi.game(n, p = 0.02),
  BA = barabasi.game(n, m = 5, directed = FALSE),
  WS = watts.strogatz.game(dim = 1, size = n, nei = 5, p = 0.05),
  Complete = make_full_graph(n),
  Star = make_star(n, mode = "undirected")
)

# SIS Simulation Function
simulate_sis <- function(graph, p0, beta, gamma, steps) {
  n_nodes <- vcount(graph)
  status <- rbinom(n_nodes, 1, p0) 
  infected_history <- numeric(steps)
  adj <- as_adj(graph) 
  
  for (t in 1:steps) {
    infected_history[t] <- sum(status) / n_nodes
    new_status <- status
    
    for (i in 1:n_nodes) {
      if (status[i] == 1) {
        # Recovery with probability gamma
        if (runif(1) < gamma) new_status[i] <- 0
      } else {
        # Infection attempt from neighbors with probability beta
        neighbors <- which(adj[i, ] == 1)
        if (length(neighbors) > 0) {
          inf_neighs <- sum(status[neighbors])
          if (runif(1) < (1 - (1 - beta)^inf_neighs)) {
            new_status[i] <- 1
          }
        }
      }
    }
    status <- new_status
  }
  return(infected_history)
}

# Run Simulations
results <- lapply(networks, simulate_sis, p0, beta, gamma, time_steps)

# Plotting results
colors <- c("blue", "red", "green", "purple", "orange")
plot(results$ER, type = "l", col = colors[1], ylim = c(0, 1), 
     xlab = "Time Step", ylab = "Proportion of Infected Nodes")
lines(results$BA, col = colors[2])
lines(results$WS, col = colors[3])
lines(results$Complete, col = colors[4])
lines(results$Star, col = colors[5])

legend("topright", legend = names(results), col = colors, lty = 1, cex = 0.8)

# Analysis of Eigenvalues and Thresholds
par(mfrow=c(2, 3)) 
results_task2 <- list()

for (name in names(networks)) {
  # Calculate lambda1 and Forecast threshold is 1/lambda1
  adj_matrix <- as_adj(networks[[name]])
  lambda1 <- eigen(adj_matrix, only.values = TRUE)$values[1]
  lambda1 <- Re(lambda1)
  
  # Determine Critical Beta
  beta_critical <- gamma / lambda1
  beta_low  <- beta_critical * 0.9
  beta_high <- beta_critical * 1.1
  if(beta_high > 1) beta_high <- 1
  
  cat(name, "Network:\n")
  cat("  Leading Eigenvalue (lambda1):", round(lambda1, 3), "\n")
  cat("  Forecast Threshold (1/lambda1):", round(1/lambda1, 4), "\n")
  cat("  Epidemic likely? ", ifelse(lambda1 > (gamma/beta), "YES", "NO"), "\n\n")
  cat(sprintf("  Critical Beta: %.4f\n", beta_critical))
  cat(sprintf("  Running Low Beta: %.4f | High Beta: %.4f\n", beta_low, beta_high))
  
  # Simulate
  sim_low  <- simulate_sis(g, p0, beta_low, gamma, time_steps)
  sim_high <- simulate_sis(g, p0, beta_high, gamma, time_steps)
  
  # Plot Comparison
  plot(sim_high, type="l", col="red", ylim=c(0,1), lwd=2,
       main=paste(name, "\nThresh Beta:", round(beta_critical, 4)),
       ylab="Infected Ratio", xlab="Time")
  lines(sim_low, col="blue", lwd=2)
  abline(h=0, lty=2, col="gray")
  if(name == "ER") legend("topright", c("Above Thresh", "Below Thresh"), col=c("red", "blue"), lty=1, cex=0.7)
}
# Reset layout
par(mfrow=c(1,1))