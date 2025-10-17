#NEW
# Install and load required packages
if (!require("igraph")) install.packages("igraph")
if (!require("parallel")) install.packages("parallel")
library(igraph)
library(parallel)

# Setting random seed for reproducibility
set.seed(42)

# Helper function for truncated beta sampling with error checking
truncated_beta <- function(alpha, beta, lower, upper) {
  tryCatch({
    if(lower >= upper) {
      warning("Lower bound >= upper bound in truncated_beta")
      return(lower)
    }
    x <- rbeta(1, alpha, beta)
    while(x < lower || x > upper) {
      x <- rbeta(1, alpha, beta)
    }
    return(x)
  }, error = function(e) {
    warning(paste("Error in truncated_beta:", e$message))
    return(lower)
  })
}

# Helper function for Dirichlet sampling
rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  x <- matrix(rgamma(n*k, alpha), ncol = k, byrow = TRUE)
  sm <- rowSums(x)
  return(x/sm)
}

# Generate synthetic network with debugging
generate_network <- function(n = 100, k = 3, theta = NULL) {
  cat(sprintf("[Debug] Generating network: n=%d, k=%d\n", n, k))
  
  if (is.null(theta)) {
    theta <- rep(1/k, k)
  }
  
  # Generate community assignments
  z <- sample(1:k, n, replace = TRUE, prob = theta)
  
  # Initialize P with assortative structure
  P <- matrix(c(
    0.30, 0.08, 0.08,
    0.08, 0.10, 0.02,
    0.08, 0.02, 0.10
  ), nrow = 3, byrow = TRUE)
  
  
  # Generate adjacency matrix
  A <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      A[i,j] <- A[j,i] <- rbinom(1, 1, P[z[i], z[j]])
    }
  }
  
  cat("[Debug] Network generation complete\n")
  return(list(A = A, z_true = z, P = P))
}
# Modified SBM Gibbs sampler with Algorithm 2 implementation
sbm_gibbs <- function(A, k, B = 2000, burnin = 500) {
  n <- nrow(A)
  cat(sprintf("[Debug] Starting Gibbs sampler: n=%d, k=%d, B=%d, burnin=%d\n", n, k, B, burnin))
  
  # Initialize parameters
  z <- sample(1:k, n, replace = TRUE)
  theta <- rep(1/k, k)
  P <- matrix(0.1, k, k)
  diag(P) <- 0.3
  tau <- 0.2  # Initial value for tau
  
  # Storage for samples
  z_samples <- matrix(0, nrow = B, ncol = n)
  P_samples <- array(0, dim = c(B, k, k))
  theta_samples <- matrix(0, nrow = B, ncol = k)
  tau_samples <- numeric(B)
  
  # Prior parameters
  alpha_theta <- rep(1, k)
  alpha_P <- beta_P <- matrix(1, k, k)
  
  # Helper function to sample tau using rejection sampling
  sample_tau <- function(P) {
    p <- min(diag(P))  # minimum of within-community probabilities
    q <- max(P[row(P) != col(P)])  # maximum of between-community probabilities
    
    if (p <= q) {
      warning("Assortativity constraint violated: p <= q")
      return(tau)  # Return current tau if constraint is violated
    }
    
    # Rejection sampling for tau
    accept <- FALSE
    max_attempts <- 1000
    attempt <- 0
    
    while (!accept && attempt < max_attempts) {
      attempt <- attempt + 1
      # Propose tau from uniform(q, p)
      tau_prop <- runif(1, q, p)
      
      # Calculate acceptance probability
      log_ratio <- -(k*(k-1)/2) * log(tau_prop) - k * log(1 - tau_prop)
      
      if (log(runif(1)) < log_ratio) {
        accept <- TRUE
        return(tau_prop)
      }
    }
    
    if (!accept) {
      warning("Maximum rejection sampling attempts reached for tau")
      return(tau)  # Return current tau if max attempts reached
    }
  }
  
  # Main Gibbs sampling loop
  for (iter in 1:(B + burnin)) {
    if (iter %% 500 == 0) {
      cat(sprintf("[Debug] Iteration %d/%d\n", iter, B + burnin))
    }
    
    # Update z
    for (i in 1:n) {
      probs <- rep(0, k)
      for (a in 1:k) {
        log_prob <- log(theta[a])
        for (j in 1:n) {
          if (j != i) {
            if (A[i,j] == 1) {
              log_prob <- log_prob + log(P[a,z[j]])
            } else {
              log_prob <- log_prob + log(1 - P[a,z[j]])
            }
          }
        }
        probs[a] <- exp(log_prob)
      }
      if (sum(probs) == 0) {
        warning("All probabilities zero in z update")
        probs <- rep(1/k, k)
      }
      z[i] <- sample(1:k, 1, prob = probs)
    }
    
    # Update theta
    n_a <- tabulate(z, k)
    theta <- rdirichlet(1, alpha_theta + n_a)
    
    # Update P with assortativity constraints
    for (a in 1:k) {
      for (b in a:k) {
        n_ab <- sum(A[z == a, z == b])
        m_ab <- sum(z == a) * sum(z == b)
        if (a == b) {
          m_ab <- m_ab - sum(z == a)
          # Within-community: truncated to (tau, 1)
          P[a,b] <- P[b,a] <- truncated_beta(
            alpha_P[a,b] + n_ab,
            beta_P[a,b] + m_ab - n_ab,
            tau, 1.0
          )
        } else {
          # Between-community: truncated to (0, tau)
          P[a,b] <- P[b,a] <- truncated_beta(
            alpha_P[a,b] + n_ab,
            beta_P[a,b] + m_ab - n_ab,
            0.001, tau
          )
        }
      }
    }
    
    # Update tau using the conditional distribution
    tau <- sample_tau(P)
    
    # Store samples after burnin
    if (iter > burnin) {
      z_samples[iter - burnin,] <- z
      P_samples[iter - burnin,,] <- P
      theta_samples[iter - burnin,] <- theta
      tau_samples[iter - burnin] <- tau
    }
  }
  
  cat("[Debug] Gibbs sampling complete\n")
  
  # Custom MAP estimation
  z_map <- rep(0, n)
  for (i in 1:n) {
    tab <- tabulate(z_samples[,i], k)
    z_map[i] <- which.max(tab)
  }
  
  return(list(
    z_map = z_map,
    z_samples = z_samples,
    P_samples = P_samples,
    theta_samples = theta_samples,
    tau_samples = tau_samples,
    P_final = P,
    theta_final = theta,
    tau_final = tau
  ))
}
# Function to calculate Adjusted Rand Index
adjustedRandIndex <- function(x, y) {
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) 
    stop("vectors must be same length")
  
  n <- length(x)
  t1 <- table(x)
  t2 <- table(y)
  t3 <- table(x, y)
  
  # Calculate the rand index components
  a <- sum(choose(t3, 2))
  b1 <- sum(choose(t1, 2))
  b2 <- sum(choose(t2, 2))
  d <- choose(n, 2)
  
  # Calculate expected index
  exp_ri <- (b1 * b2)/(d)
  
  # Calculate max index
  max_ri <- (b1 + b2)/2
  
  # Calculate adjusted rand index
  ari <- (a - exp_ri)/(max_ri - exp_ri)
  
  return(ari)
}
# Standard SBM Gibbs sampler (no assortativity)
sbm_gibbs_standard <- function(A, k, B = 2000, burnin = 500) {
  n <- nrow(A)
  cat(sprintf("[Debug] Starting standard SBM Gibbs sampler: n=%d, k=%d\n", n, k))
  
  # Initialize parameters
  z <- sample(1:k, n, replace = TRUE)
  theta <- rep(1/k, k)
  P <- matrix(0.1, k, k)
  diag(P) <- 0.3
  
  # Storage for samples
  z_samples <- matrix(0, nrow = B, ncol = n)
  P_samples <- array(0, dim = c(B, k, k))
  theta_samples <- matrix(0, nrow = B, ncol = k)
  
  # Prior parameters
  alpha_theta <- rep(1, k)
  alpha_P <- beta_P <- matrix(1, k, k)
  
  # Main Gibbs loop
  for (iter in 1:(B + burnin)) {
    if (iter %% 500 == 0) {
      cat(sprintf("[Debug] Iteration %d/%d (standard)\n", iter, B + burnin))
    }
    
    # Update z
    for (i in 1:n) {
      probs <- rep(0, k)
      for (a in 1:k) {
        log_prob <- log(theta[a])
        for (j in 1:n) {
          if (j != i) {
            if (A[i,j] == 1) {
              log_prob <- log_prob + log(P[a,z[j]])
            } else {
              log_prob <- log_prob + log(1 - P[a,z[j]])
            }
          }
        }
        probs[a] <- exp(log_prob)
      }
      if (sum(probs) == 0) {
        warning("All probabilities zero in standard z update")
        probs <- rep(1/k, k)
      }
      z[i] <- sample(1:k, 1, prob = probs)
    }
    
    # Update theta
    n_a <- tabulate(z, k)
    theta <- rdirichlet(1, alpha_theta + n_a)
    
    # Update P
    for (a in 1:k) {
      for (b in a:k) {
        n_ab <- sum(A[z == a, z == b])
        m_ab <- sum(z == a) * sum(z == b)
        if (a == b) m_ab <- m_ab - sum(z == a)
        
        P[a,b] <- P[b,a] <- rbeta(1,
                                  alpha_P[a,b] + n_ab,
                                  beta_P[a,b] + m_ab - n_ab)
      }
    }
    
    # Store after burnin
    if (iter > burnin) {
      z_samples[iter - burnin,] <- z
      P_samples[iter - burnin,,] <- P
      theta_samples[iter - burnin,] <- theta
    }
  }
  
  # MAP estimate
  z_map <- rep(0, n)
  for (i in 1:n) {
    tab <- tabulate(z_samples[,i], k)
    z_map[i] <- which.max(tab)
  }
  
  cat("[Debug] Standard SBM Gibbs sampling complete\n")
  return(list(
    z_map = z_map,
    z_samples = z_samples,
    P_samples = P_samples,
    theta_samples = theta_samples,
    P_final = P,
    theta_final = theta
  ))
}

# Function to run a single simulation with detailed error reporting
run_simulation <-function(k, use_assortative = TRUE) {
  tryCatch({
    cat(sprintf("[Debug] Starting simulation for k=%d\n", k))
    
    # Generate synthetic network
    net_data <- generate_network(n = 100, k = k)
    A <- net_data$A
    z_true <- net_data$z_true
    
    # Run Gibbs sampler
    cat(sprintf("[Debug] Running Gibbs sampler for k=%d\n", k))
    result <- if (use_assortative) sbm_gibbs(A, k) else sbm_gibbs_standard(A, k)
    
    
    # Calculate metrics
    ari <- adjustedRandIndex(result$z_map, z_true)
    
    cat(sprintf("[Debug] Simulation completed for k=%d, ARI: %.3f\n", k, ari))
    
    return(list(
      success = TRUE,
      k = k,
      ari = ari,
      z_true = z_true,
      z_est = result$z_map,
      P_final = result$P_final,
      theta_final = result$theta_final,
      tau_final = result$tau_final,
      tau_samples = result$tau_samples,
      error = NULL
    ))
  }, error = function(e) {
    cat(sprintf("[Error] Simulation failed for k=%d: %s\n", k, e$message))
    return(list(
      success = FALSE,
      k = k,
      error = e$message
    ))
  })
}

# Main simulation study function
run_simulation_study <- function(n_sims = 100, k_values = c(3, 4), use_assortative = TRUE) {
  start_time <- Sys.time()
  cat(sprintf("[%s] Starting simulation study with %d cores\n", 
              format(start_time, "%H:%M:%S"),
              min(parallel::detectCores() - 1, 6)))
  cat("=====================================\n")
  
  # Setup parallel processing
  n_cores <- min(parallel::detectCores() - 1, 6)
  cl <- parallel::makeCluster(n_cores)
  cat(sprintf("[%s] Setting up cluster environment\n", format(Sys.time(), "%H:%M:%S")))
  
  # Export necessary functions and libraries to cluster
  clusterExport(cl, varlist = c("generate_network", "sbm_gibbs", "sbm_gibbs_standard", 
                                "truncated_beta", "rdirichlet", 
                                "adjustedRandIndex", "run_simulation"
  ))
  
  
  
  # Initialize results storage
  all_results <- list()
  
  # Run simulations for each k
  for (k in k_values) {
    cat(sprintf("[%s] Starting simulations for k=%d\n", format(Sys.time(), "%H:%M:%S"), k))
    
    # Split simulations into chunks
    chunk_size <- max(1, floor(n_sims/n_cores))
    chunks <- split(1:n_sims, ceiling(seq_along(1:n_sims)/chunk_size))
    
    k_results <- list()
    
    for (i in seq_along(chunks)) {
      cat(sprintf("[%s] Processing chunk %d/%d for k=%d\n", 
                  format(Sys.time(), "%H:%M:%S"), 
                  i, length(chunks), k))
      
      # Run simulations in parallel
      chunk_results <- parLapply(cl, chunks[[i]], function(dummy, k_val, ua) {
        run_simulation(k = k_val, use_assortative = ua)
      }, k_val = k, ua = use_assortative)
      
      k_results <- c(k_results, chunk_results)
      
      cat(sprintf("[%s] Completed %d/%d simulations for k=%d\n", 
                  format(Sys.time(), "%H:%M:%S"),
                  min(i*chunk_size, n_sims), n_sims, k))
    }
    
    all_results[[as.character(k)]] <- k_results
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Calculate total time
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")
  
  cat(sprintf("[%s] All simulations completed!\n", format(end_time, "%H:%M:%S")))
  cat(sprintf("[%s] Total time: %.1f minutes\n", format(end_time, "%H:%M:%S"), total_time))
  cat("=====================================\n")
  
  # Analyze results
  cat(sprintf("[%s] Starting analysis\n", format(Sys.time(), "%H:%M:%S")))
  cat(sprintf("[%s] Processing results\n", format(Sys.time(), "%H:%M:%S")))
  
  final_results <- list()
  
  for (k in k_values) {
    cat(sprintf("[%s] Analyzing results for k=%d\n", format(Sys.time(), "%H:%M:%S"), k))
    k_results <- all_results[[as.character(k)]]
    
    # Filter successful simulations
    successful_sims <- Filter(function(x) x$success, k_results)
    cat(sprintf("[%s] Valid simulations for k=%d: %d\n", 
                format(Sys.time(), "%H:%M:%S"), 
                k, length(successful_sims)))
    
    if (length(successful_sims) > 0) {
      # Calculate summary statistics
      ari_values <- sapply(successful_sims, function(x) x$ari)
      # Safely extract tau if it exists
      tau_values <- sapply(successful_sims, function(x) if (!is.null(x$tau_final)) x$tau_final else NA)
      tau_values <- na.omit(tau_values)
      
      final_results[[as.character(k)]] <- list(
        mean_ari = mean(ari_values),
        sd_ari = sd(ari_values),
        mean_tau = if (length(tau_values) > 0) mean(tau_values) else NA,
        sd_tau = if (length(tau_values) > 0) sd(tau_values) else NA,
        n_valid = length(successful_sims)
      )
      
    }
  }
  
  # Check if any valid results
  if (length(unlist(final_results)) == 0) {
    cat(sprintf("[%s] ERROR: No valid simulations completed. Please check the parameters and try again.\n",
                format(Sys.time(), "%H:%M:%S")))
  } else {
    # Create summary data frame
    summary_df <- data.frame(
      k = k_values,
      mean_ari = sapply(final_results, function(x) x$mean_ari),
      sd_ari = sapply(final_results, function(x) x$sd_ari),
      mean_tau = sapply(final_results, function(x) x$mean_tau),
      sd_tau = sapply(final_results, function(x) x$sd_tau),
      n_valid = sapply(final_results, function(x) x$n_valid)
    )
    print(summary_df)
  }
  
  # Save raw results
  cat(sprintf("[%s] Saving raw results\n", format(Sys.time(), "%H:%M:%S")))
  file_name <- if (use_assortative) "simulation_results2.rds" else "simulation_results_standard2.rds"
  
  saveRDS(list(
    raw_results = all_results,
    summary = final_results,
    parameters = list(
      n_sims = n_sims,
      k_values = k_values,
      total_time = total_time
    )
  ), file = file_name )
  
  return(final_results)
}

# Run the simulation study
results <- run_simulation_study(n_sims = 100, k_values = c(3, 4))
results_standard <- run_simulation_study(n_sims = 100, k_values = c(3, 4), use_assortative = FALSE)



#-------------------------------------------------------------------------
# Load saved RDS files
results_assort <- readRDS("simulation_results.rds")
results_standard <- readRDS("simulation_results_standard.rds")

extract_metrics <- function(results_list, k, use_assortative = TRUE) {
  sims <- results_list$raw_results[[as.character(k)]]
  
  valid <- Filter(function(x) {
    is.list(x) && !is.null(x$success) && x$success &&
      !is.null(x$z_true) && !is.null(x$z_est)
  }, sims)
  
  if (length(valid) == 0) {
    warning(paste("No valid runs found for k =", k, "model =", if (use_assortative) "Assortative" else "Standard"))
    return(data.frame())
  }
  
  # Local ARI calc (in case it's not stored)
  compute_ari <- function(x) {
    tab <- table(x$z_true, x$z_est)
    a <- sum(choose(tab, 2))
    b1 <- sum(choose(rowSums(tab), 2))
    b2 <- sum(choose(colSums(tab), 2))
    d <- choose(length(x$z_true), 2)
    exp_ri <- (b1 * b2) / d
    max_ri <- (b1 + b2) / 2
    (a - exp_ri) / (max_ri - exp_ri)
  }
  
  data.frame(
    ARI = sapply(valid, function(x) if (!is.null(x$ari)) x$ari else compute_ari(x)),
    n_comm = sapply(valid, function(x) length(unique(x$z_est))),
    model = if (use_assortative) "Assortative SBM" else "Standard SBM",
    k = k
  )
}

metrics_a3 <- extract_metrics(results_assort, 3, TRUE)
metrics_a4 <- extract_metrics(results_assort, 4, TRUE)
metrics_s3 <- extract_metrics(results_standard, 3, FALSE)
metrics_s4 <- extract_metrics(results_standard, 4, FALSE)

all_metrics <- do.call(rbind, list(metrics_a3, metrics_a4, metrics_s3, metrics_s4))

library(ggplot2)

# Histogram of estimated communities
plot_hist <- function(df, model_type, k_val) {
  df_sub <- subset(df, model == model_type & k == k_val)
  ggplot(df_sub, aes(x = n_comm)) +
    geom_histogram(binwidth = 1,
                   fill = ifelse(model_type == "Assortative SBM", "#E9BCB7", "#A2D2FF"),
                   color = "black") +
    theme_minimal() +
    labs(title = paste(model_type, "(k =", k_val, ")"),
         x = "Estimated Number of Communities", y = "Count")
}

# Boxplot of ARI
plot_box <- function(df, model_type, k_val) {
  df_sub <- subset(df, model == model_type & k == k_val)
  ggplot(df_sub, aes(x = factor(k), y = ARI, fill = model)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.4) +
    theme_minimal() +
    scale_fill_manual(values = c("Assortative SBM" = "#E9BCB7", "Standard SBM" = "#A2D2FF")) +
    labs(title = paste(model_type, "(k =", k_val, ")"),
         x = "True k", y = "Adjusted Rand Index")
}

# Histograms
plot_hist(all_metrics, "Assortative SBM", 3)
plot_hist(all_metrics, "Standard SBM", 3)
plot_hist(all_metrics, "Assortative SBM", 4)
plot_hist(all_metrics, "Standard SBM", 4)

# Boxplots
plot_box(all_metrics, "Assortative SBM", 3)
plot_box(all_metrics, "Standard SBM", 3)
plot_box(all_metrics, "Assortative SBM", 4)
plot_box(all_metrics, "Standard SBM", 4)







#Save them

output_dir <- "."

# Save histogram plots
for (model_type in c("Assortative SBM", "Standard SBM")) {
  for (k_val in c(3, 4)) {
    p <- plot_hist(all_metrics, model_type, k_val)
    fname <- paste0("hist_", tolower(gsub(" ", "_", model_type)), "_k", k_val, ".png")
    ggsave(filename = file.path("C:/Users/farno/Documents", fname), plot = p, width = 6, height = 5, dpi = 300)
  }
}

# Save boxplot plots
for (model_type in c("Assortative SBM", "Standard SBM")) {
  for (k_val in c(3, 4)) {
    p <- plot_box(all_metrics, model_type, k_val)
    fname <- paste0("boxplot_", tolower(gsub(" ", "_", model_type)), "_k", k_val, ".png")
    ggsave(filename = file.path("C:/Users/farno/Documents", fname), plot = p, width = 6, height = 5, dpi = 300)
  }
}






#------------------------------------------------------
#NODES
library(igraph)

# Example input:
# A: adjacency matrix
# z_true: true community labels
# z_map: estimated labels from MAP (from sbm_gibbs or sbm_gibbs_standard)

# Create igraph object
g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)

# Choose force-directed layout for consistent positioning
layout <- layout_with_fr(g)

# Set consistent node size and no labels
vertex_size <- 6
edge_col <- "gray80"

get_colors <- function(labels) {
  palette <- c("#E9BCB7", "#A2D2FF", "#B5EAD7", "#FFDAC1", "#FF9AA2", "#D5AAFF")
  palette[labels]
}

# Plot TRUE communities
png("graph_true_communities.png", width = 600, height = 600)
plot(g,
     layout = layout,
     vertex.color = get_colors(z_true),
     vertex.label = NA,
     vertex.size = vertex_size,
     edge.color = edge_col,
     main = "True Communities")
dev.off()

# Plot ESTIMATED communities
png("graph_estimated_communities.png", width = 600, height = 600)
plot(g,
     layout = layout,
     vertex.color = get_colors(z_map),
     vertex.label = NA,
     vertex.size = vertex_size,
     edge.color = edge_col,
     main = "Estimated Communities (MAP)")
dev.off()







#EXTRAWORK-------------------------------------------------------------
library(igraph)

# Getting one successful simulation from results
sim_list <- results_assort$raw_results[["3"]]
sim <- sim_list[[which(sapply(sim_list, function(x) x$success))[1]]]

# Getting the true and estimated labels
z_true <- sim$z_true
z_map <- sim$z_est
n <- length(z_true)

# Rebuilding the same network as in your simulation
P <- matrix(0.1, 3, 3)
diag(P) <- 0.3
A <- matrix(0, n, n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    A[i,j] <- A[j,i] <- rbinom(1, 1, P[z_true[i], z_true[j]])
  }
}

# Create the igraph object
g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
layout <- layout_with_fr(g)

# Color palette
get_colors <- function(labels) {
  palette <- c("#E9BCB7", "#A2D2FF", "#B5EAD7", "#FFDAC1", "#FF9AA2", "#D5AAFF")
  palette[labels]
}

# Plot TRUE communities
png("graph_true_communities_k3_assortative.png", width = 600, height = 600)
plot(g,
     layout = layout,
     vertex.color = get_colors(z_true),
     vertex.label = NA,
     vertex.size = 6,
     edge.color = "gray80",
     main = "True Communities (k = 3)")
dev.off()

# Plot ESTIMATED communities
png("graph_estimated_communities_k3_assortative.png", width = 600, height = 600)
plot(g,
     layout = layout,
     vertex.color = get_colors(z_map),
     vertex.label = NA,
     vertex.size = 6,
     edge.color = "gray80",
     main = "Estimated Communities (Assortative SBM, k = 3)")
dev.off()










library(igraph)

sim_list <- results_standard$raw_results[["3"]]
sim <- sim_list[[which(sapply(sim_list, function(x) x$success))[1]]]

z_true <- sim$z_true
z_map <- sim$z_est
n <- length(z_true)

P <- matrix(0.1, 3, 3)
diag(P) <- 0.3
A <- matrix(0, n, n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    A[i,j] <- A[j,i] <- rbinom(1, 1, P[z_true[i], z_true[j]])
  }
}

g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
layout <- layout_with_fr(g)

get_colors <- function(labels) {
  palette <- c("#E9BCB7", "#A2D2FF", "#B5EAD7", "#FFDAC1", "#FF9AA2", "#D5AAFF")
  palette[labels]
}

png("graph_estimated_communities_k3_standard.png", width = 600, height = 600)
plot(g,
     layout = layout,
     vertex.color = get_colors(z_map),
     vertex.label = NA,
     vertex.size = 6,
     edge.color = "gray80",
     main = "Estimated Communities (Standard SBM, k = 3)")
dev.off()











library(igraph)

sim_list <- results_assort$raw_results[["4"]]
sim <- sim_list[[which(sapply(sim_list, function(x) x$success))[1]]]

z_true <- sim$z_true
z_map <- sim$z_est
n <- length(z_true)

P <- matrix(0.1, 3, 3)
diag(P) <- 0.3
A <- matrix(0, n, n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    A[i,j] <- A[j,i] <- rbinom(1, 1, P[z_true[i], z_true[j]])
  }
}

g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
layout <- layout_with_fr(g)

get_colors <- function(labels) {
  palette <- c("#E9BCB7", "#A2D2FF", "#B5EAD7", "#FFDAC1", "#FF9AA2", "#D5AAFF")
  palette[labels]
}

png("graph_estimated_communities_k4_assortative.png", width = 600, height = 600)
plot(g,
     layout = layout,
     vertex.color = get_colors(z_map),
     vertex.label = NA,
     vertex.size = 6,
     edge.color = "gray80",
     main = "Estimated Communities (Assortative SBM, k = 4)")
dev.off()








library(igraph)

sim_list <- results_standard$raw_results[["4"]]
sim <- sim_list[[which(sapply(sim_list, function(x) x$success))[1]]]

z_true <- sim$z_true
z_map <- sim$z_est
n <- length(z_true)

P <- matrix(0.1, 3, 3)
diag(P) <- 0.3
A <- matrix(0, n, n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    A[i,j] <- A[j,i] <- rbinom(1, 1, P[z_true[i], z_true[j]])
  }
}

g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
layout <- layout_with_fr(g)

get_colors <- function(labels) {
  palette <- c("#E9BCB7", "#A2D2FF", "#B5EAD7", "#FFDAC1", "#FF9AA2", "#D5AAFF")
  palette[labels]
}

png("graph_estimated_communities_k4_standard.png", width = 600, height = 600)
plot(g,
     layout = layout,
     vertex.color = get_colors(z_map),
     vertex.label = NA,
     vertex.size = 6,
     edge.color = "gray80",
     main = "Estimated Communities (Standard SBM, k = 4)")
dev.off()

