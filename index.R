# Functions

createInitial <- function(ie, il, im, ip, iw) {
  State <- c(rep("W", iw), rep("P", ip), rep("Mec", im), rep("LI", il), rep("E", ie))
  Counter <- rep(1, iw + ip + im + il + ie)
  init_df <- data.frame(State, Counter, stringsAsFactors = F)
}

eggsBorn <- function(n_queens, spawn_gene) {
  ##Does ovarian development matter? Can a birth rate evolve?
  ##Outputs a dataframe of new eggs
  n_spawned <- round(n_queens * spawn_gene)
  spawned_eggs <- data.frame(State = rep("E", n_spawned), Counter = rep(1, n_spawned))
}

eggsDied <- function(n_aggressor, n_eggs, eateggs_gene) {
  ##eateggs_gene represents the number of eggs each aggressor eats
  n_deadeggs <- round(n_aggressor * eateggs_gene)
  if (n_deadeggs < n_eggs) {
    dead_eggs <- sample(n_eggs, n_deadeggs)
  } else {
    dead_eggs <- 1:n_eggs
  }
  dead_eggs
}

LarvaeCensus <- function(n_forager, model_type, correction, shift, uplim, diff, a) {
  ##model_type could be used to change the model from a 'high cost of entry' model
  ##to a possible 'resource exhaustion' model? correction could be used to guarantee reasonable 
  ##numbers within the function, and cap could indicate at which point the function becomes unrealistic.
  ##This function outputs the number of the surviving larvae.
  if (n_forager < shift) {
    n_survivorlarvae <- round(correction * (n_forager ^ model_type))
  } else {
    ###I have no idea? My best idea so far is continuing 
    ###the model's derivative where it left off.
    n_survivorlarvae <- uplim - (diff * exp(-a * (n_forager - shift)))
    n_survivorlarvae <- round(n_survivorlarvae)
  }
  n_survivorlarvae
}

LarvaeDied <- function(larvae_census, n_larvae) {
  #This function outputs a vector of the larvae that have died (in the vector 1:n_larvae)
  #This output must be used in conjuction with the 'larvae' vector.
  if (larvae_census >= n_larvae) {
    dead_larvae <- numeric()
  } else {
    dead_larvae <- sample(n_larvae, n_larvae - larvae_census)
  }
  dead_larvae
}

createPreferences <- function(n_eggs, n_larvae, n_mec, sig_eq, sig_lf, sig_lq, sig_la, sig_mq) {
  ##This function must be run before assigning workers b/c it creates the decision_vector variable
  p_f <- sigmoidVector(n_larvae, sig_lf)
  p_eq <- sigmoidVector(n_eggs, sig_eq)
  p_lq <- sigmoidVector(n_larvae, sig_lq)
  p_mq <- sigmoidVector(n_mec, sig_mq)
  p_a <- sigmoidVector(n_larvae, sig_la)
  p_q <- (p_eq + p_lq + p_mq) / 3
  pref_vector <- c(p_f, p_q, p_a, p_eq, p_lq, p_mq)
}

assignWorkers <- function(n_workers, decision_vector) {
  path_assign <- sample(1:6, n_workers, replace = T)
  random_assign <- matrix(runif(3 * n_workers), ncol = 3)
  real_assign <- path_options[path_assign, 1]
  undecided <- 1:n_workers
  new_undecided <- undecided
  for (i in 1:3) {
    real_assign[new_undecided] <- path_options[path_assign[new_undecided], i]
    real_assign[new_undecided] <- as.numeric(!(random_assign[new_undecided, i] > decision_vector[real_assign[new_undecided]])) * real_assign[new_undecided]
    new_undecided <- undecided[real_assign == 0]
  }
  real_assign[new_undecided] <- 4
  worker_assign <- runif(4)
  worker_assign[1] <- sum(real_assign == 1)
  worker_assign[2] <- sum(real_assign == 2)
  worker_assign[3] <- sum(real_assign == 3)
  worker_assign[4] <- sum(real_assign == 4)
  worker_assign
}

workersDied <- function(n_workers, worker_deathrate) {
  ##Note: this function is completely separate from the queensKilled function and must be run after it
  random_workers <- runif(n_workers)
  vec_deadworkers <- random_workers < worker_deathrate
  dead_workers <- which(vec_deadworkers)
}


queensKilled <- function(n_queens, n_aggressors, n_workers, killqueen_rate) {
  ##Note: this output must be used to update the n_queen and n_workers variables and the current_df
  ##Note: this function must be performed before the eggsBorn function
  n_deadqueens <- round(n_aggressors * killqueen_rate)
  if (n_deadqueens < n_queens) {
    dead_queens <- sample(n_workers, n_deadqueens)
  } else {
    dead_queens <- sample(n_workers, n_queens)
  }
  dead_queens
}

fullModel <- function(chromosome, fixed, time_steps) {
  sig_eq <- chromosome[1:4]
  sig_lf <- chromosome[5:8]
  sig_lq <- chromosome[9:12]
  sig_la <- chromosome[13:16]
  sig_mq <- chromosome[17:20]
  killqueen_rate <- chromosome[21]
  spawn_gene <- chromosome[22]
  eateggs_gene <- chromosome[23]
  model_type <- fixed[1]
  correction <- fixed[2]
  shift <- fixed[3]
  uplim <- fixed[4]
  init_eggs <- round(fixed[5])
  init_larvae <- round(fixed[6])
  init_mec <- round(fixed[7])
  init_pupae <- round(fixed[8])
  init_workers <- round(fixed[9])
  current_df <- createInitial(init_eggs, init_larvae, init_mec, init_pupae, init_workers)
  worker_deathrate <- fixed[10]
  history <- matrix(rep(0, time_steps * 26), ncol = time_steps)
  
  # Larvae Census begins
  ## This is for reducing the number of times the Lambert-W function needs to run
  shift_slope <- correction * model_type * (shift ^ (model_type - 1))
  y_shift <- correction * (shift ^ model_type)
  diff <- uplim - y_shift
  x <- shift_slope / diff
  a <- lamW::lambertW0(x)
  # Larvae Census ends
  
  for (i in 1:time_steps) {
    eggs <- which(current_df$State == "E")
    larvae <- which(current_df$State %in% larvae_class)
    mec <- which(current_df$State == "Mec")
    pupae <- which(current_df$State == "P")
    workers <- which(current_df$State == "W")
    n_eggs <- length(eggs)
    n_larvae <- length(larvae)
    n_mec <- length(mec)
    n_pupae <- length(pupae)
    n_workers <- length(workers)
    
    # Number of ants at each stage at the beginning of the time-step
    history[1, i] <- n_eggs
    history[2, i] <- n_larvae
    history[3, i] <- n_mec
    history[4, i] <- n_pupae
    history[5, i] <- n_workers
    
    # Number of ants in the colony at the beginning of the time-step
    history[16, i] <- sum(history[1:5, i])
    
    pref_vector <- createPreferences(n_eggs, n_larvae, n_mec, sig_eq, sig_lf, sig_lq, sig_la, sig_mq)
    if (n_workers == 0) {
      worker_assign <- rep(0, 4)
    } else {
      worker_assign <- assignWorkers(n_workers, pref_vector[1:3])
    }
    
    # Numbers to determine how workers are assigned to roles
    ## 1/2/3 -> PF/PQ/PA, i.e. actual decision vector
    ## 4/5/6 -> PEQ/PLQ/PMQ, i.e. numbers for determining PQ
    history[21:26, i] <- pref_vector
    
    # Percentage of workers in each active role during the time-step
    history[6, i] <- worker_assign[1]
    history[7, i] <- worker_assign[2]
    history[8, i] <- worker_assign[3]
    
    dead_queens <- queensKilled(worker_assign[2], worker_assign[3], n_workers, killqueen_rate) 
    if (length(dead_queens) != 0) {
      deadqueens_pos <- workers[dead_queens]
      worker_assign[2] <- worker_assign[2] - length(dead_queens)
      workers <- workers[-dead_queens]
      n_workers <- length(workers)
    } else {
      deadqueens_pos <- integer(0)
    }
    spawned_eggs <- eggsBorn(worker_assign[2], spawn_gene)
    
    # Ants coming into the colony
    history[9, i] <- nrow(spawned_eggs)
    
    current_df <- rbind(current_df, spawned_eggs)
    eggs <- which(current_df$State == "E")
    n_eggs <- length(eggs)
    dead_eggs <- eggsDied(worker_assign[3], n_eggs, eateggs_gene)
    larvae_census <- LarvaeCensus(worker_assign[1], model_type, correction, shift, uplim, diff, a)
    dead_larvae <- LarvaeDied(larvae_census, n_larvae)
    dead_workers <- workersDied(n_workers, worker_deathrate)
    
    # Larvae that survive the time-step
    history[10, i] <- min(c(larvae_census, n_larvae))
    
    # Ants leaving the colony
    history[11, i] <- length(dead_eggs)
    history[12, i] <- length(dead_larvae)
    history[13, i] <- length(dead_workers)
    history[14, i] <- length(dead_queens)
    
    # Overall change in colony size within the time-step
    history[15, i] <- history[9, i] - sum(history[11:14, i])
    
    dead_ants <- c(eggs[dead_eggs], larvae[dead_larvae], workers[dead_workers], deadqueens_pos)
    dead_ants <- dead_ants[dead_ants != 0]
    dead_ants <- dead_ants[!is.na(dead_ants)]
    if (length(dead_ants) == 0) {
      current_df
    } else {
      current_df <- current_df[-dead_ants, ]
    }
    continue_pos <- which(current_df$Counter == as.numeric(growing_up[current_df$State]))
    
    # Number of ants moving from one stage to another at the end of the time-step
    ## E -> L
    history[17, i] <- length(which(current_df$State[continue_pos] == "E"))
    ## L -> Mec
    history[18, i] <- length(which(current_df$State[continue_pos] %in% larvae_class))
    ## Mec -> P
    history[19, i] <- length(which(current_df$State[continue_pos] == "Mec"))
    ## P -> W
    history[20, i] <- length(which(current_df$State[continue_pos] == "P"))
    
    current_df$State[continue_pos] <- as.character(next_step[current_df$State[continue_pos]])
    current_df$Counter[continue_pos] <- rep(1, length(continue_pos))
    if (length(continue_pos) == 0) {
      current_df$Counter <- current_df$Counter + 1
    } else {
      current_df$Counter[-continue_pos] <- current_df$Counter[-continue_pos] + 1
    }
  }
  history
}

fullModelM <- function(chromosome, fixed, time_steps = 365) {
  history <- fullModel(chromosome, fixed, time_steps)
  res <- mean(history[16, ])
}

sigmoid <- function(x, a = 0, k = 1, b = 0.1, m = 100) {
  a + (k - a) / ((1 + exp(-b * (x - m))))
}

sigmoidVector <- function(number, vector) {
  sigmoid(number, a = vector[1], k = vector[2], b = vector[3], m = vector[4])
}

path_options <- matrix(c(1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1), ncol = 3, byrow = T)

growing_up <- c(8, 2, 2, 2, 10, 3, 15, 0)
names(growing_up) <- c("E", "LI", "LII", "LIII", "LIV", "Mec", "P", "W")

next_step <- c("LI", "LII", "LIII", "LIV", "Mec", "P", "W")
names(next_step) <- c("E", "LI", "LII", "LIII", "LIV", "Mec", "P")

larvae_class <- c("LI", "LII", "LIII", "LIV")

testLS <- function(test_points, model_type, correction, change, cap) {
  output <- rep(0, test_points)
  for (i in 1:test_points) {
    output[i] <- LarvaeCensus(i, model_type, correction, change, cap)
  }
  output
}

colonialOptimizer <- function(generations, suggestions = NULL) {
  GA <- ga(type = "real-valued", fullModelM, fixed, min = min, 
           max = max, maxiter = generations, keepBest = T, suggestions = suggestions)
}

min <- c(rep(0, times = 23))
max <- c(rep(c(1, 1, 2, 100000), times = 5), 20, 6, 0)
fixed <- c(3, 0.1, 50, 15000, 0, 0, 0, 0, 100, 0.02)