################################################################################
## REVISING THE BORGATTI-EVERETT CORE-PERIPHERY MODEL
## (3) Functions
## R script written by José Luis Estévez (University of Helsinki)
## Date: Nov 20th, 2024
################################################################################

# Core-periphery functions: Adapted from netUtils package

# Standard core-periphery model as in UCInet with delta for inter-cat blocks
cp.ucinet <- function(graph,delta=NA, ...){
  A <- igraph::as_adjacency_matrix(graph, type = "both", sparse = FALSE)
  ev <- igraph::degree(graph, mode = "all", loops = FALSE)
  thresh <- unique(ev)
  optcorr <- -2
  optperm <- NULL
  
  for (tr in thresh) {
    evabs <- as.integer(ev >= tr) # 1 (core), 0 (periphery)
    # if all members are assigned to the core (or periphery), jump to next iteration
    if (all(evabs == 1 | all(evabs == 0))){
      next
    }
    
    E <- outer(evabs, evabs, "+")
    E[E == 1] <- delta
    E[E == 2] <- 1
    diag(E) <- NA
    if (sum(E, na.rm = TRUE) == 0) {
      (next)()
    }
    tmp <- suppressWarnings(graph_cor(E, A))
    if (is.na(tmp)) {
      (next)()
    }
    if (tmp > optcorr) {
      optperm <- evabs
      optcorr <- tmp
    }
  }
  return(list(vec = optperm, corr = round(optcorr,4)))
}

# Exact density blocks
cp.exactden <- function(graph, delta, ...) {
  A <- igraph::as_adjacency_matrix(graph, type = "both", sparse = FALSE)
  ev <- igraph::degree(graph, mode = "all", loops = FALSE)
  thresh <- unique(ev)
  optcorr <- -2
  optperm <- NULL
  
  for (tr in thresh) {
    evabs <- as.integer(ev >= tr)  # 1 (core), 0 (periphery)
    # if all members are assigned to the core (or periphery), jump to next iteration
    if (all(evabs == 1 | all(evabs == 0))){
      next
    }
    
    # Ensure dimensions 
    n_core <- sum(evabs == 1)
    n_periphery <- sum(evabs == 0)
    
    # Define the blocks in A
    Acc <- A[evabs == 1, evabs == 1]
    App <- A[evabs == 0, evabs == 0]
    Acp <- matrix(sort(A[evabs == 1, evabs == 0], decreasing = TRUE), nrow = n_core, byrow = TRUE)  # sorted vector
    Apc <- matrix(sort(A[evabs == 0, evabs == 1], decreasing = TRUE), ncol = n_core, byrow = FALSE)  # sorted vector
    
    # Putting all blocks back together for AA
    if (n_core == 1){
      AA <- rbind(c(Acc, Acp), cbind(Apc, App))
    } else if (n_periphery == 1){
      AA <- rbind(cbind(Acc, Acp), cbind(Apc, App))
    } else {
      AA <- matrix(0, nrow = n_core + n_periphery, ncol = n_core + n_periphery)
      AA[1:n_core, 1:n_core] <- Acc
      AA[(n_core + 1):(n_core + n_periphery), (n_core + 1):(n_core + n_periphery)] <- App
      AA[1:n_core, (n_core + 1):(n_core + n_periphery)] <- Acp
      AA[(n_core + 1):(n_core + n_periphery), 1:n_core] <- Apc
    } 
    diag(AA) <- NA
    
    # Ideal pattern blocks
    Ecc <- matrix(1, nrow = n_core, ncol = n_core)  # full core
    Epp <- matrix(0, nrow = n_periphery, ncol = n_periphery)  # empty periphery
    
    # Inter-categorical block images
    n <- length(Acp)
    ones <- round(delta * n)  # rounded to nearest integer
    zeroes <- n - ones
    Ecp <- matrix(c(rep(1, ones), rep(0, zeroes)), nrow = n_core, byrow = TRUE)
    Epc <- matrix(c(rep(1, ones), rep(0, zeroes)), ncol = n_core, byrow = FALSE)
    
    # Putting all ideal blocks back together
    if (n_core == 1){
      E <- rbind(c(Ecc, Ecp), cbind(Epc, Epp))
    } else if (n_periphery == 1){
      E <- rbind(cbind(Ecc, Ecp), c(Epc, Epp))
    } else {
      E <- matrix(0, nrow = n_core + n_periphery, ncol = n_core + n_periphery)
      E[1:n_core, 1:n_core] <- Ecc
      E[(n_core + 1):(n_core + n_periphery), (n_core + 1):(n_core + n_periphery)] <- Epp
      E[1:n_core, (n_core + 1):(n_core + n_periphery)] <- Ecp
      E[(n_core + 1):(n_core + n_periphery), 1:n_core] <- Epc
    } 
    diag(E) <- NA
    
    if (sum(E, na.rm = TRUE) == 0) {
      next
    }
    tmp <- suppressWarnings(cor(as.vector(AA),as.vector(E), use='complete.obs'))
    if (is.na(tmp)) {
      next
    }
    if (tmp > optcorr) {
      optperm <- evabs
      optcorr <- tmp
    }
  }
  return(list(vec = optperm, corr = round(optcorr, 4)))
}

# Minimum density blocks
cp.minden <- function(graph, delta=0, ...) {
  A <- igraph::as_adjacency_matrix(graph, type = "both", sparse = FALSE)
  ev <- igraph::degree(graph, mode = "all", loops = FALSE)
  thresh <- unique(ev)
  optcorr <- -2
  optperm <- NULL
  
  for (tr in thresh) {
    evabs <- as.integer(ev >= tr)  # 1 (core), 0 (periphery)
    # if all members are assigned to the core (or periphery), jump to next iteration
    if (all(evabs == 1 | all(evabs == 0))){
      next
    }
    
    # Ensure dimensions 
    n_core <- sum(evabs == 1)
    n_periphery <- sum(evabs == 0)
    
    # Define the blocks in A
    Acc <- A[evabs == 1, evabs == 1]
    App <- A[evabs == 0, evabs == 0]
    Acp <- matrix(sort(A[evabs == 1, evabs == 0], decreasing = TRUE), nrow = n_core, byrow = TRUE)
    Apc <- matrix(sort(A[evabs == 0, evabs == 1], decreasing = TRUE), ncol = n_core, byrow = FALSE)
                         
    # Putting all blocks back together for AA
    if (n_core == 1){
      AA <- rbind(c(Acc, Acp), cbind(Apc, App))
    } else if (n_periphery == 1){
      AA <- rbind(cbind(Acc, Acp), c(Apc, App))
    } else {
      AA <- matrix(0, nrow = n_core + n_periphery, ncol = n_core + n_periphery)
      AA[1:n_core, 1:n_core] <- Acc
      AA[(n_core + 1):(n_core + n_periphery), (n_core + 1):(n_core + n_periphery)] <- App
      AA[1:n_core, (n_core + 1):(n_core + n_periphery)] <- Acp
      AA[(n_core + 1):(n_core + n_periphery), 1:n_core] <- Apc
    } 
    diag(AA) <- NA
    
    # Ideal pattern blocks
    Ecc <- matrix(1, nrow = n_core, ncol = n_core)  # full core
    Epp <- matrix(0, nrow = n_periphery, ncol = n_periphery)  # empty periphery
    
    # Inter-categorical block images
    n <- length(Acp)
    ones <- ceiling(delta * n)  # rounded up
    Ecp <- t(Acp)  # use the transposed for allocating the values more easily
    Ecp[1:ones] <- 1
    Ecp <- t(Ecp)  # return to its original shape
    Epc <- Apc
    Epc[1:ones] <- 1
    
    # Putting all ideal blocks back together
    if (n_core == 1){
      E <- rbind(c(Ecc, Ecp), cbind(Epc, Epp))
    } else if (n_periphery == 1){
      E <- rbind(cbind(Ecc, Ecp),c(Epc, Epp))
    } else {
      E <- matrix(0, nrow = n_core + n_periphery, ncol = n_core + n_periphery)
      E[1:n_core, 1:n_core] <- Ecc
      E[(n_core + 1):(n_core + n_periphery), (n_core + 1):(n_core + n_periphery)] <- Epp
      E[1:n_core, (n_core + 1):(n_core + n_periphery)] <- Ecp
      E[(n_core + 1):(n_core + n_periphery), 1:n_core] <- Epc
    } 
    diag(E) <- NA
    
    if (sum(E, na.rm = TRUE) == 0) {
      next
    }
    tmp <- suppressWarnings(cor(as.vector(E), as.vector(AA), use='complete.obs'))
    if (is.na(tmp)) {
      next
    }
    if (tmp > optcorr) {
      optperm <- evabs
      optcorr <- tmp
    }
  }
  return(list(vec = optperm, corr = round(optcorr, 4)))
}

# p-core implementation with standard treatment of inter-cat blocks
cp.pcore <- function(graph, delta=NA, p, ...) {
  A <- igraph::as_adjacency_matrix(graph, type = "both", sparse = FALSE)
  ev <- igraph::degree(graph, mode = "all", loops = FALSE)
  thresh <- unique(ev)
  optcorr <- -2
  optperm <- NULL
  
  for (tr in thresh) {
    evabs <- as.integer(ev >= tr)  # 1 (core), 0 (periphery)
    # if all members are assigned to the core (or periphery), jump to next iteration
    if (all(evabs == 1 | all(evabs == 0))){
      next
    }
    
    # Ensure dimensions 
    n_core <- sum(evabs == 1)
    n_periphery <- sum(evabs == 0)
    
    # Define the blocks in A
    Acc <- A[evabs == 1, evabs == 1]
    # Two images of the block
    if(n_core > 1){
      # Core image 1, sorted by row
      Acc1 <- t(apply(Acc, 1, sort, decreasing = TRUE)) 
      Acc1[,n_core] <- NA
      # Core image 2, sorted by column
      Acc2 <- apply(Acc, 2, sort, decreasing = TRUE)
      Acc2[n_core,] <- NA
    }else{
      Acc1 <- Acc2 <- Acc
    }
    App <- A[evabs == 0, evabs == 0]
    Acp <- matrix(A[evabs == 1, evabs == 0], nrow = n_core)
    Apc <- matrix(A[evabs == 0, evabs == 1], ncol = n_core)
    
    # Putting all blocks back together for AA
    if (n_core == 1){
      AA <- rbind(c(Acc, Acp), cbind(Apc, App))
    } else if (n_periphery == 1){
      AA <- rbind(cbind(Acc, Acp), c(Apc, App))
    } else {
      AA <- matrix(0, nrow = n_core + n_periphery, ncol = n_core + n_periphery)
      AA[1:n_core, 1:n_core] <- Acc
      AA[(n_core + 1):(n_core + n_periphery), (n_core + 1):(n_core + n_periphery)] <- App
      AA[1:n_core, (n_core + 1):(n_core + n_periphery)] <- Acp
      AA[(n_core + 1):(n_core + n_periphery), 1:n_core] <- Apc
    } 
    diag(AA) <- NA
    
    # Ideal pattern blocks
    # The core images
    Ecc <- matrix(1, nrow = n_core, ncol = n_core)  # full core
    n <- nrow(Acc) - 1
    ones <- ceiling(p * n) # rounded up
    if (n_core > 1){
      Ecc1 <- Acc1
      Ecc1[,1:ones] <- 1 # first columns filled with ones
      Ecc2 <- Acc2
      Ecc2[1:ones,] <- 1 # first rows filled with ones
    }else{
      Ecc1 <- Ecc2 <- Ecc
    }
    Epp <- matrix(0, nrow = n_periphery, ncol = n_periphery)  # empty periphery
    
    # Inter-categorical block images
    Ecp <- matrix(delta,nrow = nrow(Acp),ncol = ncol(Acp))
    Epc <- t(Ecp)
    
    # Putting all ideal blocks back together
    if (n_core == 1){
      E <- rbind(c(Ecc, Ecp), cbind(Epc, Epp))
    } else if (n_periphery == 1){
      E <- rbind(cbind(Ecc, Ecp),c(Epc, Epp))
    } else {
      E <- matrix(0, nrow = n_core + n_periphery, ncol = n_core + n_periphery)
      E[1:n_core, 1:n_core] <- Ecc
      E[(n_core + 1):(n_core + n_periphery), (n_core + 1):(n_core + n_periphery)] <- Epp
      E[1:n_core, (n_core + 1):(n_core + n_periphery)] <- Ecp
      E[(n_core + 1):(n_core + n_periphery), 1:n_core] <- Epc
    } 
    diag(E) <- NA
    
    # Let's create weights
    W <- rep(1,length(AA))
    Wcc1 <- rep(0.5,length(Acc1)) # core observations are weighted by half
    Wcc2 <- rep(0.5,length(Acc2))
    
    if (sum(E, na.rm = TRUE) == 0) {
      next
    }
    # All together into dataset 
    X <- data.frame(a = c(as.vector(Acc1),as.vector(Acc2),as.vector(AA)),
                    e = c(as.vector(Ecc1),as.vector(Ecc2),as.vector(E)),
                    w = c(as.vector(Wcc1),as.vector(Wcc2),as.vector(W)))
    X <- na.omit(X) # remove NAs
    tmp <- suppressWarnings(weightedCorr(x=X$e,y=X$a,weights=X$w,method='Pearson'))
    if (is.na(tmp)) {
      next
    }
    if (tmp > optcorr) {
      optperm <- evabs
      optcorr <- tmp
    }
  }
  return(list(vec = optperm, corr = round(optcorr, 4)))
}

# p-core implementation with minimum density blocks
cp.minden.pcore <- function(graph, delta=0, p, ...) {
  A <- igraph::as_adjacency_matrix(graph, type = "both", sparse = FALSE)
  ev <- igraph::degree(graph, mode = "all", loops = FALSE)
  thresh <- unique(ev)
  optcorr <- -2
  optperm <- NULL
  
  for (tr in thresh) {
    evabs <- as.integer(ev >= tr)  # 1 (core), 0 (periphery)
    # if all members are assigned to the core (or periphery), jump to next iteration
    if (all(evabs == 1 | all(evabs == 0))){
      next
    }
    
    # Ensure dimensions 
    n_core <- sum(evabs == 1)
    n_periphery <- sum(evabs == 0)
    
    # Define the blocks in A
    Acc <- A[evabs == 1, evabs == 1]
    # Two images of the block
    if(n_core > 1){
      # Core image 1, sorted by row
      Acc1 <- t(apply(Acc, 1, sort, decreasing = TRUE)) 
      Acc1[,n_core] <- NA
      # Core image 2, sorted by column
      Acc2 <- apply(Acc, 2, sort, decreasing = TRUE)
      Acc2[n_core,] <- NA
    }else{
      Acc1 <- Acc2 <- Acc
    }
    App <- A[evabs == 0, evabs == 0]
    Acp <- matrix(sort(A[evabs == 1, evabs == 0], decreasing = TRUE), nrow = n_core, byrow = TRUE)
    Apc <- matrix(sort(A[evabs == 0, evabs == 1], decreasing = TRUE), ncol = n_core, byrow = FALSE)
    
    # Putting all blocks back together for AA
    if (n_core == 1){
      AA <- rbind(c(Acc, Acp), cbind(Apc, App))
    } else if (n_periphery == 1){
      AA <- rbind(cbind(Acc, Acp), c(Apc, App))
    } else {
      AA <- matrix(0, nrow = n_core + n_periphery, ncol = n_core + n_periphery)
      AA[1:n_core, 1:n_core] <- Acc
      AA[(n_core + 1):(n_core + n_periphery), (n_core + 1):(n_core + n_periphery)] <- App
      AA[1:n_core, (n_core + 1):(n_core + n_periphery)] <- Acp
      AA[(n_core + 1):(n_core + n_periphery), 1:n_core] <- Apc
    } 
    diag(AA) <- NA
    
    # Ideal pattern blocks
    # The core images
    Ecc <- matrix(1, nrow = n_core, ncol = n_core)  # full core
    n <- nrow(Acc) - 1
    ones <- ceiling(p * n) # rounded up
    if (n_core > 1){
      Ecc1 <- Acc1
      Ecc1[,1:ones] <- 1 # first columns filled with ones
      Ecc2 <- Acc2
      Ecc2[1:ones,] <- 1 # first rows filled with ones
    }else{
      Ecc1 <- Ecc2 <- Ecc
    }
    Epp <- matrix(0, nrow = n_periphery, ncol = n_periphery)  # empty periphery
    
    # Inter-categorical block images
    n <- length(Acp)
    ones <- ceiling(delta * n)  # rounded up
    Ecp <- t(Acp)  # use the transposed for allocating the values more easily
    Ecp[1:ones] <- 1
    Ecp <- t(Ecp)  # return to its original shape
    Epc <- Apc
    Epc[1:ones] <- 1
    
    # Putting all ideal blocks back together
    if (n_core == 1){
      E <- rbind(c(Ecc, Ecp), cbind(Epc, Epp))
    } else if (n_periphery == 1){
      E <- rbind(cbind(Ecc, Ecp),c(Epc, Epp))
    } else {
      E <- matrix(0, nrow = n_core + n_periphery, ncol = n_core + n_periphery)
      E[1:n_core, 1:n_core] <- Ecc
      E[(n_core + 1):(n_core + n_periphery), (n_core + 1):(n_core + n_periphery)] <- Epp
      E[1:n_core, (n_core + 1):(n_core + n_periphery)] <- Ecp
      E[(n_core + 1):(n_core + n_periphery), 1:n_core] <- Epc
    } 
    diag(E) <- NA
    
    # Let's create weights
    W <- rep(1,length(AA))
    Wcc1 <- rep(0.5,length(Acc1)) # core observations are weighted by half
    Wcc2 <- rep(0.5,length(Acc2))
    
    if (sum(E, na.rm = TRUE) == 0) {
      next
    }
    # All together into dataset 
    X <- data.frame(a = c(as.vector(Acc1),as.vector(Acc2),as.vector(AA)),
                    e = c(as.vector(Ecc1),as.vector(Ecc2),as.vector(E)),
                    w = c(as.vector(Wcc1),as.vector(Wcc2),as.vector(W)))
    X <- na.omit(X) # remove NAs
    tmp <- suppressWarnings(weightedCorr(x=X$e,y=X$a,weights=X$w,method='Pearson'))
    if (is.na(tmp)) {
      next
    }
    if (tmp > optcorr) {
      optperm <- evabs
      optcorr <- tmp
    }
  }
  return(list(vec = optperm, corr = round(optcorr, 4)))
}

################################################################################

# Let's create a function to simulate core-periphery networks

# cp(N,d,p,k)
# N; number of nodes
# c; number of core-nodes
# p; probability of a tie, a value between zero and one [0,1]
# k; [1,1/sqrt(p)]. If 1, Erdos-Renyi networks. If 1/sqrt(p), clear core-periphery network structure

# If a core-periphery structure looks like this:
# A_cc - A_cp
# A_pc - A_pp
# The probability of a tie in each block is given by:
# k^2p - kp
# kp   - p

random.cp <- function(N = 50, c = 25, p = 1/9, k = 1) {
  # Conditions to run the algorithm
  if (c >= N) {
    stop('c must be smaller than N')
  }
  if (p < 0 | p > 1) {
    stop('p must contain a probability: [0,1]')
  }
  if (k < 1 | k > (1/p)^(1/2)) {
    stop('k must be in the range [1,1/sqrt(p)]')
  }
  
  # Create the matrix
  vct <- c(rep(1, c), rep(0, N - c))  # memberships (1 = core, 0 = periphery)
  names <- paste(ifelse(vct == 1, 'C', 'P'), seq_along(vct), sep = '')
  mtx <- matrix(NA, nrow = N, ncol = N, dimnames = list(names, names))  # empty matrix
  
  for (i in seq_along(vct)) {
    for (j in seq_along(vct)) {
      if (i >= j) {  # leave diagonal empty
        if (vct[i] == 1 & vct[j] == 1) {  # if cc
          mtx[i, j] <- rbinom(1, 1, k^2 * p)
        } else if (vct[i] == 0 & vct[j] == 0) {  # if pp
          mtx[i, j] <- rbinom(1, 1, p)
        } else {  # if cp or pc
          mtx[i, j] <- rbinom(1, 1, k * p)
        }
      }
    }
  }
  
  # Turn into igraph format
  grp <- igraph::graph_from_adjacency_matrix(mtx, mode = 'lower', diag = FALSE)
  return(grp)
}

# Function to create networks with two core-periphery structures

# cp(N,d,p,k,connected)
# N; number of nodes
# c; number of core-nodes
# p; probability of a tie, a value between zero and one [0,1]
# k; [1,1/sqrt(p)]. If 1, Erdos-Renyi networks. If 1/sqrt(p), clear core-periphery network structure
# connected; How are the two core-periphery structures connected? 'no', 'core', 'periphery'
# Depending on the version chosen, the networks are created as follows:
# Unconnected ('no'), the default option
# k^2p -  kp  -  0    -  0
# kp   -  p   -  0    -  0
# 0    -  0   -  k^2p -  kp
# 0    -  0   -  kp   -  p
# Connected via core nodes '(core') 
# k^2p      -  kp  -  p/sqrt(c) -  0
# kp        -  p   -  0         -  0
# p/sqrt(c) -  0   -  k^2p      -  kp
# 0         -  0   -  kp        -  p
# Connected via peripheral nodes ('periphery)
# k^2p -  kp          -  0    -  0
# kp   -  p           -  0    -  p/sqrt(N-c)
# 0    -  0           -  k^2p -  kp
# 0    -  p/sqrt(N-c) -  kp   -  p

random.cp2 <- function(N = 50, c = 25, p = 1/9, k = 1, connected='no') {
  # Conditions to run the algorithm
  if (N %% 2 == 1) {
    stop('Choose an even number for N') 
  }
  if (c >= N) {
    stop('c must be smaller than N (an preferably also an even number)')
  }
  if (p < 0 | p > 1) {
    stop('p must contain a probability: [0,1]')
  }
  if (k < 1 | k > (1/p)^(1/2)) {
    stop('k must be in the range [1,1/sqrt(p)]')
  }
  # If c is an odd number, rounded it up
  if(c %% 2 == 1){
    c <- ceiling(c)
  }
  
  # Half the sample and half the number of core members
  N_2 <- N/2
  c_2 <- c/2
  
  # Create the matrix
  vct <- c(rep(1,c_2),rep(0,N_2-c_2)) # memberships (1 = core, 0 = periphery)
  names1 <- paste(ifelse(vct == 1, 'C', 'P'), seq_along(vct), sep = '')
  names2 <- paste(ifelse(vct == 1, 'C', 'P'), N_2 + seq_along(vct), sep = '')
  
  # Create on-diagonal matrices
  mtx1 <- mtx3 <- matrix(NA, nrow = N_2, ncol = N_2)
  
  for (i in seq_along(vct)) {
    for (j in seq_along(vct)) {
      if (i >= j) {  # fill only the lower part
        if (vct[i] == 1 & vct[j] == 1) {  # if cc
          mtx1[i, j] <- rbinom(1, 1, k^2 * p)
          mtx3[i, j] <- rbinom(1, 1, k^2 * p)
        } else if (vct[i] == 0 & vct[j] == 0) {  # if pp
          mtx1[i, j] <- rbinom(1, 1, p)
          mtx3[i, j] <- rbinom(1, 1, p)
        } else {  # if cp or pc
          mtx1[i, j] <- rbinom(1, 1, k * p)
          mtx3[i, j] <- rbinom(1, 1, k * p)
        }
      }
    }
  }
  
  # Create the off-diagonal matrix
  mtx2 <- matrix(0, nrow = N_2, ncol = N_2)
  
  for (i in seq_along(vct)) {
    for (j in seq_along(vct)) {
      if (connected == 'core') { # connected through cores
        if (vct[i] == 1 & vct[j] == 1) { mtx2[i,j] <- rbinom(1, 1, p/sqrt(c)) } 
      }else if (connected == 'periphery') { # connected through peripheries
        if (vct[i] == 0 & vct[j] == 0) { mtx2[i,j] <- rbinom(1, 1, p/sqrt(N-c)) } 
      }
    }
  }
  
  # Put everything into a larger matrix
  mtx <- matrix(NA, nrow = N, ncol = N,
                dimnames = list(c(names1,names2),c(names1,names2)))  # empty matrix
  mtx[1:(N_2), 1:(N_2)] <- mtx1
  mtx[(N_2+1):N, 1:(N_2)] <- mtx2
  mtx[(N_2+1):N, (N_2+1):N] <- mtx3
  diag(mtx) <- NA # remove the diagonal
  
  # Turn into igraph format
  grp <- igraph::graph_from_adjacency_matrix(mtx, mode = 'lower', diag = FALSE)
  return(grp)
}

################################################################################