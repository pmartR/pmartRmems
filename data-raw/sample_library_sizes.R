
sample_library_sizes <- function(target_lib_sz, desired_base_library_size = NULL) {
  n_cols <- length(target_lib_sz)
  empirical_scaling <- target_lib_sz / median(target_lib_sz)
  return(sample(round(desired_base_library_size * empirical_scaling), size = n_cols, replace = TRUE))
}

# Given a data matrix, target_data, simulate a similar data set, sim_data,
# with known True Positives.
# The number of True Positives, nTP, is set by the user.
# The effect size, effect_size, is set by the user or the effect_size_scaling_factor is
# set, resulting in a different random effect size for each Truely Positive row.

simulate_count_data <- function(target_data, treatment_cols = NULL, nTP = 30, effect_size = 10,
                                effect_size_scaling_factor = NULL, transfun = NULL,
                                shift = FALSE, base_data = NULL) {
  if (is.null(base_data)) {
    if (is.null(treatment_cols)) stop("Please specify treatment columns")
  }
  # browser()
  if (is.null(base_data)) {
    n_rows <- dim(target_data)[1]
    n_cols <- dim(target_data)[2]
    sim_data <- matrix(rep(-10, n_rows*n_cols), nrow = n_rows)

    ### Library size simulation:
    target_lib_sz <- apply(target_data, 2, sum)
    sim_lib_sz <- sample_library_sizes(target_lib_sz, median(target_lib_sz))

    ### Combine counts over samples to create a multinomial probability vector (non-normalized)
    mult_probs <- apply(target_data, 1, sum)

    for(i in 1:n_cols) {
      sim_data[,i] <- rmultinom(1, size = sim_lib_sz[i], prob = mult_probs)
    }

    ### add in some zeros (low count thresholding)

    # original percent zero
    opz <- mean(target_data == 0)
    # simulated percent zero
    spz <- mean(sim_data == 0)
    # fraction of simulated non-zero entries to be set to zero
    r <- (opz - spz) / (1 - spz)
    # vnzd: vector of non-zero simulated data
    sida <- as.vector(sim_data)
    snzd <- sida[sida != 0]
    # p[i] is the probability of setting a count of i to be 0
    # c is the scaling factor: p[i]=c/i
    ### scaling factor can be altered using argument 'transfun'
    if (is.null(transfun)) {
      transfun <- function(i) return(i)
    }
    i_transform <- transfun(1:max(snzd))

    sm <- 0
    for(i in 1:max(snzd)) {
      sm <- sm + mean(snzd == i)/i_transform[i]
    }

    c <- r / sm
    p <- rep(-10, max(snzd))
    for(i in 1:max(snzd)) {
      p[i] <- c/i_transform[i]
    }
    # and probabilities can't be greater than 1 so...
    p[p > 1] <- 1

    # introduce the zeros into the simulated data set
    zeros_mat <- matrix(rep(1, n_rows*n_cols), nrow = n_rows)
    for(i in 1:n_rows) {
      for(j in 1:n_cols) {
        if(sim_data[i,j] != 0) {
          zeros_mat[i,j] <- 1 - rbinom(1, 1, prob = p[sim_data[i,j]])
        }
      }
    }
    sim_data <- sim_data * zeros_mat
    sim_data_base <- sim_data

  } else {
    sim_data <- base_data
    sim_data_base <- base_data
  }

  ### Add signal
  # Number of true-positives per sample: nTP
  # Effect Size, a multiplicative fold effect: effect_size

  if (!is.null(effect_size_scaling_factor))
    effect_size <- exp(rnorm(nTP, mean = 0, sd = effect_size_scaling_factor))

  # Randomly sample which rows will be increased
  significant_rows <- sample(c(1:nrow(sim_data)), nTP)

  # apply signal to appropriate rows
  ## define groups
  g1 <- treatment_cols
  g2 <- which(!(1:ncol(sim_data) %in% treatment_cols))
  ## subset significant data
  sig_data <- sim_data[significant_rows,]
  ## arrange in descending order of rowsums
  ord <- order(rowSums(sig_data), decreasing = TRUE)
  sigdata_sort <- sig_data[ord,]
  sigrows_sort <- significant_rows[ord]
  ## add signal, alternating groups every other row
  seq1 <- seq(1, nrow(sigdata_sort), by = 2)
  seq2 <- seq(2, nrow(sigdata_sort), by = 2)
  sigdata_sort[seq1, g1] <- sigdata_sort[seq1, g1] * effect_size
  sigdata_sort[seq2, g2] <- sigdata_sort[seq2, g2] * effect_size
  ## TEST
  ### zeros are throwing off p-values, so this is an attempt to lessen their influence.
  if (shift == TRUE) {
    lg1 <- sigdata_sort[seq1, g1] == 0
    sigdata_sort[seq1, g1][lg1] <- sigdata_sort[seq1, g1][lg1] + effect_size/2
    lg2 <- sigdata_sort[seq2, g2] == 0
    sigdata_sort[seq2, g2][lg2] <- sigdata_sort[seq2, g2][lg2] + effect_size/2
    # sigdata_sort[sigdata_sort == 0] <- sigdata_sort[sigdata_sort == 0] + effect_size/2
  }
  ## END TEST
  ## place transformed rows back in sim_data
  sim_data[sigrows_sort,] <- sigdata_sort

  # sim_data[significant_rows, treatment_cols] <-
  #   round(sim_data[significant_rows, treatment_cols]*effect_size)

  # return both the simulated data and the vector of significant rows
  if (!is.null(effect_size_scaling_factor)) {
    return_me <- list(sim_data = sim_data, significant_rows = sort(significant_rows),
                      effect_size = effect_size[order(significant_rows)])
    return(return_me)
  }
  return_me <- list(sim_data = sim_data, significant_rows = sort(significant_rows),
                    sim_data_base = sim_data_base)

  return(return_me)
}



# t test
ttest_rrna <- function(edata, paired = TRUE, groups = NULL) {

  # split data into groups
  ncols <- ncol(edata)
  nrows <- nrow(edata)
  # browser()
  if (is.null(groups)) {
    df1 <- edata[,1:floor(ncols/2)]
    df2 <- edata[,floor(ncols/2):ncols]
  } else {
    df1 <- edata[,groups[[1]]]
    df2 <- edata[,groups[[2]]]
  }

  # run tests, return p-values
  if (!paired) {
    return(sapply(1:nrows, function(i) t.test(df2[i,], df1[i,])$p.value))
  } else {
    if (ncol(df1) != ncol(df2)) stop("Not paired data")
    return(sapply(1:nrows, function(i) t.test(df2[i,] - df1[i,])$p.value))
  }
}
