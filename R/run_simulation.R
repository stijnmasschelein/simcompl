#' Run data simulations and tests together
#'
#' The function allows to simulate data for a truncated dataset with truncation
#' based on the outcome of two decisions and unobservables. The function allows
#' to specify three tests to detect complementarity between the decisions.
#'
#' This function allows to simulate data of the form \deqn{y ~ N(\beta_0 +
#' \beta_1 x_1 + \beta_2 x_2 + \beta_3 x_1 x_2 - \delta_1 {x_1}^2 - \delta_2
#' {x_2}^2, \sigma)}.
#'
#' The parameters \eqn{\beta, \delta, and \sigma} can be set as parameters in
#' the function. The simulated data are truncated based on the outcome \eqn{y}
#' to simulate optimization by decision makers. The rate of survival determines
#' the extent of truncation and thus the strenght of the optimazation.
#'
#' The function allows to specify three different tests for the complementarity
#' of \eqn{x_1} and \eqn{x_2}. The first family of tests ignores the
#' optimization and runs a regression according to the equation above.
#' Traditionally, the quadratic terms are ignored. The augmented method allow
#' for the inclusion of the quadratic terms. The second family of tests assumes
#' optimal choices. As a result,
#'
#' \deqn{2 delta_1 x_1 = \beta_1 + beta_3 x_2}
#'
#' Typically, this matching approach does not estimate the \eqn{\delta_1}
#' parameter. All tests can be done from this function directly. The function
#' returns a single dataset for every test and parameter combination.
#'
#' @param obs A single numeric value or a list with the number of observations
#'   in the simulated dataset.
#' @param rate A single numeric value or a list with the survival rate. This
#'   should be a number between 0 and 1. The closer to 0 the stronger the
#'   pressure to be optimal.
#' @param sd A single numeric value or a list with the standard deviation of
#'   noise term of the performance. A higher \code{sd} means that the decisions
#'   of interest are less important for performance.
#' @param b A vector with four numeric values or a list of vectors. The
#'   contribution to performance of the two decisions of interest. A four
#'   element vector for the intercept, \code{x1, x2, x1:x2}
#' @param d A vector with two numeric values or a list of vectors. The quadratic
#'   effect of the two decisions that captures their the diminishing returns.
#' @param nsim The number of simulations.
#' @param family_method A single string or list of strings with the family and
#'   method of the tests to be run. Defaults to \code{"all"}. Each string should
#'   be of the from "'family'_'method'" where 'family' can be "interaction" or
#'   "match". 'method' is ignored when \code{family = "match"} and can be set to
#'   "traditional" (without quadratic effects) or "augmented" (with quadratic
#'   effects) when \code{family = "interaction"}
#' @param seed The seed can be set for reproducibility.
#' @param mc_cores The number of cores to be used when running different
#'   simulations on different cores. That means that the core unit for
#'   parallelizing is a simulation.
#' @return A dataframe with test statistics and the parameters for each
#'   simulation and test.
#' @export

run_sim <- function(obs = 200,
                    rate = .5,
                    sd = 1,
                    b1 = c(0, 0, 0, 0),
                    b2 = c(1, 1, 1),
                    d = c(1, 1, 1),
                    g1 = c(0, 0, 0),
                    g2 = c(0, 0, 0),
                    sd_eps = c(0, 0, 0),
                    nsim = 100,
                    family_method = "all",
                    seed = 12345,
                    show_progress = F,
                    mc_cores = 1){
  set.seed(seed)
  output_names <- c("method", "coefficient", "stat", "pvalue", "r2")

  ## Parameter cleaning
  make_list <- function(x, len){
    if (!is.list(x)){
      return(list(x))
    } else {
      return(x)
    }
  }
  obs <- make_list(obs)
  rate <- make_list(rate)
  sd <- make_list(sd)
  b1 <- make_list(b1)
  b2 <- make_list(b2)
  d <- make_list(d)
  g1 <- make_list(g1)
  g2 <- make_list(g2)
  sd_eps <- make_list(sd_eps)

  if (family_method == "all"){
    family_method <- list("match", "interaction_traditional",
                          "interaction_augmented", "interaction_moderation",
                          "interaction_moderationaugmented",
                          "interaction_control",
                          "conditional")
  } else if (family_method == "basic"){
    family_method <- list("match", "interaction_traditional")
  } else{
    family_method <- make_list(family_method)
  }

  ## The tests to be performed
  len_fm <- length(family_method)
  familys_in <- vector("list", len_fm)
  methods_in <- vector("list", len_fm)
  for (fm in 1:len_fm){
      family_method_in <- unlist(strsplit(family_method[[fm]], "_",
                                          fixed = TRUE))
      familys_in[[fm]] <- family_method_in[1]
      if (length(family_method_in) == 1){
        methods_in[[fm]] <- "traditional"
      } else if(length(family_method_in) == 2){
        methods_in[[fm]] <- family_method_in[2]
      } else if (length(family_method_in) > 2){
        stop(paste(unlist(params$family_method),
                   "contains too many underscores"))
      }
  }

  ## Variations
  parameters <- expand.grid(obs, rate, sd, b1, b2, d, g1, g2, sd_eps)
  N_variations <- nrow(parameters)
  names(parameters) <- c("obs", "rate", "sd", "b1", "b2",
                         "d", "g1", "g2", "sd_eps")

  ## Results
  N_res <- N_variations * len_fm * nsim
  result <- data.frame(method = character(0),
                       coefficient = numeric(0),
                       stat = numeric(0),
                       pvalue = numeric(0),
                       r2 = numeric(0),
                       b1 = character(0),
                       b2 = character(0),
                       d = character(0),
                       sd = character(0),
                       g1 = character(0),
                       g2 = character(0),
                       sd_eps = character(0),
                       N = numeric(0),
                       survival_rate = numeric(0),
                       sim_id = numeric(0),
                       stringsAsFactors = FALSE)

  for (r in 1:N_variations) {
    cat(paste(r, N_variations, sep = "/"))
    cat("\n")
    if (show_progress == TRUE){
      pb <- txtProgressBar(min = 0, max = nsim)
    }
    params <- lapply(parameters[r,], unlist)
    params_in <- lapply(params, function(x){
      ifelse(length(x) == 1, x, paste(x, collapse = ", "))
    })

    #### parallel ###

    result_list <- parallel::mclapply(1:nsim,
                                      run_single_sim,
                                      params,
                                      params_in,
                                      familys_in,
                                      methods_in,
                                      mc.cores = mc_cores)
    result <- rbind(result, do.call(rbind, result_list))
  }
  return(result)
}

run_single_sim <- function(x, params, params_in, familys, methods){
  result_run <- data.frame()
  dat <- do.call(create_sample, params)

  for (fm in 1:length(familys)){
    reg_list <- run_regression(dat,
                               family = familys[[fm]],
                               method = methods[[fm]])
    new_result <- format_reg(reg_list)
    for (i in 1:nrow(new_result)){
      result_run <- rbind(result_run, cbind(new_result[i,], params_in,
                                            list(sim_id = x)))
    }
  }
  return(result_run)
}


