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
#' @return A dataframe with test statistics and the parameters for each
#'   simulation and test.


run_sim <- function(obs = 200,
                    rate = .5,
                    sd = 1,
                    b = c(0, 1, 1, 1),
                    d = c(0.5, 0.5),
                    nsim = 100,
                    family_method = "all",
                    seed = 12345,
                    show_progress = F){
  set.seed(seed)
  output_names <- c("method", "coefficient", "tstat", "pvalue", "r2")

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
  b <- make_list(b)
  d <- make_list(d)

  parameters <- expand.grid(obs, rate, sd, b, d)
  N_variations <- nrow(parameters)
  names(parameters) <- c("N", "survival_rate", "sd", "b", "d")
  if (family_method == "all"){
    family_method <- list("match", "interaction_traditional",
                          "interaction_augmented")
  } else{
    family_method <- make_list(family_method)
  }

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

  N_res <- N_variations * len_fm * nsim
  result <- data.frame(method = rep("", N_res),
                       coefficient = rep(NA, N_res),
                       tstat = rep(NA, N_res),
                       pvalue = rep(NA, N_res),
                       r2 = rep(NA, N_res),
                       b = rep("", N_res),
                       d = rep("", N_res),
                       sd = rep(NA, N_res),
                       N = rep(NA, N_res),
                       survival_rate = rep(NA, N_res),
                       sim_id = rep(NA, N_res),
                       stringsAsFactors = FALSE)

  for (r in 1:N_variations) {
    cat(paste(r, N_variations, sep = "/"))
    cat("\n")

    if (show_progress == TRUE){
      pb <- txtProgressBar(min = 0, max = nsim)
    }
    params <- parameters[r,]
    b_in <- unlist(params$b)
    d_in <- unlist(params$d)
    sd_in <- unlist(params$sd)
    obs_in <- unlist(params$N)
    rate_in <- unlist(params$survival_rate)
    params_in <- c(paste(b_in, collapse = ", "),
                       paste(d_in, collapse = ", "),
                       sd_in, obs_in, rate_in)

    for (s in 1:nsim){
      dat <- create_sample(
        b = b_in, d = d_in, sd = sd_in, obs = obs_in, rate = rate_in
      )
      for (fm in 1:length(familys_in)){
        reg_list <- run_regression(dat, familys_in[[fm]],
                                   method = methods_in[[fm]])
        index_res <- len_fm * (N_variations * (s-1) + (r-1)) + fm
        result[index_res, ] <- c(format_reg(reg_list),
                                               params_in, s)
      }
      if (show_progress == TRUE){
        setTxtProgressBar(pb, s)
      }
    }
    if (show_progress == TRUE){
    close(pb)
    }
  }
  vars <- c("sd", "N", "survival_rate")
  result[vars] <- lapply(result[vars], as.numeric)
  return(result)
}
