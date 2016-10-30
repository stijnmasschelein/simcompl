nr_columns <- 3

#' Perform the traditional single equation regression approach under the
#' assumption that \code{x1, x2} are optimally chosen.
#'
#' @family tests
#' @param dat A simulated dataset with 7 columns.
#' @return A list with the regression object and the name of the method

match_reg <- function(dat){
  stopifnot(ncol(dat) == nr_columns)
  reg <- lm(x1 ~ x2, data = dat)
  return(list(reg = reg, method = "matching"))
}

#' Perform the traditional single equation regression approach under the
#' assumption that \code{x1, x2} are randomly chosen.
#'
#' @family tests
#' @param dat A simulated dataset with 7 columns.
#' @param method 'traditional' for the typical model with interaction and
#'   'augmented' for the model with interaction and quadratic effects.
#' @seealso match_reg run_regression format_reg
#'
#' @return A list with the regression object and the name of the method

interaction_reg <- function(dat, method = "traditional"){
  stopifnot(ncol(dat) == nr_columns)

  method_options <- c("traditional", "augmented")
  if (!(method %in% method_options)){
    method_options_str <- paste(method_options, collapse = ", ")
    stop("method has no valid input. It should be one of the following: ",
         method_options_str, ".")
  } else {
  reg <- switch(method,
                traditional = lm(y ~ x1*x2, data = dat),
                augmented = lm(y ~ x1*x2 + I(x1^2) + I(x2^2), data = dat))
  return(list(reg = reg, method = paste("interaction", method, sep = "_")))
  }
}

#' A general function to run all possible regression based approaches.
#'
#' @family tests
#' @param dat A simulated dataset with 7 columns.
#' @param family "match" for the approach where the optimal choice of the
#'    \code{x}'s is chosen and "interaction" where the \code{x}'s are assumed
#'    to be chosen optimally.
#' @param method "traditional" for the typical model with interaction and
#'   'augmented' for the model with interaction and quadratic effects.
#'   \code{method} is ignored when \code{family = "match"}.
#'
#' @return A list with the regression object and the name of the method

run_regression <- function(dat, family = "match", method = "traditional"){
  family_options <- c("match", "interaction")
  if (!(family %in% family_options)){
    family_options_str <- paste(family, collapse = ", ")
    stop("family has no valid input. It should be one of the following: ",
       family_options_str, ".")
  } else {
    reg_list <- switch(family,
                       match = match_reg(dat),
                       interaction = interaction_reg(dat, method = method))
    return(reg_list)
  }
}

#' Transform the regression and method list into a vector.
#'
#' @param list_reg A list with two elements. The regression object and a string
#' with the name of the methodology.
#' @return A row vector with type of analysis, coefficient, t - statistic, and
#' p - value of the complementarity test.

format_reg <- function(list_reg){
  if (is.null(list_reg)){
    stop("The argument to format_reg is NULL")
  } else {
    reg <- list_reg$reg
    method <- list_reg$method
    summ <- summary(reg)
    row_test <- length(reg$coefficients)
    results <- data.frame(method = method,
                 coefficient = reg$coefficients[row_test],
                 tstat = summ$coefficients[row_test, "t value"],
                 pvalue = summ$coefficients[row_test, "Pr(>|t|)"],
                 r2 = summ$r.squared,
                 stringsAsFactors = FALSE)
    row.names(results) <- NULL
    return(results)
  }
}
