nr_columns <- 4

#' Perform the traditional single equation regression approach under the
#' assumption that \code{x1, x2} are optimally chosen.
#'
#' @export
#' @family tests
#' @param dat A simulated dataset with 5 columns.
#' @seealso interaction_reg run_regression format_reg
#' @return A list with the regression object and the name of the method

match_reg <- function(dat){
  stopifnot(ncol(dat) >= nr_columns)
  reg <- lm(x1 ~ x2 + z, data = dat)
  return(list(reg = reg, method = "matching"))
}

#' Perform the traditional single equation regression approach under the
#' assumption that \code{x1, x2} are randomly chosen.
#'
#' @export
#' @family tests
#' @param dat A simulated dataset with 5 columns.
#' @param method 'traditional' for the typical model with interaction and
#'   'augmented' for the model with interaction and quadratic effects.
#'   'moderation' takes into account the moderation effect of z1 on y. The
#'   'moderationaugmented` option takes into account both moderation and
#'   'augmented' effects.
#' @seealso match_reg run_regression format_reg
#'
#' @return A list with the regression object and the name of the method

interaction_reg <- function(dat, method = "traditional"){
  stopifnot(ncol(dat) >= nr_columns)

  method_options <- c("traditional", "augmented", "moderation",
                      "moderationaugmented")
  if (!(method %in% method_options)){
    method_options_str <- paste(method_options, collapse = ", ")
    stop("method has no valid input. It should be one of the following: ",
         method_options_str, ".")
  } else {
  reg <- switch(method,
                traditional = lm(y ~ x1*x2, data = dat),
                augmented = lm(y ~ x1*x2 + I(x1^2) + I(x2^2), data = dat),
                moderation = lm(y ~ x1*x2 + x1*z + x2*z, data = dat),
                moderationaugmented = lm(y ~ x1*x2 + x1*z + x2*z +
                                            I(x1^2) + I(x2^2), data = dat))
  return(list(reg = reg, method = paste("interaction", method, sep = "_")))
  }
}

#' The SUR approach to detecting the complementarity
#'
#' @export
#' @family tests
#' @param dat A simulated dataset with 5 columns
#' @seealso match_reg run_regression format_reg
#'
#' @return A list with the sur object and the name of the method

sur_reg <- function(dat){
  stopifnot(ncol(dat) == nr_columns + 1)
  # should be rewritten. See packages by Hadley Wickham
  if (!require("systemfit")){
    stop("The `systemfit` is not installed. Install `systemfit` to use the sur
         option")
  }
  reg <- systemfit::systemfit(list(eq1 = x1 ~ z, eq2 = x2 ~ z),
                              data = dat, method = "SUR", maxit = 100)
  return(list(reg = reg, method = "sur"))
}

#' The conditional correlation approach to detecting the complementarity
#'
#' @export
#' @family tests
#' @param dat A simulated dataset with 5 columns
#' @seealso match_reg run_regression format_reg
#'
#' @return A list with the sur object and the name of the method

cond_reg <- function(dat){
  reg1 <- lm(x1 ~ z, data = dat)
  reg2 <- lm(x2 ~ z, data = dat)
  return(list(reg1 = reg1, reg2 = reg2, method = "conditional"))
}

#' A general function to run all possible regression based approaches.
#'
#' @export
#' @family tests
#' @param dat A simulated dataset with 7 columns.
#' @param family "match" for the approach where the optimal choice of the
#'    \code{x}'s is chosen and "interaction" where the \code{x}'s are assumed
#'    to be chosen optimally.
#' @param method "traditional" for the typical model with interaction,
#'   'augmented' for the model with interaction and quadratic effects, and
#'   'sur' for the seemingly unrelated regression method.
#'   \code{method} is ignored when \code{family = "match"} or
#'   \code{family = "sur"}.
#' @seealso match_reg interaction_regression format_reg
#'
#' @return A list with the regression object and the name of the method

run_regression <- function(dat, family = "match", method = "traditional"){
  family_options <- c("match", "interaction", "sur")
  if (!(family %in% family_options)){
    family_options_str <- paste(family, collapse = ", ")
    stop("family has no valid input. It should be one of the following: ",
       family_options_str, ".")
  } else {
    reg_list <- switch(family,
                       match = match_reg(dat),
                       interaction = interaction_reg(dat, method = method),
                       sur = sur_reg(dat))
    return(reg_list)
  }
}

#' Transform the regression and method list into a vector.
#'
#' @export
#' @param list_reg A list with two elements. The regression object and a string
#' with the name of the methodology.
#' @seealso match_reg interaction_regression run_regression
#' @return A row vector with type of analysis, coefficient, t - statistic, and
#' p - value of the complementarity test for the 'sur' (in one estimation) and
#' for the two ols regression approach with the correlation between residuals
#' ('resid')

format_reg <- function(list_reg){
  if (is.null(list_reg)){
    stop("The argument to format_reg is NULL")
  } else if(list_reg$method == "sur"){
    reg <- list_reg$reg
    cov_matrix <- reg$residCov
    cor_resid <- cov_matrix[1,2]/sqrt(cov_matrix[1,1] * cov_matrix[2,2])
    summ_reg <- summary(reg)
    cortest <- cor.test(summ_reg$residuals[,1], summ_reg$residuals[,2])
    results <- data.frame(method = "sur",
                          coefficient = cor_resid,
                          r2 = summ_reg$ols.r.squared,
                          pvalue =  cortest$p.value,
                          stat = cortest$statistic,
                          stringsAsFactors = FALSE)
  } else if(list_reg$method == "conditional"){
    reg1 <- list_reg$reg1
    reg2 <- list_reg$reg2
    summ1 <- summary(reg1)
    summ2 <- summary(reg2)
    cor_test <- cor.test(resid(reg1), resid(reg2))
    results <- data.frame(method = "conditional",
                          coefficient = cor_test$statistic,
                          r2 = mean(summ1$r.squared, summ2$r.squared),
                          pvalue = cortest$p.value,
                          stat = cortest$statistic,
                          stringsAsFactors = FALSE)
  }
    else {
    reg <- list_reg$reg
    method <- list_reg$method
    family <- strsplit(method, split = "_")[[1]][1]
    summ <- summary(reg)
    coefstring <- switch(family,
                   interaction = "x1:x2",
                   matching = "x2")
    results <- data.frame(method = method,
                 coefficient = reg$coefficients[coefstring],
                 stat = summ$coefficients[coefstring, "t value"],
                 pvalue = summ$coefficients[coefstring, "Pr(>|t|)"],
                 r2 = summ$r.squared,
                 stringsAsFactors = FALSE)
  }
  row.names(results) <- NULL
  return(results)
}
