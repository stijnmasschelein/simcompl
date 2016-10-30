#' check whether an object is numeric and of a predefined length
#'
#' @param input The object to check
#' @param len The desired length of the object
#' @param interval The interval the inputs should be in


check_numeric <- function(input, len, interval = NULL) {
  name <- deparse(substitute(input))
  if (!is.numeric(input)) {
    stop(paste(name, "is not numerical"))
  }
  if (length(input) != len) {
    stop(paste(name, "should have", len, ifelse(len == 1, "item", "items")))
  }
  if (!is.null(interval) && (any(input <= interval[1]) ||
                             any(input >= interval[2]))) {
    stop(paste(name, "does not lie in the interval"))
  }
}

#' Create a sample of surviving observations for simulation
#'
#' @param obs The number of observations in the dataset
#' @param rate The survival rate. This should be a number between 0 and 1. The
#'   closer to 0 the stronger the pressure to be optimal.
#' @param sd The standard deviation of noise term of the performance. A higher
#'   \code{sd} means that the decisions of interest are less important for
#'   performance.
#' @param b The contribution to performance of the two decisions of interest. A
#'   four element vector for the intercept, \code{x1, x2, x1:x2}
#' @param d The quadratic effect of the two decisions that captures their the
#'   diminishing returns.
#' @return A dataset with simulated \code{y} and \code{x1}, \code{x2} for
#'   surviving observations. The data has dimensions = \code{c(obs, 7)}.

create_sample <-
  function(obs = 200,
           rate = .5,
           sd = 1,
           b = c(0, 1, 1, 1)
           ,
           d = c(.5, .5)) {
    check_numeric(obs, 1)
    check_numeric(rate, 1)
    check_numeric(sd, 1, c(0, Inf))
    check_numeric(b, 4)
    check_numeric(d, 2, c(0, Inf))

    N <- obs / rate
    # For now the x's are independently generated
    x1 <- rnorm(N)
    x2 <- rnorm(N)
    x_exp  <- cbind(model.matrix( ~ (x1 + x2) ^ 2), x1 ^ 2,  x2 ^ 2)
    yhat = x_exp %*% c(b,-d)
    y = rnorm(yhat, yhat, sd)

    # stupid selection mechanism
    # select the rate highest percentile, can probably be improvemed by working
    # with a probabilitstic model. Only improvement when directly modeling the
    # censoring?
    survive <- I(y > quantile(y, (1 - rate)))
    xsurv <- cbind(x1[survive], x2[survive])
    ysurv <- y[survive]
    dt <-  as.data.frame(cbind(ysurv, xsurv))
    names(dt) <- c("y", "x1", "x2")
    return(dt)
  }
