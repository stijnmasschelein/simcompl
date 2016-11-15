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
#' @param g A vector of two coefficients that indicate the size of the
#' moderation effect of z on the relation between x and y.
#' @return A dataset with simulated \code{y} and \code{x1}, \code{x2} for
#'   surviving observations. The data has dimensions = \code{c(obs, 7)}.
#' @export

create_sample <-
  function(obs = 200,
           rate = .5,
           sd = 1,
           b1 = c(0, 0, 0, 0),
           b2 = c(1, 1, 1),
           d = c(.5, .5, .5),
           g1 = c(0, 0, 0),
           g2 = c(0, 0, 0)
           ){
    check_numeric(obs, 1)
    check_numeric(rate, 1, c(0, 1))
    check_numeric(sd, 1, c(0, Inf))
    check_numeric(b1, 4)
    check_numeric(b2, 3)
    check_numeric(d, 3, c(0, Inf))
    check_numeric(g1, 3)
    check_numeric(g2, 3)

    freq <- 1/rate
    if(ceiling(freq) != floor(freq)){
      cfreq <- c(ceiling(freq), floor(freq))
      probfreq <- c(freq - cfreq[2], cfreq[1] - freq)
      Nfunction <- function(){
        sample(size = 1, x = cfreq, prob = probfreq)
      }
    }else{
      Nfunction <- function(){freq}
    }

    xsurv <- matrix(rep(NA, 3 * obs), ncol = 3)
    zsurv <- rep(NA, obs)
    ysurv <- rep(NA, obs)

    for (i in 1:obs){
      N <- Nfunction()
      z <- replicate(N, rnorm(1))
      w <- replicate(N, rnorm(1))
      x1 <- rnorm(N)
      x2 <- rnorm(N)
      x3 <- rnorm(N)
      x_exp  <- cbind(model.matrix( ~ (x1 + x2 + x3) ^ 2), x1 ^ 2,  x2 ^ 2,
                      x3 ^ 2, z * x1, z * x2, z * x3, w * x1, w * x2, w * x3)
      yhat = x_exp %*% c(b1, b2, -d, g1, g2)
      y = rnorm(yhat, yhat, sd)
      # stupid selection mechanism
      # select the rate highest percentile, can probably be improvemed by working
      # with a probabilitstic model. Only improvement when directly modeling the
      # censoring?
      max_y = max(y)
      ysurv[i] = max_y
      xsurv[i,] = cbind(x1, x2, x3)[y == max_y, ]
      zsurv[i] = z[y == max_y]
    }
    dt <-  as.data.frame(cbind(ysurv, xsurv, zsurv))
    names(dt) <- c("y", "x1", "x2", "x3", "z")
    return(dt)
  }
