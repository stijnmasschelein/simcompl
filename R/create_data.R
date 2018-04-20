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
  if (!is.null(interval) && (any(input < interval[1]) ||
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
           rep = 1,
           rate = .5,
           sd = 1,
           sd_eps = c(0, 0, 0),
           b1 = c(0, 0, 0, 0),
           b2 = c(1, 1, 1),
           d = c(1, 1, 1),
           g1 = c(0, 0, 0),
           g2 = c(0, 0, 0)
           ){
    check_numeric(obs, 1, c(0, Inf))
    check_numeric(rep, 1, c(0, Inf))
    check_numeric(rate, 1, c(0, 1))
    check_numeric(sd, 1, c(0, Inf))
    check_numeric(sd_eps, 3, c(0, Inf))
    check_numeric(b1, 4)
    check_numeric(b2, 3)
    check_numeric(d, 3, c(0, Inf))
    check_numeric(g1, 3)
    check_numeric(g2, 3)

    # set Nfunction
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

    xsurv <- matrix(rep(NA, 3 * obs * rep), ncol = 3)
    zsurv <- rep(NA, obs * rep)
    wsurv <- rep(NA, obs * rep)
    ysurv <- rep(NA, obs * rep)
    id <- rep(NA, obs * rep)

    for (i in 1:obs){
      N <- Nfunction()
      eps <- c(0, rnorm(3, rep(0, 3), sd_eps)) # 'fixed' effects
      w <- rep(rnorm(1), N) # 'fixed' effects
      for (j in 1:rep){
        z <- rep(rnorm(1), N) # varies for identification
        x1 <- runif(N, min = -5, max = 5)
        x2 <- runif(N, min = -5, max = 5)
        x3 <- runif(N, min = -5, max = 5)
        x_exp  <- cbind(model.matrix( ~ (x1 + x2 + x3) ^ 2), x1 ^ 2,  x2 ^ 2,
                        x3 ^ 2, z * x1, z * x2, z * x3, w * x1, w * x2, w * x3)
        yhat = x_exp %*% c(b1 + eps, b2, -d/2, g1, g2)
        y = rnorm(length(yhat), yhat, sd)
        # stupid selection mechanism
        # can probably be improvemed by working with a probabilitstic model.
        # Only improvement when directly modeling the censoring?
        max_y = max(y)
        count = rep * (i - 1) + j
        ysurv[count] = max_y
        xsurv[count,] = cbind(x1, x2, x3)[y == max_y, ]
        zsurv[count] = z[y == max_y]
        wsurv[count] = w[y == max_y]
        id[count] = i
      }
    }
    dt <-  as.data.frame(cbind(ysurv, xsurv, zsurv, wsurv, id))
    names(dt) <- c("y", "x1", "x2", "x3", "z", "w", "id")
    return(dt)
  }
