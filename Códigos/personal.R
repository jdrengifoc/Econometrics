library(geometry) # Dot product.
library(tseries) # Augmented Dickey Fuller.
library(urca)
# AR(p) -------------------------------------------------------------------
#'AR(p) with standard error series with the parametrization and traditional transformed
#' series.
#' 
#' @param p Order of the AR. Must be a integer number.
#' @param n Sample size of the series.
#' @param phi Real coefficients of the AR. Predefined, are random numbers
#' between -1 and 1. Must be of length p.
#' @param x0 p initial values of the series. Predefined, are random numbers
#' between -1 and 1. Must be of length p.
#' @param alpha Tendency real parameter. Real number.
#' @return List that has the original AR(\code{p}) series with parametrization
#' \code{phi} and \code{alpha}, as well as an associated series. This series are:
#' x.tend (tendency stationary), x.int (differentiation stationary) and x.int.tend
#' (tendency and differentiation stationary)
#' @example 
#' # AR(1) with 1000 sample sizes.
#' p <- 1; n <- 1000
#' phi <- c(0.2, 0.4); alpha <- 0.03
#' x0 <- c(0.3, 0.5)
#' ar <- AR(p, n, phi, x0, alpha)
#' plot(ar$x)
#' plot(ar$x.int)
#' @example 
#' # AR(4) with 1000 sample sizes and default parameters.
#' ar <- AR(4, 1000)
#' plot(ar$x)
#' plot(ar$x.tend)
AR <- function(p, n, 
               phi = runif(p, min = -1, max = 1),
               x0 = runif(p, min = -1, max = 1),
               alpha = 0.01, delta = 0){
  # Initiallize series.
  x <- c(x0, rep(NA, n-p))
  
  # AR simulation for each observation.
  for (i in (p+1):n){
    I <- (i-p):(i-1)
    if (p==1){
      x[i] <- delta + x[I]*phi + rnorm(1)
    } else if (p > 1){
      x[i] <- delta + dot(x[I], phi) + rnorm(1)
    } else {
      warning("p must be a positive integer.")
    }
  }
  # Organize results.
  results <- list(x = x,
                  phi = phi,
                  alpha = alpha,
                  x.int = diffinv(x),
                  x.tend = alpha*1:n + x,
                  x.int.tend = alpha*1:(n+1) + diffinv(x)
                  )
  return(results)
}

# Game zone ---------------------------------------------------------------


# Parameters
p <- 2; n <- 1000
# Default
phi <- c(0.2, 0.4); alpha <- 0.01
x0 <- c(0.3, 0.5)
ar <- AR(p, n, delta = 2)

# Default R x <- arima.sim(list(order = c(1,0,0), ar = 0.2),n = 100)

# Estacionario.
plot(ar$x, type = 'l')
# Estacionario en tendencia.
plot(ar$x.tend, type = 'l')
# Estacionario en diferencia.
plot(ar$x.int, type = 'l')
# Estacionario en diferencia y en tendencia.
plot(ar$x.int.tend, type = 'l')
acf(ar$x)
pacf(ar$x)
adf.test(x)


# Pruebas de raices unitarias ---------------------------------------------

#' Augmented dickey Fuller
#' Standard: AR(1) k/lags = 0
#' Always use trend and intercept.
adf.test(ar$x.tend + 20, k = 0)

# Supports trend and intercept non by default.
summary(ur.df(ar$x + 20, lags = 0, selectlags = "BIC"))
summary(ur.df(ar$x + 20, type = "trend", lags = 0, selectlags = "BIC"))
  
#' Phillips & Perron
#' lshort is TRUE, then the truncation lag parameter is set to 
#' trunc(4*(n/100)^0.25), otherwise trunc(12*(n/100)^0.25) is used. 
pp.test(ar$x, lshort = T)


#' KPSS
#' lshort is TRUE, then the truncation lag parameter is set to 
#' trunc(4*(n/100)^0.25), otherwise trunc(12*(n/100)^0.25) is used.
kpss.test(ar$x, null = "Level", lshort = T)
#' null supports trend
kpss.test(ar$x, null = "Trend", lshort = F)
kpss.test(ar$x.tend, null = "Level", lshort = F)
kpss.test(ar$x.tend, null = "Trend", lshort = F)

# Also allows user lags.
summary(ur.kpss(ar$x.tend, type="tau", use.lag = 1))

## Para la presentación

ar <- AR(p = 2, n = 1000, x0 = c(14, 14.1), phi = c(0.2, -0.6), alpha = 0.02, delta = 20)
data <- ar$x.tend
plot(data, type = 'l')
plot(ar$x, type = 'l')
adf.test(data)
pp.test(data)
kpss.test(data)
kpss.test(data, null = "Trend")


# In Development ----------------------------------------------------------

# Un modelo ARMA(p,q) es de la forma
# xt = ϕ1xt−1 + ϕ2xt−2 + ⋯ + ϕpxt-p + εt − θ1εt−1 − θ2εt−2 − ⋯ − θqεt−q
# phi(L)x_t = tetha(L)eps_t

pp <- function(p){
  return(exists(p))
}

ARMA <- function(p, q, n, params = list(phi = runif(p, min = -1, max = 1), 
                                        tetha = runif(q, min = -1. max = 1)),
                 x0 = runif(p+q, min = -1, max = 1),
                 alpha = 0.01, delta = 0){
  # Initiallize series.
  x <- c(x0, rep(NA, n-p))
  
  # AR simulation for each observation.
  for (i in (p+1):n){
    I <- (i-p):(i-1)
    if (p==1){
      x[i] <- delta + x[I]*phi + rnorm(1)
    } else if (p > 1){
      x[i] <- delta + dot(x[I], phi) + rnorm(1)
    } else {
      warning("p must be a positive integer.")
    }
  }
  # Organize results.
  results <- list(x = x,
                  phi = phi,
                  alpha = alpha,
                  x.int = diffinv(x),
                  x.tend = alpha*1:n + x,
                  x.int.tend = alpha*1:(n+1) + diffinv(x)
  )
  return(results)
}
