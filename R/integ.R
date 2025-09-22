###
# This file is part of the package "comphy"
###

#' Numerical integration using the trapezoid or simpson's rule
#'
#' Computes the definite integral of \eqn{f(x)} between \eqn{a}
#' and \eqn{b}, using one of the three numerical integration
#' Newton-Cotes rules, trapezoid, Simpson's 1/3 or Simpson's 3/8.
#' 
#' The default method is Simpson's 1/3 rule. For this method to
#' be applied correctly, the number of regular intervals must be
#' even. If this is not the case, the area corresponding to the
#' last interval will be calculated with the trapezoid rule with
#' a warning being posted.
#' 
#' When using the Simpson's 3/8 rule, the number of regular
#' intervals must be a multiple of 3. If this is not the case,
#' the last or last two intervals will be computed with the
#' trapezoid rule.
#' 
#' @param x A vector of real numbers. Grid points corresponding 
#' to the tabulated (known) values of the function. These must be
#' equally spaced (regular grid).
#' @param f A vector of real numbers. Tabulated (known) values 
#' of the function, corresponding to the grid \code{x}.
#' @param scheme A character indicating the integration rule to
#'        follow. Possible values are "trap" (trapezoid rule),
#'        "sim13" (Simpson's 1/3 rule), and "sim38" (Simpson's
#'        3/8 rule). Default scheme is "sim13".
#' @return A real number, corresponding to the numeric approximation
#' of the definite integral of \eqn{f(x)}.
#'
#' @examples 
#' # Tabulated values: f(x) = x^2
#' x <- seq(0,2,length.out=21) # number of intervals is even
#' f <- x^2
#' 
#' # Integral between 0 and 2
#' # The correct result is 2^3/3=8/3=2.6666...
#' nvalue <- numint_reg(x,f) # Defaul nethod simpson's 1/3
#' print(nvalue)
#' 
#' # If the number of intervals is not even,
#' # a warning is issued
#' y <- seq(0,2,length.out=22)
#' g <- y^2
#' nvalue <- numint_reg(y,g)
#' print(nvalue)
#' 
#' @export
numint_reg <- function(x,f,scheme="sim13") {
  # Basic checks
  if (length(x) != length(f)) 
    stop("x and f must have the same length.")
  
  # Check grid regularity
  h <- x[2] - x[1]
  if (any(abs(diff(x) - h) > .Machine$double.eps^0.5)) 
    stop("Grid 'x' must be regular.")
  
  # Check scheme
  if (!(scheme %in% c("trap", "sim13", "sim38"))) {
    stop("Scheme must be one of 'trap', 'sim13', or 'sim38'.")
  }
  
  # n is the number of elementary intervals
  # n = n. of points - 1
  n <- length(x) - 1
  
  # Trapezoidal rule
  if (scheme == "trap") {
    # Weights
    w <- c(1,rep(2,length.out=n-1),1)
    
    # Convolution and sum
    tmp <- w*f
    nvalue <- (h/2)*sum(tmp)
  }
  
  # Simpson's 1/3 rule (requires even number of intervals)
  else if (scheme == "sim13") {
    # Know what to do if n is not even
    if (n %% 2 == 0) {
      # Weights
      w <- c(1,rep(c(4,2),length.out=n-2),4,1)
      
      # Convolution and sum
      tmp <- w*f
      nvalue <- (h/3)*sum(tmp)
    } else {
      # Apply Simpson's 1/3 to everything but the last interval
      m <- n - 1
      
      # Weights
      w <- c(1,rep(c(4,2),length.out=m-2),4,1)
      
      # Convolution and sum
      tmp <- w*f[1:n]
      nvalue <- (h/3)*sum(tmp)
      
      # Add last contribution with trapezoid rule
      nvalue <- nvalue + (h/2)*(f[n]+f[n+1])
      
      # Warning as not pure Simpson
      warning("Last contribution to integral is from trapezoid rule as the number of intervals is not even.")
    }
  }
  
  # Simpson's 3/8 rule (requires n divisible by 3)
  else if (scheme == "sim38") {
    # Know what to do if n is not divisible by 3
    if (n %% 3 == 0) {
      # Weights
      w <- c(1,rep(c(3,3,2),length.out=n-3),3,3,1)
      
      # Convolution and sum
      tmp <- w*f
      nvalue <- (3*h/8)*sum(tmp)
    } else {
      # Apply Simpson's 3/8 to everything but the last interval
      if ((n+1)%%3 == 0) {
        m <- n - 2
        
        # Weights
        w <- c(1,rep(c(3,3,2),length.out=m-3),3,3,1)
        
        # Convolution and sum
        tmp <- w*f[1:(m+1)]
        nvalue <- (3*h/8)*sum(tmp)
        
        # Add last contributions with trapezoid rule
        nvalue <- nvalue + (h/2)*(f[n-1] + 2*f[n] + f[n+1]) 
      }
      if ((n+2)%%3 == 0) {
        m <- n - 1
        
        # Weights
        w <- c(1,rep(c(3,3,2),length.out=m-3),3,3,1)
        
        # Convolution and sum
        tmp <- w*f[1:(m+1)]
        nvalue <- (3*h/8)*sum(tmp)
        
        # Add last contributions with trapezoid rule
        nvalue <- nvalue + (h/2)*(f[n] + f[n+1])
      }
      
      # Warning as not pure Simpson 3/8
      warning("Last contribution to integral is from trapezoid rule as the number of intervals is not a multiple of 3.")
    }
  }
  
  return(nvalue)
}


#' Numerical integration using \eqn{n}-point Gaussian quadrature.
#'
#' Computes the definite integral of \eqn{f(x)} between \eqn{a}
#' and \eqn{b}, using the method of Gaussian quadrature. The
#' default number of points is
#'
#' @param f A function to integrate.
#' @param a Lower bound of integration (default -1).
#' @param b Upper bound of integration (default 1).
#' @param n Number of quadrature points (default 5).
#'
#' @return A list with three elements. The first
#'         is a numeric vector containing the nodes
#'         of the quadrature. The second is a numeric
#'         vector containing the corresponding weight.
#'         The third is a real number corresponding to
#'         the approximate value of the integral.
#'         
#' @examples
#' # Integral in [-1,1] of 2x-1.
#' # Value is -2 and n=1 is enough for exact result
#' 
#' # Define the function
#' f <- function(x) {ff <- 2*x-1; return(ff)}
#' 
#' # 1-point quadrature
#' ltmp <- Gquad(f,-1,1,n=1)
#' 
#' # The only zero is x1=0
#' print(ltmp$xt)
#' 
#' # The only weight is w1=2
#' print(ltmp$wt)
#' 
#' # Quadrature gives exact integral
#' print(ltmp$itg)
#' 
#' # 2-point quadrature
#' ltmp <- Gquad(f,-1,1,n=2)
#' print(ltmp) # Same result but more zeros and weights
#' 
#' # Default, n=5, is accurate enough
#' ltmp <- Gquad(exp,-1,1)
#' print(ltmp$itg)
#' 
#' # Different extremes of integration
#' ltmp <- Gquad(exp,1,4)
#' print(ltmp$itg)
#' 
#' @export
Gquad <- function(f,a=-1,b=1,n=5) {
  # The lowest quadrature is linear
  if (n < 1) stop("n must be a positive integer.")
  
  # n=1 must be treated separately
  if (n == 1) {
    x <- 0
    w <- 2
    xt <- (b-a)/2*x+(b+a)/2
    wt <- (b-a)/2*w
    
    # Compute weighted sum
    itg <- sum(wt*f(xt))
    
    # List containing zeros, weights, and integral
    ltmp <- list(xt=xt,wt=wt,itg=itg)
    
    return(ltmp)
  }
  
  # Golub-Welsch algorithm to avoid using a Legendre R package
  i <- 1:(n-1)
  a_diag <- rep(0,length.out=n)          # Diagonal entries (Legendre: all 0)
  b_sub <- i/sqrt(4*i^2-1)               # Sub-diagonal entries
  T <- diag(a_diag)                      # Build Jacobi matrix
  T[cbind(i+1,i)] <- b_sub               # Lower diagonal
  T[cbind(i,i+1)] <- b_sub               # Upper diagonal
  
  eig <- eigen(T,symmetric=TRUE)
  x <- eig$values                        # Nodes in [-1,1]
  V <- eig$vectors
  w <- 2 * (V[1,])^2                     # Weights
  
  # Transform from [-1,1] to [a,b]
  xt <- (b-a)/2*x+(b+a)/2
  wt <- (b-a)/2*w
  
  # Compute weighted sum
  itg <- sum(wt*f(xt))
  
  # List containing zeros, weights, and integral
  ltmp <- list(xt=xt,wt=wt,itg=itg)
  
  return(ltmp)
}



###--------------------------------------------------------------
### Auxiliary functions (not for users)
###--------------------------------------------------------------

# Add here ...
