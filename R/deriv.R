###
# This file is part of the package "comphy"
###

#' First derivative for an irregular grid
#'
#' Computes the first derivative at a given point, for a function
#' known only through its values tabulated on an irregular grid,
#' where the distance between successive points of the variable
#' is in general not constant. The algorithm is based on the
#' divided differences (see \code{\link{divdif}}).
#' 
#' This numerical derivative should be used only when the 
#' function is known at specific points. When the analytic form
#' of the function is available and the grid of values of the
#' independent variable can be arbitrarily chosen, then it is
#' better to compute the derivative using other more appropriate
#' and faster methods. 
#' 
#' @param x0 A vector of real numbers. These are the values 
#' where the first derivative needs to be calculated. These 
#' values need to be within the interval defined by the 
#' tabulated grid, \code{x}.
#' @param x A vector of real numbers. Grid points corresponding 
#' to the tabulated (known) values of the function.
#' @param f A vector of real numbers. Tabulated (known) values 
#' of the function, corresponding to the grid \code{x}.
#' @return A vector of real numbers. These are the numeric 
#' approximations to the first derivative of the function at
#' all values in \code{x0}. The first derivative is exact when
#' the function is a polynomial of degree \eqn{n}, where
#' \eqn{n} is less than the number of tabulated values.
#'
#' @examples 
#' # Tabulated values: f(x) = 2*x^2-1
#' x <- c(0,1,3,7)
#' f <- 2*x^2-1
#' 
#' # The derivative needs to be computed at three values
#' x0 <- c(1.1,4,6.5)
#' 
#' # First derivatives
#' f1 <- deriv_irr(x0,x,f)
#' print(f1)
#' 
#' @export
deriv_irr <- function(x0,x,f) {
  # Polynomial's degree
  n <- length(x)-1
  
  # Calculate divided differences
  P <- divdif(x,f)
  
  # Derivatives of the approximating polynomial
  ders <- c()
  for (iD in 1:length(x0)) {
  
    # Vector of differences (x0-x[i])
    vdiff <- x0[iD]-x
  
    # Main sum
    extT <- P[1,2]
    if (n >= 2) {
      for (i in 3:(n+1)) {
        extT <- extT+P[1,i]*aprod(vdiff,i)
      }
    }
    ders <- c(ders,extT)
  }
  
  return(ders)
}

#' Forward differences
#'
#' Computes forward differences of all orders up to n, based
#' on n+1 tabulated points on a regular grid.
#' 
#' The forward difference of first order is
#' \deqn{
#'  f(x_i+h)-f(x_i)
#' }
#' Forward differences of higher orders follow from this one,
#' where the function \eqn{f} is replaced by the forward
#' difference of previous order. All values are contained in a
#' \eqn{(n+1)\times(n+1)} upper triangular matrix.
#' 
#' @param f A vector of real numbers. Tabulated (known) values 
#' of the function, corresponding to a regular grid.
#' @return An upper triangular matrix with \eqn{n+1} rows and
#' \eqn{n+1} columns. The first column includes the tabulated 
#' values of the function. The second column includes the eqn{n}
#' forward differences of first order and a zero. The third
#' column includes the \eqn{n-1} forward differences of second
#' order and two zeros. And so on.
#'
#' @examples 
#' # Tabulated values: f(x) = x^3+x^2-x-1
#' x <- c(0,1,2,3)
#' f <- x^3+x^2-x-1
#' 
#' # Triangular matrix with forward differences
#' F <- forwdif(f)
#' print(F)
#' 
#' @export
forwdif <- function(f) {
  # The grid over which f are given is assumed regularly spaced.
  
  # Polynomial degree
  n <- length(f)-1
  
  # Triangular matrix of differences
  F <- matrix(rep(0,times=(n+1)*(n+1)),nrow=n+1,ncol=n+1)
  F[,1] <- f
  
  # Main loop
  for (i in 2:(n+1)) {
    F[1:(n-i+2),i] <- F[2:(n-i+3),i-1]-F[1:(n-i+2),i-1]
  }
  
  return(F)
}


#' Backward differences
#'
#' Computes backward differences of all orders up to n, based
#' on n+1 tabulated points on a regular grid.
#' 
#' The backward difference of first order is
#' \deqn{
#'  f(x_i)-f(x_i-h)
#' }
#' Backward differences of higher orders follow from this one,
#' where the function \eqn{f} is replaced by the backward
#' difference of previous order. All values are contained in a
#' \eqn{(n+1)\times(n+1)} lower triangular matrix.
#' 
#' @param f A vector of real numbers. Tabulated (known) values 
#' of the function, corresponding to a regular grid.
#' @return A lower triangular matrix with \eqn{n+1} rows and
#' \eqn{n+1} columns. The first column includes the tabulated 
#' values of the function. The second column includes a zero
#' and the eqn{n} backward differences of first order. The third
#' column includes two zeros and the \eqn{n-1} forward 
#' differences of second order. And so on.
#'
#' @examples 
#' # Tabulated values: f(x) = x^3+x^2-x-1
#' x <- c(0,1,2,3)
#' f <- x^3+x^2-x-1
#' 
#' # Triangular matrix with backward differences
#' B <- backdif(f)
#' print(B)
#' 
#' @export
backdif <- function(f) {
  # The grid over which f are given is assumed regularly spaced.
  
  # Polynomial degree
  n <- length(f)-1
  
  # Triangular matrix of differences
  B <- matrix(rep(0,times=(n+1)*(n+1)),nrow=n+1,ncol=n+1)
  B[,1] <- f
  
  # Main loop
  for (i in 2:(n+1)) {
    B[i:(n+1),i] <- B[i:(n+1),i-1]-B[(i-1):n,i-1]
  }
  
  return(B)
}


#' First derivative on a regular grid
#'
#' Computes the first derivative of a function at selected points using the 
#' forward difference, backward difference, or centred difference. A regularly
#' spaced grid with corresponding values of the function must be available,
#' as well as a subset of the same grid points at which the derivative must
#' be calculated. For forward and backward differences, the last, respectively 
#' the derivative cannot be calculated at the first or last grid point. For
#' centred difference it cannot be calculated at both first and last grid point. 
#'
#' @param x0 A numeric vector. Values at which the derivative is computed.
#'           Must be an exact subset of \code{x}. Approximated values of
#'           \code{x} will not be accepted.
#' @param x A numeric vector. Regular grid points where the function is tabulated.
#' @param f A numeric vector. Tabulated values of the function at grid \code{x}.
#' @param scheme A one-letter character indicating which difference to use. 
#'               Possible values are "c", "f", "b" for centred, forward and
#'               backward, respectively.
#'
#'
#' @examples
#' x <- seq(0, 1, length.out = 7)
#' f <- x^3 + x^2 - x - 1
#' x0 <- c(0.25, 0.5, 0.75)
#' deriv_reg(x0, x, f)
#'
#' @export
deriv_reg <- function(x0,x,f,scheme="c") {
  # x and f must have same length
  stopifnot(length(x) == length(f))
  
  # The grid x must be regularly spaced.
  # Common practice to use sqrt(eps) rather than eps itself.
  n <- length(x)-1
  h <- x[2]-x[1]
  if (any(abs(diff(x)-h) > .Machine$double.eps^0.5)) 
    stop("Grid 'x' must be regular.")
  
  # Stop unknown values for scheme
  if (scheme != "c" & scheme != "f" & scheme !="b") {
    stop("Wrong letter for scheme to use (c,f,b).")
  }
  
  # Verify requested values in x0 exist in x
  idx <- vapply(x0,function(xi) which.min(abs(x-xi)),integer(1))
  closest_x <- x[idx]
  if (any(abs(closest_x-x0) > .Machine$double.eps^0.5)) {
    stop("One or more values of x0 do not coincide with values of x.")
  }
  
  # Use central difference
  if (scheme == "c") {
    F <- forwdif(f)
    B <- backdif(f)
    
    # Build central differences
    CC <- F[,2]+B[,2]
    
    # Derivatives
    ders <- c(NA,CC[2:n]/(2*h),NA)
  }
  
  # Use forward difference
  if (scheme == "f") {
    F <- forwdif(f)
    
    # Derivatives
    ders <- c(F[1:n,2]/h,NA)
  }
  
  # Use backward difference
  if (scheme == "b") {
    B <- backdif(f)
    
    # Derivatives
    ders <- c(NA,B[2:(n+1),2]/h)
  }
  
  # Only requested values
  ders <- ders[idx]
  
  return(ders)
}


###--------------------------------------------------------------
### Auxiliary functions (not for users)
###--------------------------------------------------------------

# Function to generate the elements for the sum
# in the derivative on irregular grid, function
aprod <- function(vdiff,i) {
  # Number of elements
  m <- i-1
  
  # Combinations without repetitions
  M <- combn(1:m,i-2)
  
  # Sum
  innT <- 0
  for (j in 1:length(M[1,])) {
    innT <- innT+prod(vdiff[M[,j]])
  }
  
  return(innT)
}


# Permutations of a vector
getPerms <- function(x) {
  if (length(x) == 1) {
    return(x)
  }
  else {
    res <- matrix(nrow = 0, ncol = length(x))
    for (i in seq_along(x)) {
      res <- rbind(res, cbind(x[i], Recall(x[-i])))
    }
    return(res)
  }
}