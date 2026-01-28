###
# This file is part of the package "comphy"
###


#' Multilinear Least Squares
#'
#' Find the parameters, \eqn{a_1,\dots,a_{m+1}}, of the linear
#' model with \eqn{m} parameters, using the least squares technique
#' on a group of \eqn{n} data points in the \eqn{m} dimensional
#' Cartesian space.
#' 
#' The linear model with \eqn{m} parameters has the following analytic
#' form:
#' \deqn{
#'  y = a_1x_1 + a_2x_2 + \dots + a_mx_m + a_{m+1} 
#' }
#' The \eqn{n} data points are contained in a matrix or data frame with
#' \eqn{m+1} columns and \eqn{n} rows. The first \eqn{m} elements of
#' each row contain the coordinates of a data point; the last element 
#' contains the corresponding value of the linear fitting, \eqn{y_i}.
#' The least squares procedure is carried out as solution of a matrix
#' equation, via the \code{\link[base]{solve}} function.
#' 
#' @param x A \eqn{n \times (m+1)} matrix or data frame where the first 
#'          \eqn{m} elements of each row contain the coordinates of a
#'          data point, and the last element contain the value
#'          corresponding to the linear model.
#' @param intercept A logical variable. It indicates whether to omit
#'                  or keep the constant \eqn{a_{m+1}} in the model.
#'                  The default is \code{intercept=TRUE} as removing
#'                  the constant is not advisable, unless there is
#'                  an absolute certainty (for example in mechanistic
#'                  models) that it has to be removed.
#' @param tol A real number. The solution of a linear system can be
#'            compromised when the condition number of the matrix of
#'            coefficients is particularly high (ill-conditioned
#'            matrices). \code{tol} is the reciprocal of the 
#'            condition number. For values of tol smaller than 1e-17,
#'            ill-conditioning is deemed to be sever enough not to
#'            guarantee an accurate solution. For such values the
#'            function stops execution, returning an error message.
#'            In fact, the solution can still be accurate,
#'            notwithstanding ill-conditioning, and the user can 
#'            force the calculation of a solution using a value of
#'            \code{tol} smaller than 1e-17. Default is \code{NULL},
#'            corresponding to a \code{tol=1e-17}.
#' @return A vector of length \eqn{m} containing the \eqn{m}
#'         numeric values of the estimated linear model parameters.
#'         The function also prints out the numerical value of the
#'         sum of squared residuals. If more than one solution is
#'         possible (infinite-solutions case) the function returns
#'         a \code{NULL} and prints out a related message.
#'         
#' @examples 
#' # 5 points exactly on y = 2x_1 - x_2 + 3
#' p1 <- c(0,1,2)
#' p2 <- c(1,0,5)
#' p3 <- c(1,1,4)
#' p4 <- c(0,2,1)
#' p5 <- c(2,0,7)
#' 
#' # Assemble points in a single matrix for input
#' x <- matrix(c(p1,p2,p3,p4,p5),ncol=3,byrow=TRUE)
#' 
#' # Find the least squares estimate of a_1,a_2,a_3
#' a <- solveLS(x)
#' print(a)
#' 
#' @export
solveLS <- function(x,intercept=TRUE,tol=NULL) {
  # Check input
  ans <- (is.matrix(x) | is.data.frame(x))
  if (!ans) {
    msg <- "Input has to be either a matrix or a data frame\n"
    warning(msg)
    return(NULL)
  }
  
  # Number of data points (n) and parameters (m+1)
  tmp <- dim(x)
  n <- tmp[1]
  m <- tmp[2]-1
  
  # Build A and y matrices
  ones <- matrix(rep(1,times=n),ncol=1)
  A <- as.matrix(x[,1:m])
  colnames(A) <- NULL
  rownames(A) <- NULL
  if (intercept) A <- cbind(A,ones)
  y <- matrix(x[,m+1],ncol=1)
  colnames(y) <- NULL
  rownames(y) <- NULL
  
  # Build F and g matrices to get solution through solve
  F <- t(A) %*% A
  g <- t(A) %*% y
  
  # Check whether F is singular
  d <- det(F)
  if (abs(d) < 1e-6) {
    msg <- paste0("There are infinite solutions to ",
                  "this least squares fitting.\n")
    warning(msg)
    return(NULL)
  }
  
  # Change tolerance, if required
  if (is.null(tol)) tol <- 1e-17
  
  # Solution
  a <- solve(F,g,tol=tol)
  
  # Print out the sum of squared residuals:
  eps <- y-A %*% a
  d <- sum(eps^2)
  message(sprintf("Sum of squared residuals: %f", d))
  
  # Reshape for output
  a <- as.vector(a)
  
  return(a)
}

#' Polynomial Least Squares
#'
#' Find the parameters, \eqn{a_1,\dots,a_{m+1}}, of the polynomial
#' model of degree \eqn{m} (1D function), using the least squares 
#' technique on a group of \eqn{n} data points.
#' 
#' The polynomial model has the following analytic form:
#' \deqn{
#'  y = a_1x^m + a_2x^{m-1} + \dots + a_mx + a_{m+1} 
#' }
#' The \eqn{n} data points are contained in a matrix or data frame 
#' with \eqn{2} columns, containing the coordinates of each data
#' point, and \eqn{n} rows. The least squares procedure is carried 
#' out as solution of a matrix equation, via the 
#' \code{\link{solveLS}} function.
#' 
#' @param pts A \eqn{n \times 2} matrix or data frame where each row
#'            contains the coordinates of a data point used for
#'            regression.
#' @param m An integer. The degree of the polynomial to be used as
#'          model for the regression.
#' @param tol A real number. The solution of a linear system can be
#'            compromised when the condition number of the matrix of
#'            coefficients is particularly high (ill-conditioned
#'            matrices). \code{tol} is the reciprocal of the 
#'            condition number. For values of tol smaller than 1e-17,
#'            ill-conditioning is deemed to be sever enough not to
#'            guarantee an accurate solution. For such values the
#'            function stops execution, returning an error message.
#'            In fact, the solution can still be accurate,
#'            notwithstanding ill-conditioning, and the user can 
#'            force the calculation of a solution using a value of
#'            \code{tol} smaller than 1e-17. Default is \code{NULL},
#'            corresponding to a \code{tol=1e-17}.
#' @return A named list with two elements:
#'  \describe{
#'   \item{a}{A vector of length \eqn{m} containing the \eqn{m}
#'            numeric values of the estimated polynomial's
#'            coefficients.If more than one solution is possible,
#'            (infinite-solutions case) the function returns a
#'            \code{NULL} and prints out a related message.} 
#'  \item{SSE}{A real number. The numerical value of the sum of 
#'             squared residuals.}
#' } 
#'         
#' @examples 
#' # 21 points close to the quadratic x^2 - 5*x + 6
#' x <- seq(-2,5,length=21)
#' set.seed(7766)
#' eps <- rnorm(21,mean=0,sd=0.5)
#' y <- x^2-5*x+6+eps
#' 
#' # Data frame
#' pts <- data.frame(x=x,y=y)
#' 
#' # Regression
#' ltmp <- polysolveLS(pts,m=2)
#' print(names(ltmp))
#' print(ltmp$a)
#' print(ltmp$SSE)
#' 
#' @export
polysolveLS <- function(pts,m,tol=NULL) {
  # Check input
  ans <- (is.matrix(pts) | is.data.frame(pts))
  if (!ans) {
    msg <- "Input has to be either a matrix or a data frame\n"
    warning(msg)
    return(NULL)
  }
  
  # Number of data points (n)
  tmp <- dim(pts)
  n <- tmp[1]
  
  # Independent variable
  x <- pts[[1]]
  
  # Dependent variable
  y <- pts[[2]]
  
  # Prepare for multilinear regression
  M <- matrix(x^m,ncol=1)
  if (m > 1) {
    for (i in (m-1):1) M <- cbind(M,matrix(x^i,ncol=1))
  }
  M <- cbind(M,matrix(y,ncol=1))
  colnames(M) <- NULL
  rownames(M) <- NULL
  
  ### Multilinear regression (solveLS code) ###
  
  # Build A and y matrices
  ones <- matrix(rep(1,times=n),ncol=1)
  A <- as.matrix(M[,1:m])
  A <- cbind(A,ones)
  y <- matrix(M[,m+1],ncol=1)
  
  # Build F and g matrices to get solution through solve
  F <- t(A) %*% A
  g <- t(A) %*% y
  
  # Check whether F is singular
  d <- det(F)
  if (abs(d) < 1e-6) {
    msg <- paste0("There are infinite solutions to ",
                  "this least squares fitting.\n")
    message(msg)
    return(NULL)
  }
  
  # Change tolerance, if required
  if (is.null(tol)) tol <- 1e-17
  
  # Solution
  a <- solve(F,g,tol=tol)
  
  # Sum of squared residuals
  eps <- y-A %*% a
  d <- sum(eps^2)
  
  # Reshape for output
  a <- as.vector(a)
   
  return(list(a=a,SSE=d))
}


#' Find optimal polynomial model
#'
#' \code{which_poly} tries polynomial regression with polynomials
#' from degree 0 (a constant) to degree 6, on data provided. It
#' then outputs values of the variance of the residuals for each 
#' degree and displays a plot of the same versus the degree number, 
#' in an effort to suggest the degree of the best polynomial for the 
#' regression. The regression coefficients can then be calculated
#' with the function \code{\link{polysolveLS}}.
#' 
#' The ability of a polynomial regression to account for most data
#' variability, without including data noise is reflected in how
#' the variance,
#' \deqn{
#'   \sigma_e^2=(\sum_{i=1}^n \epsilon_i^2)/(n-m-1)
#' }
#' drops with the increasing degree of the polynomial used to 
#' perform the regression. A sudden drop, followed by values slowly
#' decreasing, or alternating slightly increasing and decreasing
#' behaviour, indicates that the degree corresponding to the sudden
#' drop belongs to the polynomial modelling most data variability
#' and neglecting data noise. As polynomial regression is normally
#' used with polynomials of degree up to 4 or 5, a default set of
#' polynomials up to degree 6 is here tried out. Degrees higher than
#' 6 can be forced by the user, but the risk with higher degrees is
#' that the system of normal equations connected with regression
#' becomes severly ill conditioned. In such situations the user 
#' should change the tolerance (\code{tol}) to values smaller than
#' the default \code{1e-17}. 
#' 
#' @param pts A \eqn{n \times 2} matrix or data frame where each row
#'            contains the coordinates of a data point used for
#'            regression.
#' @param mmax An integer. The highest degree of the polynomial to 
#'             be used to calculate the variance of the residuals.
#'             The default value is 6.
#' @param plt A logical variable to command the display of the
#'                plot of the variance vs the polynomials' degree.
#'                The default is \code{plt=TRUE}.
#' @param tol A real number. The solution of a linear system can be
#'            compromised when the condition number of the matrix of
#'            coefficients is particularly high (ill-conditioned
#'            matrices). \code{tol} is the reciprocal of the 
#'            condition number. For values of tol smaller than 1e-17,
#'            ill-conditioning is deemed to be sever enough not to
#'            guarantee an accurate solution. For such values the
#'            function stops execution, returning an error message.
#'            In fact, the solution can still be accurate,
#'            notwithstanding ill-conditioning, and the user can 
#'            force the calculation of a solution using a value of
#'            \code{tol} smaller than 1e-17. Default is \code{NULL},
#'            corresponding to a \code{tol=1e-17}.
#' @return A data frame with two columns, the first named \code{m}
#'         and including the degrees of all polynomials tested. The
#'         second called \code{sige} and including the value of the
#'         variances corresponding to all values of \code{m}. The
#'         function also displays a plot of \code{sige} vs \code{m},
#'         by default.
#'         
#' @examples 
#' # 21 points close to the quadratic x^2 - 5*x + 6
#' x <- seq(-2,5,length=21)
#' set.seed(7766)
#' eps <- rnorm(21,mean=0,sd=0.5)
#' y <- x^2-5*x+6+eps
#' 
#' # Data frame
#' pts <- data.frame(x=x,y=y)
#' 
#' # Try function without plot
#' ddd <- which_poly(pts,plt=FALSE)
#' print(ddd)
#' 
#' # Try function with plot and extending
#' # highest polynomials' degree to 10
#' ddd <- which_poly(pts,mmax=10)
#' 
#' @export
which_poly <- function(pts,mmax=6,plt=TRUE,tol=NULL) {
  # Vector of degrees
  dvec <- 0:mmax
  
  # Vector of variances
  vvec <- c()
  
  # Separate data points coordinates
  x <- pts[[1]]
  y <- pts[[2]]
  
  # Number of data points
  n <- length(x)
  
  # Polynomial with degree 0 is a special case
  eps <- y-mean(y)
  vvec <- c(vvec,sum(eps^2)/(n-1))
  
  # Loop to test remaining polynomials
  for (i in 1:mmax) {
    ltmp <- polysolveLS(pts,i)
    vvec <- c(vvec,ltmp$SSE/(n-i-1))
  }
  
  # Output
  ddd <- data.frame(m=dvec,sige=vvec)
  
  # Plot
  if (plt) {
    plot(ddd,type="b",ylab="Variance")
  }
  
  return(ddd)
}
  


###--------------------------------------------------------------
### Auxiliary functions (not for users
###--------------------------------------------------------------
