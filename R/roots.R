###
# This file is part of the package "comphy"
###


#' Bisection method for roots
#'
#' Find the zero of a function of one variable, \eqn{f(x)}, given a starting
#' value close to the zero, or an interval including the zero.
#' 
#' Finding the zero of \eqn{f(x)} is equivalent to finding the roots of the
#' equation:
#' \deqn{
#'  f(x) = 0 
#' }
#' The algorithm used is based on the bisection method that needs an initial
#' interval within which the root is supposed to reside. When multiple roots
#' are involved, the algorithm will only find one among those inside the
#' chosen interval. The algorithm can be started also with just one value,
#' x0, supposedly close to the wanted root. In this case, an interval
#' is selected so that the function at the extremes of the interval has
#' opposite signs. If such an interval is not found, the function dumps a
#' warning message and returns NULL. The bisection method has a slow 
#' convergence rate and it does not converge at all in specific situations.
#' 
#' @param fn A function of one variable. If the function includes variables,
#'           these will have to be passed as additional variables, using the
#'           same names in the function's definition (see examples).
#' @param x0 A numeric variable. The starting value to find the initial
#'           interval, when a search interval is not provided. The default
#'           value is `x0=0`.
#' @param lB A numeric variable indicating the lower (left) extreme of 
#'           the search interval. If not given, this number will be selected
#'           starting from `x0` and in small steps `eps` of values
#'           smaller than `x0`, until a value of `lB` is found for
#'           which the function `f` has sign opposite of the sign it
#'           has at `rB`. Default is for `lB` not to be entered (`lB=NULL`).
#' @param rB Same as `lB`, but corresponding to the upper (right) extreme of
#'           the search interval. Default is for `rB` not to be entered 
#'           (`rB=NULL`).
#' @param tol A real number, in general a small number. The width of the
#'            smallest interval containing the zero of the function just
#'            before the algorithm stops. This means that the largest error
#'            \eqn{|x-x_t|} between the numerical value of the root found, 
#'            \eqn{x}, and its correct value, \eqn{x_t}, is tol. Default
#'            value is 1e-9.
#' @param imax A positive integer. The maximum number of bisections of the
#'             interval, while searching the zero of the function. The
#'             default value is 1e6, although convergence is normally 
#'             obtained with a number of bisections much smaller than `imax`.
#'             `imax` is important to stop search in those cases in which the
#'             function has no zeros in the search interval provided.
#' @param eps A real number. The step size needed for the selection of a
#'            search interval, when this is not provided. In such a 
#'            situation, symmetric intervals with increasing width around
#'            `x0` are considered where the left and right extremes are 
#'            `x0-i*eps` and `x0+i*eps`, respectively, where `i` is a positive
#'            integer, progressively increasing from 1 to the maximum
#'            allowed value `imax`. Search for the selected interval stops
#'            when is the signs of the function `f` calculated at the extremes 
#'            are opposite. If the search interval is not found, a warning
#'            message is printed and NULL is returned. Default value is
#'            0.1.
#' @param message A logical variable to state whether messages about the 
#'                root and the largest error have to be printed. The default
#'                is for the messages to be printed (message=TRUE).
#' @param logg A logical variable to state whether information on the series
#'             of bisected intervals is printed (TRUE) or not (FALSE). Default
#'             is for such information not to be printed (FALSE).
#' @param ... Parameters passed to function `fn`, if needed.
#' 
#' @return A numeric value, the zero of the function (or, equivalently, the
#'         root of the equation \eqn{f(x)=0}).
#'         
#' @examples 
#' # The quadratic equation x^2-5*x+6=0 has two roots, 2 and 3
#' ff <- function(x) return(x^2-5*x+6)
#' 
#' # Find root 2, starting from a single point close to 2
#' x0 <- 1
#' x <- roots_bisec(ff,x0=1)
#' print(x)
#' 
#' # Find root 3, using an interval (no message printing)
#' x <- roots_bisec(ff,lB=2.8,rB=4,message=FALSE)
#' print(x)
#' 
#' # Function with a parameter f(x) = exp(x) - k
#' ff <- function(x,k=2) return(exp(x)-k)
#' 
#' # Solution of exp(x)=3 is log(3)
#' x <- roots_bisec(ff,k=3)
#' print(log(3))
#' 
#' @export
roots_bisec <- function(fn,x0=0,lB=NULL,rB=NULL,
                        tol=1e-9,imax=1e6,eps=0.1,
                        message=TRUE,logg=FALSE,...) {
  # Check input
  ans <- is.function(fn)
  if (!ans) {
    msg <- paste0("Argument 'fn' must be a function.\n")
    warning(msg)
    
    return(NULL)
  }
  
  # If either of lB and rB is NULL, select two
  # values around x0 such that fn(lB) and fn(rB)
  # have different signs
  if (is.null(lB) | is.null(rB)) {
    sgno <- 1
    icyc <- 1
    lB <- x0
    rB <- x0
    while (sgno >= 0 & icyc < imax) {
      lB <- lB-eps
      rB <- rB+eps
      sgno <- sign(fn(lB,...))*sign(fn(rB,...))
      icyc <- icyc+1
    }
    
    # If failed to find interval, return with message
    if (icyc >= imax) {
      msg <- "Failed to find starting interval within given imax.\n"
      warning(msg)
      
      return(NULL)
    }
  }
  
  # Root search using bisection
  
  # This part is to swap lB and rB if rB < lB for some reason
  if (rB < lB) {
    A <- lB
    lB <- rB
    B <- A
  }
  
  # Can't give an interval of length 0
  if (lB == rB) {
    msg <- "Use different values of lB and rB.\n"
    warning(msg)
    
    return(NULL)
  }
  
  # Output starting search interval
  if (message) {
    msg <- sprintf("Searching interval: [%f,%f].\n",lB,rB)
    cat(msg)
  }
  
  # Create log (if required)
  if (logg) {
    llB <- lB
    rrB <- rB
  }
  
  # Actual search ...
  a <- lB
  b <- rB
  icyc <- 1
  while (abs(b-a) > tol & icyc < imax) {
    c <- (a+b)/2
    icyc <- icyc+1
    if (fn(a,...)*fn(c,...) == 0 | fn(b,...)*fn(c,...) == 0) {
      a <- c
      b <- c
    } else if (fn(a,...)*fn(c,...) < 0) {
      b <- c
    } else if (fn(b,...)*fn(c,...) < 0) {
      a <- c
    }
    
    # Store extremes (if required)
    if (logg) {
      llB <- c(llB,a)
      rrB <- c(rrB,b)
    }
  }
  err <- abs(lB-rB)/(2^icyc)
  
  # When a,b and c coincide
  if (abs(a-c) <= tol | abs(b-c) <= tol) err <- tol
  
  # Creates series-of-intervals dataframe (if required)
  if (logg) {
    derr <- abs(llB-rrB)
    cc <- 0.5*(llB+rrB)
    llg <- data.frame(Left=llB,Right=rrB,Root=cc,Difference=derr)
  }
  
  # Can't give an interval of length 0
  if (icyc >= imax) {
    msg <- "No root found between the given values of lB and rB: "
    msg <- paste0(msg,sprintf("[%f,%f]\n",lB,rB))
    warning(msg)
    
    # Dump intervals (if required)
    if (logg) print(llg)
    
    return(NULL)
  }
  
  if (message) {
    msg <- sprintf("The root is %f. The error is less than %e.\n",c,err)
    cat(msg)
  }
  
  # Dump intervals (if required)
  if (logg) print(llg)
  
  return(c)
}


#' Newton method for roots
#'
#' Find the zero of a function of one variable, \eqn{f(x)}, given a starting
#' value close to the zero, using Newton method.
#' 
#' Finding the zero of \eqn{f(x)} is equivalent to finding the roots of the
#' equation:
#' \deqn{
#'  f(x) = 0 
#' }
#' The algorithm used is based on Newton method that needs an initial guess,
#' x0, and the analytic expression of the function's first derivative. The 
#' method has a much faster convergence rate than both the bisection and secant
#' methods, but it does not converge when the initial guess or any other
#' subsequent approximations accidentally coincide with an optimal point of the
#' function, i.e. a point at which the first derivative is zero. The algorithm
#' can also potentially be stuck in an endless loop of repeating values for
#' special combinations of functions and initial guess.
#' 
#' @param f0 A function of one variable. If the function includes variables,
#'           these will have to be passed as additional variables, using the
#'           same names in the function's definition (see examples).
#' @param f1 A function equal to the first derivative of `f0`. Parameters that
#'           are potentially included in `f0`, must be also included in `f1`.
#' @param x0 A numeric variable. The initial guess starting Newton's algorithm.
#' @param tol A real small number. The smallest difference
#'            between the new zero's approximation and the previous one, above 
#'            which the algorithm keeps working. As soon as the difference is 
#'            less than `tol`, the algorithm stops and the current approximation
#'            is returned as the final approximation to the function's root.
#'            Default value is 1e-9.
#' @param ftol A real small number. When `ftol` is not NULL (default value),
#'             Newton's algorithm stops when \eqn{|f(x)| < ftol}. This
#'             parameter essentially introduces a different stopping criterion.
#' @param imax A positive integer. The maximum number of iterations of the
#'             algorithm. The default value is 1e6, although convergence is 
#'             normally obtained with a number of iterations much smaller than 
#'             imax. imax is important to stop search in those cases in which 
#'             the algorithm gets stuck in endless loops (non-convergence).
#' @param message A logical variable to state whether messages about the 
#'                root and the error have to be printed. The default
#'                is for the messages to be printed (`message=TRUE`).
#' @param logg A logical variable to state whether information on the series
#'             of approximating roots is printed (TRUE) or not (FALSE). Default
#'             is for such information not to be printed (FALSE).
#' @param ... Parameters passed to the two functions `f0` and `f1`, if any.
#' 
#' @return A numeric value, the zero of the function (or, equivalently, the
#'         root of the equation \eqn{f(x)=0}).
#'         
#' @examples 
#' # The quadratic equation x^2-5*x+6=0 has two roots, 2 and 3
#' f0 <- function(x) return(x^2-5*x+6)
#' 
#' # First derivative
#' f1 <- function(x) return(2*x-5)
#' 
#' # Find root 2, starting from a single point close to 2
#' x0 <- 1
#' x <- roots_newton(f0,f1,x0=1)
#' print(x)
#' 
#' # Find root 3 (no message printing)
#' x <- roots_newton(f0,f1,x0=4,message=FALSE)
#' print(x)
#' 
#' # Function with a parameter f(x) = exp(kx) - 2
#' f0 <- function(x,k=2) return(exp(k*x)-2)
#' 
#' # First derivative (it includes the parameter)
#' f1 <- function(x,k=2) return(k*exp(k*x))
#' 
#' # Solution of exp(2x)-2=0 is log(2)/2
#' x <- roots_newton(f0,f1,k=2)
#' print(log(2)/2)
#' 
#' @export
roots_newton <- function(f0,f1,x0=0,tol=1e-9,imax=1e6,ftol=NULL,
                        message=TRUE,logg=FALSE,...) {
  # Check input
  ans <- is.function(f0)
  if (!ans) {
    msg <- paste0("Argument 'f0' has to be a function.\n")
    warning(msg)
    
    return(NULL)
  }
  ans <- is.function(f1)
  if (!ans) {
    msg <- paste0("Argument 'f1' must be a function, the derivative of f0.\n")
    warning(msg)
    
    return(NULL)
  }
  
  # If x0 yields f0(x0)=0 stop (arbitrary number, very small)
  if (abs(f0(x0,...)) < 1e-15) {
    xr <- x0
    if (message) {
      msg <- sprintf("The root is %f. The error is less than 1e-15.\n",xr)
      cat(msg)
    }
    
    return(xr)
  }
  
  # If first derivative is close to zero (or zero!),
  # stop and ask to use a different initial guess
  if (abs(f1(x0,...)) <= 1e-15) {
    msg <- sprintf("The initial guess, x0, is not suitable to start\n")
    msg <- paste0(msg,sprintf("  the algorithm as the first derivative is zero.\n"))
    warning(msg)
    
    return(NULL)
  }
  
  # Create log of approximating roots, if rquired
  if (logg) {
    lxr <- x0
    ltol <- NA
  }
  
  # Stopping criterion depends on
  # whether ftol is used or not
  x1 <- x0 - f0(x0,...)/f1(x0,...)
  ncyc <- 1
  if (is.null(ftol)) {
    while (abs(x1-x0) >= tol & ncyc <= imax) {
      if (logg) {
        lxr <- c(lxr,x1)
        ltol <- c(ltol,abs(x0-x1))
      }
      x0 <- x1
      if (abs(f1(x0,...)) < 1e-15) {
        msg <- "The algorithm does not converge with this choice\n"
        msg <- paste0(msg,"of initial guess.\n")
        stop(msg)
      }
      x1 <- x0 - f0(x0,...)/f1(x0,...)
      ncyc <- ncyc+1
    }
  } else {
    while (abs(f0(x0,...)-f0(x1,...)) >= ftol & ncyc <= imax) {
      if (logg) {
        lxr <- c(lxr,x1)
        ltol <- c(ltol,abs(x0-x1))
      }
      x0 <- x1
      if (logg) lxr <- c(lxr,x0)
      if (abs(f1(x0,...)) < 1e-15) {
        msg <- "Algorithm diverges.\n"
        msg <- paste0(msg,sprintf("Current guess of root is %\f.\n",x0))
        warning(msg)
        
        return(x0)
      }
      x1 <- x0 - f0(x0,...)/f1(x0,...)
      ncyc <- ncyc+1
    }  
  }
  
  # Adjust roots in a data frame
  if (logg) {
    llg <- data.frame(Root=lxr,Shift=ltol)
  }
  
  # If ncyc > imax stop wasn't for convergence
  if (ncyc > imax) {
    msg <- "The algorithm has stopped because maximum number of iterations\n"
    msg <- paste0(msg,sprintf("  was reached. Current guess of root is %f.\n",x1))
    msg <- paste0(msg,"  In this case, he algorithm is likely not to converge.\n")
    warning(msg)
    if (logg) print(llg)
    xr <- x1
    
    return(xr)
  }
  
  # Algorithm has converged
  xr <- x0
  if (message) {
    if (is.null(ftol)) {
      msg <- sprintf("The root is %f. The error is less than %e.\n",xr,tol)
    }
    if (!is.null(ftol)) {
      msg <- sprintf("The root is %f. The error is less than %e.\n",xr,ftol)
    }
    cat(msg)
  }
  
  # Dump intervals (if required)
  if (logg) print(llg)
  
  return(xr)
}


#' Secant method for roots
#'
#' Find the zero of a function of one variable, \eqn{f(x)}, given a starting
#' value close to the zero, using the secant method.
#' 
#' Finding the zero of \eqn{f(x)} is equivalent to finding the roots of the
#' equation:
#' \deqn{
#'  f(x) = 0 
#' }
#' The algorithm used is essentially a reworking of Newton-Raphson, where the
#' first derivative is replaced by a finite difference computed with values
#' x0 and x1. Thus two values, x0 and x1, needs to be selected to start the
#' procedure. The convergence for this method is in general achieved faster 
#' than with the bisection method and slightly less fast than with Newton-Raphson. 
#' The algorithm can fail to converge when the secant in one of the iterations
#' is parallel to the x axis.
#' 
#' @param fn A function of one variable. If the function includes variables,
#'           these will have to be passed as additional variables, using the
#'           same names in the function's definition (see examples).
#' @param x0,x1 Two numeric variables. The initial guesses starting the algorithm.
#' @param ftol A real small number. The algorithm stops when 
#'             \eqn{|f(x)| < ftol}. Default value is 1e-09.
#' @param imax A positive integer. The maximum number of iterations of the
#'             algorithm. The default value is 1e6, although convergence is 
#'             normally obtained with a number of iterations much smaller than 
#'             imax. imax is important to stop search in those cases in which 
#'             the algorithm gets stuck in endless loops (non-convergence).
#' @param message A logical variable to state whether messages about the 
#'                root and the error have to be printed. The default
#'                is for the messages to be printed (message=TRUE).
#' @param logg A logical variable to state whether information on the series
#'             of approximating roots is printed (TRUE) or not (FALSE). Default
#'             is for such information not to be printed (FALSE).
#' @param ... Parameters passed to function `fn`, if needed.
#' 
#' @return A numeric value, the zero of the function (or, equivalently, the
#'         root of the equation \eqn{f(x)=0}).
#'         
#' @examples 
#' # The quadratic equation x^2-5*x+6=0 has two roots, 2 and 3
#' fn <- function(x) return(x^2-5*x+6)
#' 
#' # Find root 2, starting from two points at the left of 2
#' x0 <- 0
#' x1 <- 1
#' x <- roots_secant(fn,x0,x1)
#' print(x)
#' 
#' # Find root 3 (no message printing)
#' x0 <- 5
#' x1 <- 4
#' x <- roots_secant(fn,x0,x1,message=FALSE)
#' print(x)
#' 
#' # Function with a parameter f(x) = exp(kx) - 2
#' fn <- function(x,k=2) return(exp(k*x)-2)
#' 
#' # Solution of exp(2x)-2=0 is log(2)/2
#' x0 <- 0
#' x1 <- 1
#' x <- roots_secant(fn,x0,x1,k=2)
#' print(log(2)/2)
#' 
#' @export
roots_secant <- function(fn,x0,x1,imax=1e6,ftol=1e-09,
                         message=TRUE,logg=FALSE,...) {
  # Check input
  ans <- is.function(fn)
  if (!ans) {
    msg <- paste0("Argument 'fn' has to be a function.\n")
    warning(msg)
    
    return(NULL)
  }
  
  # Stop if the two initial values are too close
  if (abs(x0-x1) < 1e-15) {
    msg <- paste0("x0 and x1 have to be different.\n")
    warning(msg)
    
    return(NULL)
  }
  
  # Swap x0 with x1 if required
  if (abs(fn(x0,...)) < abs(fn(x1,...))) {
    tmp <- x0
    x0 <- x1
    x1 <- tmp
  }
  
  # Create log of approximating roots, if required
  if (logg) {
    lx0 <- x0
    lx1 <- x1 
  }
  
  # Generate xr, holding final root value
  xr <- x1
  
  # Main loop
  ncyc <- 1
  while(abs(fn(x1,...)) >= ftol & ncyc <= imax) {
    # Stop if secant is horizontal
    if (abs(fn(x0,...)-fn(x1,...)) < 1e-15) {
      # Adjust series of x0,x1 in a data frame
      if (logg) {
        llg <- data.frame(x0=lx0,x1=lx1)
        print(llg)
      }
      msg <- sprintf("The secant at cycle %d is parallel.\n",ncyc)
      msg <- paste0(msg,"Method has failed to converge.\n")
      warning(msg)
      
      return(NULL)
    }
    x2 <- x1-fn(x1,...)*(x0-x1)/(fn(x0,...)-fn(x1,...))
    x0 <- x1
    x1 <- x2
    xr <- x1
    ncyc <- ncyc+1
    if (logg) {
      lx0 <- c(lx0,x0)
      lx1 <- c(lx1,x1)
    }
  }
  
  # Adjust series of x0,x1 in a data frame
  if (logg) {
    llg <- data.frame(x0=lx0,x1=lx1)
  }
  
  # If algorithm has not converged return NULL
  if (ncyc >= imax) {
    if (logg) print(llg)
    msg <- "The method has failed to converge.\n"
    warning(msg)
    
    return(NULL)
  }
  
  # Algorithm has converged
  if (message) {
    msg <- sprintf("The root is %f. The error is less than %e.\n",xr,ftol)
    cat(msg)
  }
  
  # Dump values if required
  if (logg) {
    print(llg)
  }
  
  return(xr)
}

###--------------------------------------------------------------
### Auxiliary functions (not for users)
###--------------------------------------------------------------
