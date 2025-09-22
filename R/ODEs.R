###
# This file is part of the package "comphy"
###

#' Euler method for systems of ODEs
#'
#' Solves a system of \eqn{m} first-order ODEs using the explicit Euler method.
#' 
#' The method is accurate and stable when the stepsize \code{h} is
#' relatively small. The local error is \eqn{O(h^2)}, while the global
#' error is \eqn{O(h)}. Other numerical methods are generally used to 
#' calculate solutions with a higher accuracy.
#'
#' @param f A function of the form \code{f(t,y)} returning a numeric vector.
#'        It must be defined before using \code{EulerODE}. This function is
#'        the right hand side of the ODE, i.e. the gradient of the ODE system.
#' @param y0 A numeric vector with initial values (length = \eqn{m}).
#' @param h Step size.
#' @param t0 Initial time.
#' @param tf Final time.
#' @param ... Other parameters potentially needed by the gradient
#'        function.
#'
#' @return A list with elements \code{t} (time points) and \code{y} 
#'         (solution matrix). The first row of the matrix contains the
#'         initial values of \code{y} at time \code{t0}. Each column of the
#'         matrix contains the numerical solution for each one of the \eqn{m}
#'         functions of the system of ODEs.
#'         
#' @examples
#' 
#' # IVP: \eqn{dy/dt=6-2y,\ y(0)=0}.
#' # Define gradient
#' f <- function(t,y) {dy <- 6-2*y; return(dy)}
#' 
#' # Solution interval
#' t0 <- 0
#' tf <- 2
#' 
#' # Initial condition
#' y0 <- 0
#' 
#' # Step
#' h <- 0.1
#' 
#' # Numerical solution
#' ltmp <- EulerODE(f,t0,tf,y0,h)
#' 
#' # Print grid
#' print(ltmp$t)
#' 
#' # Print numerical solution
#' print(ltmp$y)
#' 
#' # Example with two ODEs. 
#' # \eqn{dy_1/dt=y_1+2y_2}
#' # \eqn{dy_2/dt=(3/2)y_1-y_2}
#' # \eqn{y_1(0)=1, y_2(0)=-2}
#' 
#' # Define gradient
#' dy <- function(t,y) {
#'   dy1 <- y[1]+2*y[2] 
#'   dy2 <- 1.5*y[1]-y[2] 
#'   return(c(dy1,dy2))
#' }
#' 
#' # Solution interval
#' t0 <- 0
#' tf <- 2
#' 
#' # Initial conditions
#' y0 <- c(1,-2)
#' 
#' # Step
#' h <- 0.1
#' 
#' # Numerical solution
#' ltmp <- EulerODE(dy,t0,tf,y0,h)
#' 
#' # Print grid
#' print(ltmp$t)
#' 
#' # Print numerical solution y1
#' print(ltmp$y[,1])
#' 
#' # Print numerical solution y2
#' print(ltmp$y[,2])
#'  
#' @export
EulerODE <- function(f,t0,tf,y0,h,...) {
  # Number of steps
  nsteps <- ceiling((tf-t0)/h)
  
  # Preallocate time and solution array
  t <- numeric(nsteps+1)
  y <- matrix(0,nrow=nsteps+1,ncol=length(y0))
  
  # Set initial conditions
  t[1] <- t0
  y[1,] <- y0
  
  # Euler iteration
  for (i in 1:nsteps) {
    t[i+1] <- t[i]+h
    y[i+1,] <- y[i,] + h*f(t[i],y[i,],...)
  }
  
  return(list(t=t,y=y))
}


#' Heun method for systems of ODEs
#'
#' Solves a system of \eqn{m} first-order ODEs using the Heun method
#' (also known as the improved Euler method).
#'
#' The method improves upon the Euler method by using an average of the 
#' slopes at the beginning and end of each time step. It is more accurate,
#' with local error \eqn{O(h^3)} and global error \eqn{O(h^2)}.
#'
#' @param f A function of the form \code{f(t,y)} returning a numeric vector.
#'        It must be defined before using \code{HeunODE}. This function is
#'        the right hand side of the ODE, i.e. the gradient of the ODE system.
#' @param y0 A numeric vector with initial values (length = \eqn{m}).
#' @param h Step size.
#' @param t0 Initial time.
#' @param tf Final time.
#' @param ... Other parameters potentially needed by the gradient
#'        function.
#'
#' @return A list with elements \code{t} (time points) and \code{y} 
#'         (solution matrix). The first row of the matrix contains the
#'         initial values of \code{y} at time \code{t0}. Each column of the
#'         matrix contains the numerical solution for each one of the \eqn{m}
#'         functions of the system of ODEs.
#'         
#' @examples
#' 
#' # IVP: \eqn{dy/dt=6-2y,\ y(0)=0}.
#' # Define gradient
#' f <- function(t,y) {dy <- 6-2*y; return(dy)}
#' 
#' # Solution interval
#' t0 <- 0
#' tf <- 2
#' 
#' # Initial condition
#' y0 <- 0
#' 
#' # Step
#' h <- 0.1
#' 
#' # Numerical solution
#' ltmp <- HeunODE(f,t0,tf,y0,h)
#' 
#' # Print grid
#' print(ltmp$t)
#' 
#' # Print numerical solution
#' print(ltmp$y)
#' 
#' # Example with two ODEs. 
#' # \eqn{dy_1/dt=y_1+2y_2}
#' # \eqn{dy_2/dt=(3/2)y_1-y_2}
#' # \eqn{y_1(0)=1, y_2(0)=-2}
#' 
#' # Define gradient
#' dy <- function(t,y) {
#'   dy1 <- y[1]+2*y[2] 
#'   dy2 <- 1.5*y[1]-y[2] 
#'   return(c(dy1,dy2))
#' }
#' 
#' # Solution interval
#' t0 <- 0
#' tf <- 2
#' 
#' # Initial conditions
#' y0 <- c(1,-2)
#' 
#' # Step
#' h <- 0.1
#' 
#' # Numerical solution
#' ltmp <- HeunODE(dy,t0,tf,y0,h)
#' 
#' # Print grid
#' print(ltmp$t)
#' 
#' # Print numerical solution y1
#' print(ltmp$y[,1])
#' 
#' # Print numerical solution y2
#' print(ltmp$y[,2])
#'
#' @export
HeunODE <- function(f,t0,tf,y0,h,...) {
  # Number of steps
  nsteps <- ceiling((tf-t0)/h)
  
  # Preallocate time and solution array
  t <- numeric(nsteps+1)
  y <- matrix(0,nrow=nsteps+1,ncol=length(y0))
  
  # Set initial conditions
  t[1] <- t0
  y[1,] <- y0
  
  # Heun iteration
  for (i in 1:nsteps) {
    t[i+1] <- t[i]+h
    yp <- y[i,]+h*f(t[i],y[i,],...)            # Predictor
    y[i+1,] <- y[i,]+h/2*(f(t[i],y[i,],...)+
                            f(t[i+1],yp,...))  # Corrector
  }
  
  return(list(t=t,y=y))
}


#' Runge-Kutta 4th order method for systems of ODEs
#'
#' Solves a system of \eqn{m} first-order ODEs using the classical 
#' fourth-order Runge-Kutta method.
#'
#' This method achieves high accuracy by evaluating the gradient 
#' function four times per step. It has local error \eqn{O(h^5)} 
#' and global error \eqn{O(h^4)}. It is one of the most widely used 
#' methods for solving initial value problems numerically.
#'
#' @param f A function of the form \code{f(t,y)} returning a numeric vector.
#'        It must be defined before using \code{RK4ODE}. This function is
#'        the right hand side of the ODE, i.e. the gradient of the ODE system.
#' @param y0 A numeric vector with initial values (length = \eqn{m}).
#' @param h Step size.
#' @param t0 Initial time.
#' @param tf Final time.
#' @param ... Other parameters potentially needed by the gradient
#'        function.
#'
#' @return A list with elements \code{t} (time points) and \code{y} 
#'         (solution matrix). The first row of the matrix contains the
#'         initial values of \code{y} at time \code{t0}. Each column of the
#'         matrix contains the numerical solution for each one of the \eqn{m}
#'         functions of the system of ODEs.
#'         
#' @examples
#' 
#' # IVP: \eqn{dy/dt=6-2y,\ y(0)=0}.
#' # Define gradient
#' f <- function(t,y) {dy <- 6-2*y; return(dy)}
#' 
#' # Solution interval
#' t0 <- 0
#' tf <- 2
#' 
#' # Initial condition
#' y0 <- 0
#' 
#' # Step
#' h <- 0.1
#' 
#' # Numerical solution
#' ltmp <- RK4ODE(f,t0,tf,y0,h)
#' 
#' # Print grid
#' print(ltmp$t)
#' 
#' # Print numerical solution
#' print(ltmp$y)
#' 
#' # Example with two ODEs. 
#' # \eqn{dy_1/dt=y_1+2y_2}
#' # \eqn{dy_2/dt=(3/2)y_1-y_2}
#' # \eqn{y_1(0)=1, y_2(0)=-2}
#' 
#' # Define gradient
#' dy <- function(t,y) {
#'   dy1 <- y[1]+2*y[2] 
#'   dy2 <- 1.5*y[1]-y[2] 
#'   return(c(dy1,dy2))
#' }
#' 
#' # Solution interval
#' t0 <- 0
#' tf <- 2
#' 
#' # Initial conditions
#' y0 <- c(1,-2)
#' 
#' # Step
#' h <- 0.1
#' 
#' # Numerical solution
#' ltmp <- RK4ODE(dy,t0,tf,y0,h)
#' 
#' # Print grid
#' print(ltmp$t)
#' 
#' # Print numerical solution y1
#' print(ltmp$y[,1])
#' 
#' # Print numerical solution y2
#' print(ltmp$y[,2])
#'
#' @export
RK4ODE <- function(f,t0,tf,y0,h,...) {
  # Number of steps
  nsteps <- ceiling((tf-t0)/h)
  
  # Preallocate time and solution array
  t <- numeric(nsteps+1)
  y <- matrix(0,nrow=nsteps+1,ncol=length(y0))
  
  # Set initial conditions
  t[1] <- t0
  y[1,] <- y0
  
  # RK4 iteration
  for (i in 1:nsteps) {
    t[i+1] <- t[i]+h
    
    k1 <- f(t[i], y[i,],...)
    k2 <- f(t[i] + h/2, y[i,] + h/2 * k1,...)
    k3 <- f(t[i] + h/2, y[i,] + h/2 * k2,...)
    k4 <- f(t[i] + h,   y[i,] + h   * k3,...)
    
    y[i+1,] <- y[i,] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
  }
  
  return(list(t=t, y=y))
}


#' Solves a second-order BVP using the shooting method
#'
#' Solves \eqn{y^{('')} = f(t,y,y^{'})} on the interval \eqn{[t_0,t_f]} with
#' boundary conditions \eqn{y(t_0) = y_0}, \eqn{y(t_f) = y_f}, using the
#' shooting method and a root finder. The associated IVP is solved using
#' 4th-order Runge-Kutta with \code{RK4ODE}. The second initial value is
#' found using the bisection method with \code{roots_bisec}.
#' 
#' It is important to consider the uniqueness of the solution of a BVP.
#' If the BVP admits infinitely many solutions (a family of solutions),
#' \code{BVPshoot2} will find only one of them, depending on what initial
#' condition for the first derivative of the associated IVP, was found
#' using the bisection method.
#'
#' @param f A function of the form \code{f(t,y,dy)} returning a numeric 
#'        scalar. This defines the second-order ODE.
#' @param t0 Initial time.
#' @param tf Final time.
#' @param y0 Boundary value at \eqn{t_0}, i.e. \eqn{y(t_0)=y_0}.
#' @param yf Boundary value at \eqn{t_f}, i.e. \eqn{y(t_f)=y_f}.
#' @param h Step size for the RK4 integration.
#' @param s_guess A numeric starting guess for \eqn{y'(t_0)} (default is 1),
#'        or a 2-element numeric vector giving a bracketing interval.
#' @param tol A numeric value that tells the bisection algorithms when to
#'        stop. Default is 1e-9.
#' @param ... Additional arguments passed to \code{f}.
#'
#' @return A list with elements \code{t} (time grid) and \code{y} (solution 
#'         matrix), where \code{y[,1]} contains \eqn{y(t)} and \code{y[,2]} 
#'         its derivative.
#'         
#' @examples
#'         
#' # y"+ y = 9*sin(2*t); y(0)=-1, y(3*pi/2)=0
#' # Unique solution: y(t) = 3*sin(2*t) - cos(t)
#'
#' # Define y"=f(t,y,y')
#' f <- function(t,y,dy) {ff <- -9*sin(2*t)-y; return(ff)}
#' 
#' # Solution interval
#' t0 <- 0
#' tf <- 3*pi/2
#' 
#' # Boundary values
#' y0 <- -1
#' yf <- 0
#' 
#' # Step size
#' h <- 0.01
#' 
#' # Solution
#' ltmp <- BVPshoot2(f,t0,tf,y0,yf,h)
#' 
#' # Check
#' # Number of steps
#' n <- length(ltmp$t)-1
#' print(c(ltmp$t[1],ltmp$t[n+1]))
#' print(c(ltmp$y[1,1],ltmp$y[n+1,1]))
#' 
#' @export
BVPshoot2 <- function(f,t0,tf,y0,yf,h,s_guess=1,tol=1e-9,...) {
  
  # System of first-order ODEs for RK4ODE
  system <- function(t,y) {
    dy1 <- y[2]
    dy2 <- f(t,y[1],y[2],...)
    return(c(dy1,dy2))
  }
  
  # Shooting function: difference from target boundary value at tf
  shootF <- function(s) {
    y_init <- c(y0,s)
    sol <- RK4ODE(system,t0,tf,y_init,h)
    return(tail(sol$y[,1],1)-yf)
  }
  
  # Determine root-finding mode based on length of s_guess
  if (length(s_guess) == 2) {
    s_star <- roots_bisec(fn=shootF,lB=s_guess[1],rB=s_guess[2],
                          message=FALSE,tol=tol)
  } else if (length(s_guess) == 1) {
    s_star <- roots_bisec(fn=shootF,x0=s_guess,
                          message=FALSE,tol=tol)
  } else {
    stop("s_guess must be either a single number or a numeric vector of length 2.")
  }
  
  # Check for failure
  if (is.null(s_star)) {
    warning("Shooting method failed: root-finding did not converge.")
    return(NULL)
  }
  
  # Final solution with best initial slope
  y_init <- c(y0,s_star)
  out <- RK4ODE(system,t0,tf,y_init,h)
  
  return(out)
}


#' Linear shooting method for second-order linear BVPs
#'
#' Solves a second-order linear boundary value problem using the linear
#' shooting method and superposition of two initial value problems.
#' 
#' If the solution of the associated homogeneous IVP is very small (close
#' to zero) at the second boundary (tf), the solution becomes unstable and
#' the function stops with a warning. Other methods must be used in those
#' cases.
#'
#' @param f A function of the form \eqn{f(t,y,y')} representing 
#'        the second-order ODE: \eqn{y''=f(t,y,y')}.
#' @param t0 Initial time.
#' @param tf Final time.
#' @param y0 Boundary value at \code{t0}, i.e. \eqn{y(t_0) = y_0}.
#' @param yf Boundary value at \code{tf}, i.e. \eqn{y(t_f) = y_f}.
#' @param h Step size.
#' @param ... Optional parameters passed to the gradient function 
#'        \code{f}.
#'
#' @return A list with elements \code{t} (time points) and \code{y} 
#'         (solution matrix). The first column of matrix \code{y} is
#'         the solution, \eqn{y(t)}, the second is its first derivative,
#'         \eqn{y'(t)}.
#'
#' @examples
#'
#' # Solve: y'' - (3/x)y' + (4/x^2)y = x
#' # with y(1) = 0, y(2) = 4*(log(2) + 1)
#' # Exact solution: y(x) = x^2*(log(x) - 1) + x^3
#'
#' # Gradient
#' f <- function(x,y,dy,...) {
#'   (3/x)*dy-(4/x^2)*y+x
#' }
#'
#' t0 <- 1
#' tf <- 2
#' y0 <- 0
#' yf <- 4*(log(2)+1)
#' h <- 0.01
#'
#' ltmp <- BVPlinshoot2(f,t0,tf,y0,yf,h)
#'
#' # Checks
#' n <- length(ltmp$t)-1
#' print(c(ltmp$t[1],ltmp$t[n+1]))
#' print(c(ltmp$y[1,1],ltmp$y[n+1,1]))
#'
#' @export
BVPlinshoot2 <- function(f,t0,tf,y0,yf,h,...) {
  # Define system of first-order ODEs
  system_u <- function(t,y) {
    dy1 <- y[2]
    dy2 <- f(t,y[1],y[2],...)
    
    return(c(dy1,dy2))
  }
  
  # Define system of ODEs for homogeneous part (v)
  fhom <- function(t,y,dy,...) {
    f(t,y,dy,...)-f(t,0,0,...)
  }
  system_v <- function(t,y) {
    dy1 <- y[2]
    dy2 <- fhom(t,y[1],y[2],...)
    
    return(c(dy1,dy2))
  }
  
  # Solve u'' = f, with u(t0) = y0, u'(t0) = 0
  y_init_u <- c(y0,0)
  sol_u <- RK4ODE(system_u,t0,tf,y_init_u,h)
  
  # Solve v'' = f_hom, with v(t0) = 0, v'(t0) = 1
  y_init_v <- c(0,1)
  sol_v <- RK4ODE(system_v,t0,tf,y_init_v,h)
  
  # Combine solution: y(t) = u(t) + s * v(t)
  u <- sol_u$y[,1]
  v <- sol_v$y[,1]
  
  # Combine solution: y'(t) = u'(t) + s * v'(t)
  du <- sol_u$y[,2]
  dv <- sol_v$y[,2]
  
  # Check v(b) neq 0
  utf <- tail(u,1)
  vtf <- tail(v,1)
  if (abs(vtf) < 1e-12) 
    stop("Homogeneous solution vanishes at tf: cannot divide.")
  
  # Determine appropriate s and thus solutions (combined into a matrix)
  s <- (yf-utf)/vtf
  y_combined <- u+s*v
  dy_combined <- du+s*dv
  
  return(list(t=sol_u$t,y=cbind(y_combined,dy_combined)))
}


#' Sturm–Liouville eigenproblem with homogeneous Dirichlet boundary conditions
#'
#' Solves
#' \deqn{-\frac{d}{dx}\left(p(x)\,y'(x)\right) + q(x)\,y(x) = \lambda\, w(x)\,y(x)}
#' on \eqn{[a,b]} with \eqn{y(a)=0} and \eqn{y(b)=0}. The equation is discretised
#' on the interior nodes of a \strong{uniform} grid and assembled into
#' matrices \code{K} and \code{W} so that \code{K u = lambda W u}.
#' The problem is reduced to a symmetric standard eigenproblem and solved.
#'
#' Coefficients may be given as functions or numeric vectors:
#' \itemize{
#'   \item \code{p}: function on midpoints or numeric vector of length \code{length(x)-1} (midpoints).
#'   \item \code{q}, \code{w}: functions on nodes or numeric vectors of length \code{length(x)} (nodes).
#' }
#'
#' Homogeneous Dirichlet conditions are enforced by construction: unknowns
#' are interior only; the returned full eigenfunctions have zero endpoints.
#'
#' @param p Function \code{p(x)} or numeric vector at midpoints.
#' @param q Function \code{q(x)} or numeric vector at nodes.
#' @param w Function \code{w(x)} or numeric vector at nodes.
#' @param x Numeric grid including endpoints (\code{x[1]=a}, \code{x[n+1]=b}); must be uniform.
#' @param nev Integer number of eigenpairs to return (smallest); default all interior modes.
#' @param normalize Logical; if \code{TRUE}, scale interior eigenvectors so that
#'   \eqn{\sum_i h\,w_i\,u_i^2 = 1}. Default \code{TRUE}.
#' @param return_matrices Logical; if \code{TRUE}, also return \code{K} and \code{W}. Default \code{FALSE}.
#' @param check_inputs Logical; run basic checks (uniform grid, positivity of \code{p}, \code{w}). Default \code{TRUE}.
#' @param tol_uniform Tolerance for uniform‑grid check. Default \code{1e-12}.
#' @return A list with
#' \itemize{
#'   \item \code{values}: eigenvalues (ascending).
#'   \item \code{vectors_interior}: interior eigenvectors (matrix \code{(n-1) x k}).
#'   \item \code{vectors_full}: full eigenfunctions with zero endpoints (matrix \code{(n+1) x k}).
#'   \item \code{x}, \code{h}, \code{nev_used}.
#'   \item \code{K}, \code{W} if \code{return_matrices=TRUE}.
#' }
#' @examples
#' # p=1, q=0, w=1 on [0, pi]  -> eigenvalues ~ 1^2, 2^2, 3^2, ...
#' a <- 0; b <- pi; n <- 200
#' x <- seq(a, b, length.out = n+1)
#' pfun <- function(s) 1            # scalars are accepted; will be replicated
#' qfun <- function(s) 0
#' wfun <- function(s) 1
#' ep <- EPSturmLiouville2(pfun, qfun, wfun, x, nev = 4, normalize = TRUE)
#' round(ep$values, 3)              # ~ c(1, 4, 9, 16)
#'
#' @export
EPSturmLiouville2 <- function(
    p, q, w, x,
    nev = NULL,
    normalize = TRUE,
    return_matrices = FALSE,
    check_inputs = TRUE,
    tol_uniform = 1e-12
) {
  # ----- grid bookkeeping -----
  x <- as.numeric(x)
  if (length(x) < 3L) 
    stop("x must have at least 3 points (two boundaries + one interior).")
  dx <- diff(x); h <- dx[1L]
  if (max(abs(dx-h)) > tol_uniform)
    stop("EPSturmLiouville2 requires a uniform grid x.")
  n_plus_1 <- length(x); n <- n_plus_1-1L; m <- n-1L
  if (m < 1L) stop("Not enough interior points.")
  ii <- 2L:n
  xi <- x[ii]
  xmid <- (x[-1L] + x[-n_plus_1])/2  # length n
  
  # ----- helpers to evaluate coefficients (functions or vectors) -----
  # If a function returns length 1, replicate to requested length.
  eval_on <- function(f_or_v,where,name) {
    if (is.function(f_or_v)) {
      val <- f_or_v(where)
      if (length(val) == 1L) val <- rep(val,length(where))
      if (length(val) != length(where))
        stop(sprintf("%s(x) must return length %d.", name, length(where)))
      as.numeric(val)
    } else if (is.numeric(f_or_v)) {
      as.numeric(f_or_v)
    } else stop(sprintf("%s must be a function or numeric vector.",name))
  }
  
  # p at midpoints (length n)
  if (is.function(p)) {
    p_half_all <- eval_on(p,xmid,"p")
  } else {
    if (length(p) != length(xmid))
      stop("If numeric, p must have length length(x)-1 (midpoint values).")
    p_half_all <- as.numeric(p)
  }
  
  # q,w at nodes (length n+1); take interior slice below
  if (is.function(q)) {
    q_all <- eval_on(q,x,"q")
  } else {
    if (length(q) != length(x)) 
      stop("If numeric, q must have length length(x) (node values).")
    q_all <- as.numeric(q)
  }
  if (is.function(w)) {
    w_all <- eval_on(w,x,"w")
  } else {
    if (length(w) != length(x)) 
      stop("If numeric, w must have length length(x) (node values).")
    w_all <- as.numeric(w)
  }
  
  # interior-aligned arrays
  qi <- q_all[ii]                  # length m
  wi <- w_all[ii]                  # length m
  p_plus  <- p_half_all[ii]        # p_{i+1/2}, i=2..n
  p_minus <- p_half_all[ii-1L]     # p_{i-1/2}, i=2..n
  
  # ----- basic checks -----
  if (check_inputs) {
    if (any(!is.finite(qi))) stop("Non-finite q(x) on interior grid.")
    if (any(!is.finite(wi))) stop("Non-finite w(x) on interior grid.")
    if (any(!is.finite(p_plus)) || any(!is.finite(p_minus)))
      stop("Non-finite p(x) at midpoints.")
    if (any(p_plus <= 0) || any(p_minus <= 0))
      stop("p(x) must be positive on (a,b).")
    if (any(wi <= 0))
      stop("w(x) must be positive on (a,b).")
  }
  
  # ----- assemble K and W on interior -----
  K <- matrix(0.0,m,m)
  diag(K) <- (p_plus+p_minus)/(h*h)+qi
  if (m > 1L) {
    K[cbind(1:(m-1L),2:m)] <- -p_plus[-m] /(h*h)   # upper
    K[cbind(2:m,1:(m-1L))] <- -p_minus[-1L]/(h*h)  # lower
  }
  W <- diag(wi,nrow=m,ncol=m)
  
  # ----- symmetric reduction and eigen solve -----
  inv_sqrt_w <- 1/sqrt(wi)
  A <- K
  A <- inv_sqrt_w*A
  A <- t(inv_sqrt_w*t(A))
  ee <- eigen(A,symmetric=TRUE)
  vals <- Re(ee$values); vecs <- Re(ee$vectors)
  ord <- order(vals); vals <- vals[ord]; vecs <- vecs[,ord,drop=FALSE]
  
  # take first nev if requested
  if (!is.null(nev)) {
    if (nev < 1L) stop("nev must be >= 1.")
    nev <- min(nev,length(vals))
    vals <- vals[seq_len(nev)]
    vecs <- vecs[,seq_len(nev),drop = FALSE]
  } else {
    nev <- length(vals)
  }
  
  # recover interior eigenvectors u = W^{-1/2} z
  Ui <- inv_sqrt_w*vecs
  
  # normalise in discrete L2(w) if requested
  if (isTRUE(normalize)) {
    norms2 <- colSums(matrix(wi,nrow=m,ncol=nev)*(Ui*Ui))*h
    Ui <- sweep(Ui,2L,sqrt(norms2),"/")
  }
  
  # pad boundary zeros and assert homogeneous BCs
  Ufull <- matrix(0.0,n_plus_1,nev)
  Ufull[ii,] <- Ui
  if (check_inputs) {
    if (!all(Ufull[1,] == 0) || !all(Ufull[n_plus_1,] == 0))
      stop("Homogeneous Dirichlet endpoints must be zero (internal check failed).")
  }
  
  out <- list(
    values = vals,
    vectors_interior = Ui,
    vectors_full = Ufull,
    x = x,
    h = h,
    nev_used = nev
  )
  if (isTRUE(return_matrices)) { out$K <- K; out$W <- W }
  class(out) <- c("comphy_sl_ep", class(out))
  out
}





###--------------------------------------------------------------
### Auxiliary functions (not for users)
###--------------------------------------------------------------

