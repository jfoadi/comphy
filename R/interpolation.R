###
# This file is part of the package "comphy"
###


#' 1D linear interpolation
#'
#' Classic linear interpolation between two tabulated (known) points of a
#' one-variable function.
#' 
#' @param x0 A vector of real numbers. These are the grid points chosen for
#' the interpolation. All points of this grid need to be within the tabulated grid.
#' @param x A vector of real numbers. Grid points corresponding to the tabulated
#' (known) values of the function.
#' @param f A vector of real numbers. Tabulated (known) values of the function, 
#' corresponding to the grid \code{x}.
#' @return A vector of real numbers. These are the actual interpolated values
#' (calculated using linear interpolation), corresponding to all values of 
#' the grid \code{x0}.
#' @examples 
#' # Tabulated values: f(x) = 2*x^2-1
#' x <- c(0,1,3,7)
#' f <- 2*x^2-1
#' 
#' # Grid for interpolation
#' x0 <- seq(0,7,length=501)
#' 
#' # Interpolated points
#' f <- linpol(x,f,x0)
#' print(f)
#' 
#' @export 
linpol <- function(x,f,x0) {
  # Check sizes
  n <- length(x)
  if (length(f) != n) {
    stop("Incorrect size for f.")
  }
  m <- length(x0)
  
  # Check interpolating grid is within tabulated interval
  xmin <- min(x)
  xmax <- max(x)
  idx <- which(x0 < xmin | x0 > xmax)
  if (length(idx) > 0) {
    stop("Interpolating grid cannot overflow tabulated interval.")
  }
  
  # Sort x and f
  idx <- order(x)
  x <- x[idx]
  f <- f[idx]
  
  fint <- c()
  for (i in 1:m) {
    itmp <- which(x <= x0[i])
    idx1 <- itmp[length(itmp)]
    itmp <- which(x >= x0[i])
    idx2 <- itmp[1]
    
    # Tabulated interval points
    x1 <- x[idx1]
    x2 <- x[idx2]
    
    # Tabulated function points
    f1 <- f[idx1]
    f2 <- f[idx2]
    
    # First difference (x0-x2)
    dx1 <- x0[i]-x2
  
    # Second difference (x1-x0)
    dx2 <- x1-x0[i]
  
    # Bottom difference
    dx <- x1-x2
  
    # Interpolations
    if (abs(dx) < .Machine$double.eps) {
      fint <- c(fint,f1)
    } else {
      fint <- c(fint,(dx1*f1+dx2*f2)/dx)
    }
  }
  
  return(fint)
}


#' Neville-Aitken algorithm for polynomial interpolation
#'
#' Hierarchical series of linearly-interpolated \eqn{P_{ij}} values calculated using
#' Neville-Aitken's algorithm. In the \eqn{P_{ij}} expression, j is the level of the
#' algorithm and i the leftmost grid-point of the tabulated function points.
#' 
#' @param x A vector of real numbers. Grid points corresponding to the tabulated
#' (known) values of the function.
#' @param f A vector of real numbers. Tabulated (known) values of the function, 
#' corresponding to the grid \code{x}.
#' @param x0 A vector of real numbers. These are the grid points chosen for
#'           the interpolation. All points of this grid need to be 
#'           within the tabulated grid.
#' @return An upper triangular matrix of size n containing the
#'         linearly-interpolated values. \code{P[i,j]} is zero for 
#'         \eqn{i+j > n+1}.
#'            
#' @examples 
#' # Tabulated values: f(x) = x^3-2*x^2+3*x-1
#' x <- c(0.1,0.4,0.6,0.8,0.9)
#' f <- x^3-2*x^2+3*x-1
#' 
#' # Interpolation point
#' x0 <- 0.75
#' 
#' # Upper-triangular matrix of N-A values
#' P <- nevaitpol(x,f,x0)
#' 
#' # From level 4 onward the interpolated value
#' # does not change because f(x) is a 3rd-degree polynomial
#' print(P)
#' 
#' @export 
nevaitpol <- function(x,f,x0) {
  # Number of interpolating points
  n <- length(x)
  
  # Create interpolated-values matrix
  P <- matrix(data=rep(0,times=n*n),ncol=n,nrow=n)
  
  # Next three lines can be used in a future version of comphy
  ## Create Up and down corrections matrix
  #C <- matrix(data=rep(0,times=n*n),ncol=n,nrow=n)
  #D <- matrix(data=rep(0,times=n*n),ncol=n,nrow=n)
  
  # Change n into n-1 to align formula with the one in book
  n <- n-1
  
  # First column is filled with function values   
  P[,1] <- f
  
  # N-A algorithm
  for (j in 2:(n+1)) {
    for (i in 1:(n-j+2)) {
      P[i,j] <- ((x0-x[i+j-1])*P[i,j-1]+(x[i]-x0)*P[i+1,j-1])/(x[i]-x[i+j-1])
    }
  }
  
  # Next twelve lines can be used in a future version of comphy
  ## First set of corrections are the differences 
  ## with the tabulated point
  #C[,1] <- P[,1]-x
  #D[,1] <- P[,1]-x
  #
  ## All other corrections
  #for (j in 2:(n+1)) {
  #  for (i in 1:(n-j+2)) {
  #    C[i,j] <- P[i,j]-P[i,j-1]
  #    D[i,j] <- P[i,j]-P[i+1,j-1]
  #  }
  #}
  
  # Next line can be used in a future version of comphy
  #return(list(P=P,C=C,D=D))
  
  return(P)
}


#' Divided differences
#'
#' Calculation of all the n*(n+1)/2 divided differences related to n tabulated
#' points of a function. The values returned fill half of a n X n matrix, the other
#' half being filled with zeros.
#' 
#' @param x A vector of real numbers. Grid points corresponding to the tabulated
#' (known) values of the function.
#' @param f A vector of real numbers. Tabulated (known) values of the function, 
#' corresponding to the grid \code{x}.
#' 
#' @return A matrix of size n X n, where n is the length of \code{x}. Each column
#' of this matrix contains the divided differences at a specified level. Thus, column
#' 1 contains the level 1 values, i.e. the n tabulated points, column 2 contains the
#' n-1 divided differences calculated with the adjacent couples of grid points,
#' column 3 contains all n-2 level 3 divided differences, and so on. In each column
#' the remaining slots (no slots in the first column, one slot in the second column,
#' two slots in the third column, etc) are filled with zeros.
#' @examples 
#' # Tabulated values: f(x)=x^3-4x^2+3x+2
#' x <- c(-1,1,2,4)
#' f <- c(-6,2,0,14)
#' 
#' # Matrix filled with divided differences and zeros
#' P <- divdif(x,f)
#' print(P)
#' 
#' # Add two tabulated points to previous set
#' x <- c(x,0,3)
#' f <- c(f,2,2)
#' 
#' # New divided differences appear, but
#' # the old ones are unchanged
#' P <- divdif(x,f)
#' print(P)
#' 
#' @export
divdif <- function(x,f) {
  # Number of interpolating points
  n <- length(x)
  
  # Create interpolated-values matrix
  P <- matrix(data=rep(0,times=n*n),ncol=n,nrow=n)
  
  # Change n into n-1 to align formula with the one in book
  n <- n-1
  
  # First column is filled with function values   
  P[,1] <- f
  
  # N-A algorithm
  for (j in 2:(n+1)) {
    for (i in 1:(n-j+2)) {
      P[i,j] <- (P[i+1,j-1]-P[i,j-1])/(x[i+j-1]-x[i])
    }
  }
  
  return(P)
}


#' Approximating polynomial for divided differences
#'
#' Calculation of a polynomial of order n via divided differences. All n tabulated
#' points provided (with n greater or equal than 2) are used by default for the 
#' calculation, but the option is available to use only np points, where np must be 
#' greater or equal than 2. In case only part of the n available tabulated points is
#' used (np < n), the first two points are fixed to be equal to the smallest and largest
#' tabulated grid x points; the remaining np-2 points are selected randomly among the 
#' n-2 remaining ones. 
#' 
#' @param x A vector of real numbers. Grid points corresponding to the tabulated
#' (known) values of the function.
#' @param f A vector of real numbers. Tabulated (known) values of the function, 
#' corresponding to the grid \code{x}.
#' @param x0 A vector of real numbers. These are the grid points chosen for
#' the interpolation.
#' @param np An integer. The number of known points used for the
#'           interpolation. np > 2 because the smallest and largest
#'           value of \code{x} have to be always among the known
#'           points. Aside from the points at the extremes of the
#'           interpolation interval, the other points are chosen
#'           randomly.
#' 
#' @return A named list of length 3 and names \code{x}, \code{f} and \code{f0}. 
#' \describe{
#'    \item{\code{x}}{Tabulated grid points used for the interpolation.}
#'    \item{\code{f}}{Tabulated function points used for the interpolation. They 
#'                    correspond to \code{x}.}
#'    \item{\code{f0}}{Interpolated values. They correspond to the input vector
#'                     \code{x0}.}
#' }
#' @examples
#' # Tabulated grid points for function sin(x)
#' x <- seq(0,3*pi/2,length=20)
#' f <- sin(x)
#' 
#' # Grid of interpolated points
#' x0 <- seq(0,3*pi/2,length=200)
#' 
#' # Interpolation using all 20 tabulated points
#' ltmp <- polydivdif(x0,x,f)
#' plot(ltmp$x,ltmp$f,pch=16)
#' points(x0,ltmp$f0,type="l")
#' 
#' # Interpolation using only five points (dangerous!)
#' ltmp <- polydivdif(x0,x,f,np=5)
#' points(ltmp$x,ltmp$f,col=2,cex=1.5)
#' points(x0,ltmp$f0,type="l",col=2)
#' 
#' @export
polydivdif <- function(x0,x,f,np=length(x)) {
  # Size of interpolating grid, x0
  m <- length(x0)
  
  # Full Size of tabulated grid, x
  n <- length(x)
  
  # Select appropriate points for interpolation. Keep first
  # and last point and a random selection in the middle
  idx1 <- which(x == min(x))
  idx2 <- which(x == max(x))
  ridx <- 1:n
  if (np >= 3 & np < n) {
    ridx <- ridx[-c(idx1,idx2)]
    jidx <- sample(ridx,size=np-2,replace=FALSE)
    idx <- c(idx1,jidx,idx2)
  } else if (np == 2) {
    idx <- c(idx1,idx2)
  } else if (np == length(x)) {
    idx <- ridx
  } else {
    stop("np is an integer between 2 and the number of points x.")
  }
  x <- x[idx]
  f <- f[idx]
  
  # New, possibly reduced, size of tabulated grid, x 
  n <- length(x)
  
  # Divided differences
  P <- divdif(x,f)
  
  # Output 
  f0 <- c()
  
  # Contributions at all points of interpolating grid, x0
  for (ix in 1:m) {
    xp <- x0[ix]
    psum <- (xp-x[n-1])*P[1,n]+P[1,n-1]
    for (i in (n-1):2) {
      psum <- psum*(xp-x[i-1])+P[1,i-1]
    }
    f0 <- c(f0,psum)
  }
  
  return(list(x=x,f=f,f0=f0))
}


#' Degree of best-interpolating polynomial
#' 
#' The degree is chosen making use of divided differences. 
#' As more and more points are used for the interpolation, 
#' the components of columns of divided differences
#' corresponding to higher orders, are smaller and smaller. 
#' They are exactly zero when the function to interpolate is 
#' a polynomial of degree, say, \eqn{k}. More specifically, 
#' all divided differences of order \eqn{k+1} and above are 
#' exactly zero. The criterion used suggests a value \eqn{k} 
#' if the average of the absolute value of all divided 
#' differences of order \eqn{k+1} is less than a given 
#' threshold \code{thr} (default 1e-6).
#' 
#' The divided differences depend on the specific points 
#' selected to calculate the interpolated curve. To avoid 
#' potential bias that might occur when the tabulated 
#' points used are not distributed uniformly, several random 
#' selections of tabulated points are performed (default
#' \code{ntrial=30}) and the highest \eqn{k} is returned.
#' 
#' @param x A vector of real numbers. Grid points corresponding to the tabulated
#' (known) values of the function.
#' @param f A vector of real numbers. Tabulated (known) values of the function, 
#' corresponding to the grid \code{x}.
#' @param thr A real number. This is the threshold to decide
#'            when a column in the triangular matrix of
#'            divided differences has a small-enough average 
#'            value (small than \code{thr}). Default is
#'            \code{thr=1e-6}.
#' @param ntrial A positive integer to decide how many random
#'               selections of the provided known (tabulated)
#'               points have to be carried out. Default is
#'               \code{ntrial=30}.
#' @return An integer corresponding to the best interpolating
#'         polynomial's degree.
#' 
#' @examples
#' # Tabulated grid points for function cos(x)
#' x <- seq(0,3*pi/2,length=20)
#' f <- cos(x)
#' 
#' # Suggested polynomial degree (default ntrial)
#' k <- decidepoly_n(x,f)
#' print(k)
#' 
#' # Increase number of random selections (ntrial=50)
#' k <- decidepoly_n(x,f,ntrial=50)
#' print(k)
#' 
#' @export
decidepoly_n <- function(x,f,thr=1e-06,ntrial=30) {
  # Length of tabulated points
  n <- length(x)
  
  # Sort x from smallest to largest
  odx <- order(x)
  x <- x[odx]
  f <- f[odx]
  
  # If only two points, it's a straight line
  if (n < 2) {
    stop("At least 2 tabulated points needed!")
  }
  if (n == 2) {
    dpoly <- 1
    return(dpoly)
  }
  
  # Irregular grid for measuring errors (two times 
  # the tabulated points)
  x0 <- stats::runif(n=2*n,min=x[1],max=x[n])
  x0 <- sort(x0)
  
  # Try this incremental scheme a few times
  max_dpoly <- 0
  for (itrial in 1:ntrial) {
    # One point at a time is added (randomly) until 
    # error is smaller than thr
    cthr <- 10*thr  # Bigger value than given threshold
    idx <- c(1,n)
    while (cthr > thr) {
      idx <- c(idx,sample((1:n)[-idx],1))
      if (length(idx) == n) {
        dpoly <- n-1
        break
      }
      P <- divdif(x[idx],f[idx])
      dd <- 0
      for (i in 1:2*n) {
        DD <- abs(P[1,length(idx)])
        for (j in idx) {
          DD <- DD*(x0[i]-x[j])
        }
        DD <- abs(DD)
        if (DD > dd) dd <- DD
      }
      if (dd < cthr) cthr <- dd
      if (cthr <= thr) dpoly <- length(idx)-2
    }
    if (dpoly > max_dpoly) max_dpoly <- dpoly
  }
  dpoly <- max_dpoly
  
  return(dpoly)
}