###
# This file is part of the package "comphy"
###


#' Gaussian Elimination
#'
#' Solution of a system of \eqn{n} equations in \eqn{n} unknowns,
#' using Gaussian elimination.
#' 
#' The linear system to solve is \eqn{Ax=b}, where \eqn{A} is
#' the \eqn{n\times n} matrix of coefficients of the \eqn{n}
#' unknowns in the \eqn{n\times 1} vector \eqn{x}, and \eqn{b}
#' is the \eqn{n\times 1} vector of known numbers. Gaussian
#' elimination consists of a series of so-called row operations
#' that transform \eqn{A} in an upper-triangular matrix. The
#' system corresponding to the transformed matrix can be solved
#' very quickly.
#' 
#' @param M The \eqn{n\times (n+1)} augmented matrix of
#'        coefficients corresponding to the system of \eqn{n}
#'        linear equations in \eqn{n} unknowns, \eqn{Ax=b}. 
#' @return A vector of length \eqn{n} containing the \eqn{n}
#'         numeric solutions for the \eqn{n} unknowns. If the
#'         system has no solutions or an infinite number of
#'         solutions, the function returns NULL and dumps a
#'         warning message.
#'         
#' @examples 
#' # System of three equations in three unknowns
#' #
#' # 3x_1 + x_2 +  x_3 = 6
#' #  x_1 - x_2 + 2x_3 = 4
#' # -x_1 + x_2 +  x_3 = 2
#' 
#' # Augmented matrix M=(A|b)
#' M <- matrix(c(3,1,-1,1,-1,1,1,2,1,6,4,2),ncol=4)
#' 
#' # Solution via Gauss elimination
#' x <- gauss_elim(M)
#' print(x)
#' 
#' @export
gauss_elim <- function(M) {
# Number of rows
m <- nrow(M)

# Number of columns
n <- ncol(M)

# One column more than # of rows
if (n != (m+1)) {
  warning("Input is not an augmented matrix.\n")
  return(NULL)
}

# Reduction to upper-triangular
M <- transform_upper(M)

# If any of the elements on the diagonal is zero,
# the system has no solutions or infinite solutions
flag <- FALSE
for (i in 1:m) {
  if (abs(M[i,i]) < 1e-6) {
    flag <- TRUE
  }
  if (flag) {
    message("This system has no solution or infinite solutions.\n")
    return(NULL)
  }
}

# Solution
x <- rep(0,length=m)
x[m] <- M[m,n]/M[m,m]
for (i in (m-1):1) {
  x[i] <- (M[i,n]-sum(x[(i+1):m]*M[i,(i+1):m]))/M[i,i]
}

return(x)
}


#' Transform to upper triangular
#' 
#' Transform an \eqn{n\times n} matrix to upper-triangular
#' form, using a series of row operations.
#' 
#' The algorithm used for the transformation is Gauss 
#' elimination, which makes use of row operations. If the input
#' matrix has \eqn{n+1} columns, the transformed 
#' \eqn{n\times (n+1)} matrix can be used to find the solution
#' of the associated system of linear equations.
#' 
#' @param M An \eqn{n\times n} or \eqn{n\times (n+1)} matrix.
#' @return The transformed \eqn{n\times n} or 
#'         \eqn{n\times (n+1)} matrix.
#'         
#' @examples 
#' # 3X3 matrix
#' #
#' # [ 3  1  1
#' #   1 -1  2
#' #  -1  1  1]
#' 
#' # Input matrix
#' A <- matrix(c(3,1,2,1,1,-1,1,1,-1),ncol=3)
#' 
#' # Upper-triangular matrix
#' U <- transform_upper(A)
#' print(U)
#' 
#' @export
transform_upper <- function(M) {
  # Number of rows
  n <- nrow(M)

  # Main loop
  for (i in 1:(n-1)) {
    # Swap rows to have largest values on diagonal
    idx <- which(abs(M[i:n,i]) == max(abs(M[i:n,i])))
    idx <- idx[length(idx)] # If more than one, pick the last
    idx <- i+idx-1
    N <- M[i,]
    M[i,] <- M[idx,]
    M[idx,] <- N
    
    # Next, perform row operations for column i
    for (j in (i+1):n) {
      if (abs(M[j,i]) > 1e-6) {
        cff <- -M[j,i]/M[i,i]
        M[j,] <- cff*M[i,]+M[j,]
      }
    }
  }
  
  return(M)
}


#' LU decomposition
#' 
#' Transform an \eqn{n\times n} matrix into a product of a 
#' lower-triangular and upper-triangular matrices, using the
#' Crout (method="crout" - default) or Doolittle (method=
#' "doolittle") method.
#' 
#' The "crout" method returns the upper triangular matrix, U, 
#' with ones on its diagonal. The "doolittle" method returns
#' the lower triangular matrix, L, with ones on its diagonal.
#' 
#' Some matrices do not have an LU decomposition unless a row
#' permutation is done to the matrix. In this function, the 
#' order of such a permutation is included in the named vector
#' \code{ord}, returned as part of the output. When the vector
#' is equal to 1,2,...,n (first n numbers, naturally ordained),
#' this means that there was no need of permuting the original
#' matrix to carry out the LU decomposition.
#' 
#' @param A An \eqn{n\times n} matrix.
#' @param method A character string. This calls two different
#'        procedures for the decomposition. Only "crout" 
#'        (default) and "doolittle" are recognised methods. A
#'        different character string forces the function to
#'        return NULL.
#' @return A named list with the lower triangular, L,
#'         upper triangular, U, matrices, and with a vector,
#'         ord, containing the permutation needed to achieve
#'         the LU factorisation.
#'         
#' @examples 
#' # 3X3 matrix
#' #
#' # [ 3  1  1
#' #   1 -1  2
#' #  -1  1  1]
#' 
#' # Input matrix
#' A <- matrix(c(3,1,-1,1,-1,1,1,2,1),ncol=3)
#' 
#' # LU decomposition
#' ltmp <- LUdeco(A)
#' print(ltmp$L)
#' print(ltmp$U)
#' print(ltmp$ord) # No permutation needed
#' 
#' # The product is the original matrix, A
#' print(ltmp$L%*%ltmp$U)
#' 
#' # Singular matrix with LU decomposition
#' A <- matrix(c(1,0,0,0,1,1,1,0,0),ncol=3)
#' print(det(A))
#' ltmp <- LUdeco(A,"doolittle")
#' print(ltmp$L)
#' print(ltmp$U)
#' print(ltmp$ord) # No permutation needed
#' 
#' # The product is the original matrix, A
#' print(ltmp$L%*%ltmp$U)
#' 
#' # Singular matrix without LU decomposition
#' A <- matrix(c(1,0,0,0,0,0,0,0,0),ncol=3)
#' ltmp <- LUdeco(A)
#' print(ltmp)
#' #
#' @export
LUdeco <- function(A,method="crout") {
  # Quit if wrong method
  if (method != "crout" & method != "doolittle") {
    msg <- "Not a recognised decomposition method.\n"
    warning(msg)
    
    return(NULL)
  }
  
  # Number of rows/columns
  n <- nrow(A)
  
  # Initialise order of rows
  rowsord <- 1:n
  
  ######### Crout #########
  if (method == "crout") {
    # Initial matrices (and first n + n -1 components)
    L <- matrix(rep(0,times=n*n),ncol=n)
    L[,1] <- A[,1]
    U <- diag(n)
    U[1,2:n] <- A[1,2:n]/A[1,1]
    
    # Remaining components
    for (j in 2:n) {
      
      # First set of relations for L
      for (i in j:n) {
        sss <- sum(L[i,1:(j-1)]*U[1:(j-1),j])
        L[i,j] <- A[i,j]-sss
      }
      
      # Second set of relations for U
      # Here i refers to the columns of U
      if (j < n) {
        for (i in (j+1):n) {
          sss <- sum(L[j,1:(i-1)]*U[1:(i-1),i])
          U[j,i] <- (A[j,i]-sss)/L[j,j]
        }
      }
    }
    
    # Check a permutation is needed
    nPerm <- FALSE
    nInf <- sum(is.infinite(L))
    if (nInf > 0) nPerm <- TRUE
    nInf <- sum(is.infinite(U))
    if (nInf > 0) nPerm <- TRUE
    nNan <- sum(is.nan(L))
    if (nNan > 0) nPerm <- TRUE
    nNan <- sum(is.nan(U))
    if (nNan > 0) nPerm <- TRUE
    if (nPerm) {
      for (i in 1:(n-1)) {
        # Swap rows to have largest values on diagonal
        idx <- which(abs(A[i:n,i]) == max(abs(A[i:n,i])))
        idx <- idx[length(idx)] # If more than one, pick the last
        idx <- i+idx-1
        N <- A[i,]
        A[i,] <- A[idx,]
        A[idx,] <- N
        idxN <- rowsord[i]
        rowsord[i] <- idx
        rowsord[idx] <- idxN
      }
      
      # Re-do LU decomposition on new A
      # Initial matrices (and first n + n -1 components)
      L <- matrix(rep(0,times=n*n),ncol=n)
      L[,1] <- A[,1]
      U <- diag(n)
      U[1,2:n] <- A[1,2:n]/A[1,1]
      
      # Remaining components
      for (j in 2:n) {
        
        # First set of relations for L
        for (i in j:n) {
          sss <- sum(L[i,1:(j-1)]*U[1:(j-1),j])
          L[i,j] <- A[i,j]-sss
        }
        
        # Second set of relations for U
        # Here i refers to the columns of U
        if (j < n) {
          for (i in (j+1):n) {
            sss <- sum(L[j,1:(i-1)]*U[1:(i-1),i])
            U[j,i] <- (A[j,i]-sss)/L[j,j]
          }
        }
      }
    }
    
    # Check again if NaN or Inf are present
    nPerm <- FALSE
    nInf <- sum(is.infinite(L))
    if (nInf > 0) nPerm <- TRUE
    nInf <- sum(is.infinite(U))
    if (nInf > 0) nPerm <- TRUE
    nNan <- sum(is.nan(L))
    if (nNan > 0) nPerm <- TRUE
    nNan <- sum(is.nan(U))
    if (nNan > 0) nPerm <- TRUE
    if (nPerm) {
      msg <- paste("The input matrix is singular",
                   "and it does not have\n",
                   "an LU decomposition.\n")
      message(msg)
      return(NULL)
    }
  }
  
  ######### Doolittle #########
  if (method == "doolittle") {
    # Initial matrices (and first n + n -1 components)
    U <- matrix(rep(0,times=n*n),ncol=n)
    U[1,] <- A[1,]
    L <- diag(n)
    L[2:n,1] <- A[2:n,1]/A[1,1]
    
    # Remaining components
    for (j in 2:n) {
      
      # First set of relations for U
      # Here i refers to the columns of U
      for (i in j:n) {
        sss <- sum(L[j,1:(j-1)]*U[1:(j-1),i])
        U[j,i] <- A[j,i]-sss
      }
      
      # Second set of relations for L
      if (j < n) {
        for (i in (j+1):n) {
          sss <- sum(L[i,1:(i-1)]*U[1:(i-1),j])
          L[i,j] <- (A[i,j]-sss)/U[j,j]
        }
      }
    }
    
    # Check a permutation is needed
    nPerm <- FALSE
    nInf <- sum(is.infinite(L))
    if (nInf > 0) nPerm <- TRUE
    nInf <- sum(is.infinite(U))
    if (nInf > 0) nPerm <- TRUE
    nNan <- sum(is.nan(L))
    if (nNan > 0) nPerm <- TRUE
    nNan <- sum(is.nan(U))
    if (nNan > 0) nPerm <- TRUE
    if (nPerm) {
      for (i in 1:(n-1)) {
        # Swap rows to have largest values on diagonal
        idx <- which(abs(A[i:n,i]) == max(abs(A[i:n,i])))
        idx <- idx[length(idx)] # If more than one, pick the last
        idx <- i+idx-1
        N <- A[i,]
        A[i,] <- A[idx,]
        A[idx,] <- N
        idxN <- rowsord[i]
        rowsord[i] <- idx
        rowsord[idx] <- idxN
      }
      
      # Re-do LU decomposition on new A
      # Initial matrices (and first n + n -1 components)
      U <- matrix(rep(0,times=n*n),ncol=n)
      U[1,] <- A[1,]
      L <- diag(n)
      L[2:n,1] <- A[2:n,1]/A[1,1]
      
      # Remaining components
      for (j in 2:n) {
        
        # First set of relations for U
        # Here i refers to the columns of U
        for (i in j:n) {
          sss <- sum(L[j,1:(j-1)]*U[1:(j-1),i])
          U[j,i] <- A[j,i]-sss
        }
        
        # Second set of relations for L
        if (j < n) {
          for (i in (j+1):n) {
            sss <- sum(L[i,1:(i-1)]*U[1:(i-1),j])
            L[i,j] <- (A[i,j]-sss)/U[j,j]
          }
        }
      }
    }
    
    # Check again if NaN or Inf are present
    nPerm <- FALSE
    nInf <- sum(is.infinite(L))
    if (nInf > 0) nPerm <- TRUE
    nInf <- sum(is.infinite(U))
    if (nInf > 0) nPerm <- TRUE
    nNan <- sum(is.nan(L))
    if (nNan > 0) nPerm <- TRUE
    nNan <- sum(is.nan(U))
    if (nNan > 0) nPerm <- TRUE
    if (nPerm) {
      msg <- paste("The input matrix is singular",
                   "and it does not have\n",
                   "an LU decomposition.\n")
      message(msg)
      return(NULL)
    }
  }
  
  return(list(L=L,U=U,ord=rowsord))
}

#' Tridiagonal linear system
#'
#' Solution of a system of \eqn{n} equations in \eqn{n} unknowns,
#' where the coefficients form a tridiagonal matrix.
#' 
#' The linear system to solve is \eqn{Ax=b}, where \eqn{A} is
#' the \eqn{n\times n} matrix of coefficients of the \eqn{n}
#' unknowns in the \eqn{n\times 1} vector \eqn{x}, and \eqn{b}
#' is the \eqn{n\times 1} vector of known numbers. Matrix \eqn{A}
#' is a tridiagonal matrix. This means that \eqn{A} is a sparse
#' matrix with non-zero elements on the main diagonal and the
#' two diagonals adjacent to the main diagonal. The special form
#' of the matrix of coefficients makes it possible to solve the
#' related system using a fast algorithm, here the Thomas
#' algorithm.
#' 
#' @param M The \eqn{n\times (n+1)} augmented matrix of
#'        coefficients corresponding to the system of \eqn{n}
#'        linear equations in \eqn{n} unknowns, \eqn{Ax=b}. 
#' @return A vector of length \eqn{n} containing the \eqn{n}
#'         numeric solutions for the \eqn{n} unknowns. If the
#'         system has no solutions or an infinite number of
#'         solutions, the function returns NULL and dumps a
#'         warning message.
#'         
#' @examples 
#' # System of four equations in four unknowns
#' #
#' # 2x_1 +  x_2               = 1
#' # 2x_1 + 3x_2 +  x_3        = 2
#' #         x_2 + 4x_3 + 2x_4 = 3
#' #                x_3 + 3x_4 = 4
#' 
#' # Augmented matrix M=(A|b)
#' M <- matrix(c(2,2,0,0,1,3,1,0,0,1,4,1,0,0,2,3,1,2,3,4),
#'             ncol=5)
#' 
#' # Solution via Thomas algorithm
#' x <- solve_tridiag(M)
#' print(x)
#' 
#' @export
solve_tridiag <- function(M) {
  # Number of rows
  m <- nrow(M)
  
  # Number of columns
  n <- ncol(M)
  
  # One column more than # of rows
  if (n != (m+1)) {
    warning("Input is not an augmented matrix.\n")
    return(NULL)
  }
  
  # If any of the elements on the diagonal is zero,
  # the system has no solutions or infinite solutions
  flag <- FALSE
  for (i in 1:m) {
    if (abs(M[i,i]) < 1e-6) {
      flag <- TRUE
    }
    if (flag) {
      message("This system has no solution or infinite solutions.\n")
      return(NULL)
    }
  }
  
  # Save the three diagonal on alpha, beta, gamma vectors
  alpha <- c(0,diag(M[2:m,]))
  beta <- diag(M)
  gamma <- c(diag(M[1:(m-1),2:m]),0)
  b <- M[,n]
  
  # Main algorithm
  for (i in 1:(m-1)) {
    ss <- alpha[i+1]/beta[i]
    beta[i+1] <- beta[i+1]-ss*gamma[i]
    b[i+1] <- b[i+1]-ss*b[i]
  }
  
  # If any of the elements on the diagonal is zero,
  # the system has no solutions or infinite solutions
  idx <- which(abs(beta) < 1e-6)
  if (length(idx) > 0) {
    message("This system has no solution or infinite solutions.\n")
    return(NULL)
  }
  
  # Solution by back-substitution
  x <- rep(0,length=m)
  x[m] <- b[m]/beta[m]
  for (i in (m-1):1) {
   x[i] <- (b[i]-gamma[i]*x[i+1])/beta[i] 
  }
  
  return(x)
}

#' Determinant of a square matrix
#' 
#' Calculates the determinant of a square matrix of size
#' \eqn{n}, using the reduction of the matrix in upper
#' triangular form.
#' 
#' @param A An \eqn{n\times n} square matrix of. 
#' @return A real number corresponding to the determinant
#'         of A.
#'         
#' @examples 
#' # Identity matrix of size 10
#' A <- diag(10)
#' print(condet(A))
#' 
#' # Random matrix with integer elements
#' A <- matrix(sample(-5:5,size=25,replace=TRUE),ncol=5)
#' print(condet(A))
#' 
#' @export
condet <- function(A) {
  # Check the matrix is a square matrix
  m <- nrow(A)
  n <- ncol(A)
  if (n != m) {
    warning("Input is not a square matrix.\n")
    return(NULL)
  }
  
  # Transform to upper triangular
  rowsord <- 1:n
  for (i in 1:(n-1)) {
    # Swap rows to have largest values on diagonal
    idx <- which(abs(A[i:n,i]) == max(abs(A[i:n,i])))
    idx <- idx[length(idx)] # If more than one, pick the last
    idx <- i+idx-1
    N <- A[i,]
    A[i,] <- A[idx,]
    A[idx,] <- N
    idxN <- rowsord[i]
    rowsord[i] <- rowsord[idx]
    rowsord[idx] <- idxN
    
    # Next, perform row operations for column i
    for (j in (i+1):n) {
      if (abs(A[j,i]) > 1e-6) {
        cff <- -A[j,i]/A[i,i]
        A[j,] <- cff*A[i,]+A[j,]
      }
    }
  }
  
  # Number of row permutations
  ipa <- oddity(rowsord)
  
  # Calculate determinant as product of diagonal
  # elements with oddity
  ddd <- ipa*prod(diag(A))
  
  return(ddd)
}

#' Parity of a permutation
#' 
#' Given a permutation of the integers from 1 to \eqn{n}, this
#' function calculates its parity (+1 or -1), i.e. the number
#' of swapping that take the permutation back to the natural
#' ordering \eqn{1,2,\dots,n}.
#' 
#' @param x A vector containing a permutation of the first
#'          \eqn{n} integers, \eqn{1,2,\dots,n}. 
#' @return A real number equal to +1 or -1, indicating the
#'         parity of the given permutation.
#'         
#' @examples 
#' # Identity permutation (10 elements)
#' x <- 1:10
#' print(oddity(x))
#' 
#' # One swap
#' x[2] <- 5
#' x[5] <- 2
#' print(oddity(x))
#' 
#' @export
oddity <- function(x) {
  # Length of vector
  n <- length(x)
  
  # Comparing loop
  iperm <- 0
  for (i in 1:(n-1)) {
    idx <- which(x[(i+1):n] == i)
    idx <- idx+i
    if (length(idx) > 0) {
      iperm <- iperm+1
      N <- x[i]
      x[i] <- x[idx]
      x[idx] <- N
    }
  }
  
  # Parity
  ipa <- 1
  if (iperm %% 2 != 0) ipa <- -1
  
  return(ipa)
}


#' The Jacobi method
#' 
#' Implementation of the Jacobi iterative method to solve a 
#' system \eqn{Ax=b} of \eqn{n} linear equations in \eqn{n} unknowns.
#' 
#' The Jacobi method guarantees a finite solution for linear systems
#' characterised by a diagonally dominant matrix \eqn{A} of 
#' coefficients. This means that each element on its diagonal must be,
#' in absolute value, larger than the sum of the absolute value of
#' all the elements in the corresponding row. 
#' 
#' @param A The \eqn{n \times n} matrix of coefficients of the 
#'          unknowns in the linear system.
#' @param b A vector of \eqn{n} constants representing the right-hand 
#'          side of the linear system. This function does not work out
#'          solutions of homogeneous systems, where the b is a vector
#'          of zeros (null vector). Therefore input with b equal to a
#'          null vector is rejected.
#' @param x0 A vector of \eqn{n} starting numeric values for the 
#'           iterations. If no values are entered for x0, a column of 
#'           zero will be adopted by default.
#' @param tol A real number indicating the threshold under which the
#'            relative increment from one solution approximation to
#'            the next is small enough to stop iteration. The default
#'            value is \code{tol=1e-6}.
#' @param nmax An integer. The maximum number of iterations allowed,
#'             if convergence according to the criterion is not reached.
#' @param ddominant A logical variable. If FALSE, the method is applied 
#'                  also if the matrix of coefficients is not
#'                  diagonally dominant (default is TRUE).
#' @return A numeric vector of length \eqn{n} with values approximating
#'         the system's solution.
#'
#' @examples 
#' # Simple system with solution 1,2,3
#' A <- matrix(c(3,1,2,-1,-4,2,1,1,7),ncol=3)
#' b <- c(4,-4,27)
#' 
#' # Solution
#' x <- PJacobi(A,b)
#' print(x)
#' 
#' # Start from a different point
#' x0 <- c(-1,2,8)
#' x <- PJacobi(A,b,x0)
#' print(x)
#' 
#' @export
PJacobi <- function(A,b,x0=NULL,tol=1e-6,nmax=100000,
                    ddominant=TRUE) {
  # Checks
  ans <- is.matrix(A)
  if (!ans) {
    msg <- "A has to be an n X n matrix."
    warning(msg)
    return(NULL)
  }
  tmp <- dim(A)
  n <- tmp[1]
  ans <- length(b) == n
  if (!ans) {
    msg <- "b has to be a vector of length n."
    warning(msg)
    return(NULL)
  }
  ans <- sum(abs((b))) > tol
  if (!ans) {
    msg <- "b has to be different from a null vector."
    warning(msg)
    return(NULL)
  }
  if (is.null(x0)) {
    x0 <- rep(0,times=n)
  }
  ans <- length(x0) == n
  if (!ans) {
    msg <- "x0 has to be a vector of length n."
    warning(msg)
    return(NULL)
  }
  
  # Turn b and x0 into appropriate matrices
  b <- matrix(b,ncol=1)
  x0 <- matrix(x0,ncol=1)
  
  # Check that A is diagonally dominant
  
  # Pivoting if not correct order
  # Swap rows to have largest values on diagonal
  for (i in 1:(n-1)) {
    idx <- which(abs(A[i:n,i]) == max(abs(A[i:n,i])))
    idx <- idx[length(idx)] # If more than one, pick the last
    idx <- i+idx-1
    N <- A[i,]
    A[i,] <- A[idx,]
    A[idx,] <- N
    N <- b[i]
    b[i] <- b[idx]
    b[idx] <- N
  }
  
  # Important check: determinant must be non zero
  tmp <- det(A)
  if (abs(tmp) < 1e-9) {
    msg <- paste("The determinant of this matrix is zero.\n",
                 "Abandoning solution process ...\n")
    warning(msg)
    return(NULL)
  }
  
  # Each a_{ii} greater than sum of rest
  ans <- TRUE
  for (i in 1:n) {
    ff <- abs(A[i,i])
    ss <- sum(abs(A[1,]))-ff
    if (ff <= ss) ans <- FALSE
  }
  if (!ans & ddominant) {
    msg <- paste0("The matrix of coefficients is not diagonally",
                 " dominant.\nNot attempting solution.\n")
    warning(msg)
    return(NULL)
  }
  if (!ans & !ddominant) {
    msg <- paste("The matrix of coefficients is not diagonally",
                 "dominant.\nAttempting solution anyway...\n\n")
    warning(msg)
  }
  
  # Proceed
  
  # Elements on the diagonal
  vtmp <- diag(A)
  
  # If one of them is zero, stop
  ans <- abs(vtmp) < 1e-9
  if (sum(ans) != 0) {
    msg <- paste0("The matrix of coefficients appears to have\n",
                 "one of the elements on the diagonal equal to\n",
                 "zero. Algorithm cannot proceed.\n")
    warning(msg)
    return(NULL)
  }
  
  # Inverse diagonal matrix
  Dinv <- diag(1/vtmp)
  
  # Iterations
  x <- x0
  
  for (icyc in 1:nmax) {
    # Copy previous results, for increment calculation
    xold <- x
    
    # Vectorial Increment
    Deltax <- Dinv %*% (b - A %*% x)
    
    # Update values
    x <- x + Deltax
    
    # Initial increment
    if (icyc == 1) {
      Mold <- max(abs(x-x0))
      if (max(abs(x-x0) < tol)) {
        msg <- paste0("Number of cycles needed to converge: 0\n",
                     "The starting value was the solution.\n")
        message(msg)
        break
      }
    }
    
    # Start measuring convergence after the first cycle
    # just to avoid the starting solution 0
    if (icyc > 1) {
      # Measure the increment
      Mnew <- max(abs(x-xold))
      
      # The procedure is not converging (it shouldn't,
      # but you never know ...)
      if (Mnew > 1e+100) {
        msg <- paste0("The increment is getting larger and larger.\n",
                     "Max value of increment at cycle ",icyc,": ",
                     Mnew,"\n")
        message(msg)
        return(NULL)
      }
      
      eps <- abs(Mold-Mnew)/Mold
      if (abs(eps-1) < tol) {
        msg <- paste("Number of cycles needed to converge:",icyc,"\n")
        msg <- paste0(msg,"Last relative increment: ",abs(1-eps),"\n")
        message(msg)
        break
      }
    }
  }
      
  # Transform x into a numeric vector
  x <- as.vector(x)
  
  return(x)
}


#' The Gauss-Seidel algorithm
#' 
#' Implementation of the Gauss-Seidel iterative method to solve a 
#' system \eqn{Ax=b} of \eqn{n} linear equations in \eqn{n} unknowns.
#' 
#' Gauss-Seidel is a variant of the Jacobi method, as it guarantees 
#' a finite solution for linear systems characterised by a diagonally 
#' dominant matrix \eqn{A} of coefficients. This means that each 
#' element on its diagonal must be, in absolute value, larger than the 
#' sum of the absolute value of all the elements in the corresponding 
#' row. Gauss-Seidel differs from Jacobi as it promises to converge
#' faster than Jacobi. 
#' 
#' @param A The \eqn{n \times n} matrix of coefficients of the 
#'          unknowns in the linear system.
#' @param b A vector of \eqn{n} constants representing the right-hand 
#'          side of the linear system. This function does not work out
#'          solutions of homogeneous systems, where the b is a vector
#'          of zeros (null vector). Therefore input with b equal to a
#'          null vector is rejected.
#' @param x0 A vector of \eqn{n} starting numeric values for the 
#'           iterations. If no values are entered for x0, a column of 
#'           zero will be adopted by default.
#' @param tol A real number indicating the threshold under which the
#'            relative increment from one solution approximation to
#'            the next is small enough to stop iteration. The default
#'            value is \code{tol=1e-6}.
#' @param nmax An integer. The maximum number of iterations allowed,
#'             if convergence according to the criterion is not reached.
#' @param ddominant A logical variable. If FALSE, the method is applied 
#'                  also if the matrix of coefficients is not
#'                  diagonally dominant (default is TRUE).
#' @return A numeric vector of length \eqn{n} with values approximating
#'         the system's solution.
#'
#' @examples 
#' # Simple system with solution 1,2,3
#' A <- matrix(c(3,1,2,-1,-4,2,1,1,7),ncol=3)
#' b <- c(4,-4,27)
#' 
#' # Solution
#' x <- GSeidel(A,b)
#' print(x)
#' 
#' # Start from a different point
#' x0 <- c(-1,2,8)
#' x <- GSeidel(A,b,x0)
#' print(x)
#' 
#' @export
GSeidel <- function(A,b,x0=NULL,tol=1e-6,nmax=100000,
                    ddominant=TRUE) {
  # Checks
  ans <- is.matrix(A)
  if (!ans) {
    msg <- "A must be an n X n matrix."
    warning(msg)
    return(NULL)
  }
  tmp <- dim(A)
  n <- tmp[1]
  ans <- length(b) == n
  if (!ans) {
    msg <- "b must be a vector of length n."
    warning(msg)
    return(NULL)
  }
  ans <- sum(abs((b))) > tol
  if (!ans) {
    msg <- "b must be different from a null vector."
    warning(msg)
    return(NULL)
  }
  if (is.null(x0)) {
    x0 <- rep(0,times=n)
  }
  ans <- length(x0) == n
  if (!ans) {
    msg <- "x0 must be a vector of length n."
    warning(msg)
    return(NULL)
  }
  
  # Turn b and x0 into vector (in case they are not)
  b <- as.vector(b)
  x0 <- as.vector(x0)
  
  # Check that A is diagonally dominant
  
  # Pivoting if not correct order
  # Swap rows to have largest values on diagonal
  for (i in 1:(n-1)) {
    idx <- which(abs(A[i:n,i]) == max(abs(A[i:n,i])))
    idx <- idx[length(idx)] # If more than one, pick the last
    idx <- i+idx-1
    N <- A[i,]
    A[i,] <- A[idx,]
    A[idx,] <- N
    N <- b[i]
    b[i] <- b[idx]
    b[idx] <- N
  }
  
  # Important check: determinant must be non zero
  tmp <- det(A)
  if (abs(tmp) < 1e-9) {
    msg <- paste("The determinant of this matrix is zero.\n",
                 "Abandoning solution process ...\n")
    warning(msg)
    return(NULL)
  }
  
  # Each a_{ii} greater than sum of rest
  ans <- TRUE
  for (i in 1:n) {
    ff <- abs(A[i,i])
    ss <- sum(abs(A[1,]))-ff
    if (ff <= ss) ans <- FALSE
  }
  if (!ans & ddominant) {
    msg <- paste0("The matrix of coefficients is not diagonally",
                  " dominant.\nNot attempting solution.\n")
    warning(msg)
    return(NULL)
  }
  if (!ans & !ddominant) {
    msg <- paste("The matrix of coefficients is not diagonally",
                 "dominant.\nAttempting solution anyway...\n\n")
    warning(msg)
  }
  
  # If one of the items on the diagonal is zero, stop
  vtmp <- diag(A)
  ans <- abs(vtmp) < 1e-9
  if (sum(ans) != 0) {
    msg <- paste0("The matrix of coefficients appears to have\n",
                  "one of the elements on the diagonal equal to\n",
                  "zero. Algorithm cannot proceed.\n")
    warning(msg)
    return(NULL)
  }
  
  # Proceed
  
  # Vector of results
  x <- x0
  
  # Main iteration
  for (icyc in 1:nmax) {
    # Copy previous results, for increment calculation
    xold <- x
    
    # Just one loop as the operations over j
    # are taken care of with the "-i" trick
    for (i in 1:n) {
      x[i] <- b[i]/A[i,i] - sum((A[i,-i]/A[i,i])*x[-i])
    }
    
    # Initial increment
    if (icyc == 1) {
      Mold <- max(abs(x-x0))
      if (Mold < tol) {
        msg <- paste0("Number of cycles needed to converge: 0\n",
                      "The starting value was the solution.\n")
        message(msg)
        break
      }
    }
    
    # Stop when increment is less than tolerance (after cycle 1)
    if (icyc > 1) {
      Deltax <- abs(x-xold)
      Mnew <- max(Deltax)
      
      # The procedure is not converging (it shouldn't,
      # but you never know ...)
      if (Mnew > 1e+100) {
        msg <- paste0("The increment is getting larger and larger.\n",
                     "Max value of Delta_X at cycle ",icyc,": ",
                     Mnew,"\n")
        message(msg)
        return(NULL)
      }
      
      eps <- abs(Mold-Mnew)/Mold
      if (abs(eps-1) < tol) {
        msg <- paste("Number of cycles needed to converge:",icyc,"\n")
        msg <- paste0(msg,"Last relative increment: ",abs(1-eps),"\n")
        message(msg)
        break
      }
    }
  }
  
  return(x)
}


#' Ill-conditioned sampling
#' 
#' Random sampling, based on the uniform distribution, of the 
#' right-hand side, \eqn{b}, of a linear system and of a perturbation, 
#' \eqn{\Delta b}, so that the solution of \eqn{Ax=b} is very 
#' different from the solution of \eqn{Ax=b+\Delta b}.
#' 
#' The degree of ill-conditioning of a system is not only measured 
#' by the matrix's condition number, but also from the solution
#' relative error. If \eqn{\Delta x} is the difference between the
#' solution, \eqn{x}, of the system related to \eqn{b} and the
#' solution, \eqn{x'}, of the system related to \eqn{b'= b-\Delta b},
#' then the ratio of the norm of \eqn{\Delta x} and the norm of
#' \eqn{x}, is the solution relative error. Norms are Frobenius norms. 
#' This function returns a named list with b and Db the chosen 
#' \eqn{b} and \eqn{\Delta b}, based on random sampling of a 
#' specified region.
#' 
#' @param A The \eqn{n \times n} matrix of coefficients of the 
#'          unknowns in the linear system.
#' @param bmax A numeric number providing the interval, (0,bmax),
#'             in which the \eqn{n} uniformly random components of
#'             \eqn{b} are selected. Default value is 100.
#' @param Dbmax A numeric number providing the interval, (0,Dbmax),
#'             in which the \eqn{n} uniformly random components of
#'             \eqn{\Delta b} are selected. Default value is 1.
#' @param ncyc An integer indicating the number of uniform random
#'             selection of the \eqn{n} components of \eqn{b} and the
#'             \eqn{n} components of \eqn{\Delta b}. The higher this
#'             number, the higher the chance of getting a high 
#'             relative solution number, but the longer the execution
#'             time of the function. Default is 100000.
#' @param iseed An integer. The seed starting random generation. If
#'              a value is provided, the (pseudo-)random generation 
#'              will reproduce exactly the same \eqn{b} and
#'              \eqn{\Delta b}. Default is NULL, which means that the
#'              seed will be randomly chosen at every execution of 
#'              the function.
#' @return A named list with names b, a vector equal of the right-hand
#'         side of the linear system, and Db, a vector equal to the
#'         perturbations, \eqn{\Delta b}, to be applied to \eqn{b}.
#'
#' @examples 
#' # This is a simple but ill-conditioned matrix
#' A <- matrix(c(2,1,1.99,1),ncol=2)
#' 
#' # Select b and Db randomly, starting with iseed=2341
#' ltmp <- illcond_sample(A,iseed=2341)
#' names(ltmp)
#' 
#' # b and b'
#' b <- ltmp$b
#' Db <- ltmp$Db
#' b2 <- b-Db
#' 
#' # Solution for b
#' x <- solve(A,b)
#' print(x)
#' 
#' # Solution for b'
#' x2 <- solve(A,b2)
#' print(x2)
#' 
#' # Difference
#' Dx <- x-x2
#' 
#' # Solution relative error (Frobenius norm)
#' print(norm(Dx,"F")/norm(x,"F"))
#' 
#' # Upper limit
#' Ainv <- solve(A)
#' print(norm(A,"F")*norm(Ainv,"F")*norm(Db,"F")/norm(b,"F"))
#' 
#' @export
illcond_sample <- function(A,
                           bmax=100,Dbmax=1,
                           ncyc=100000,
                           iseed=NULL) {
  # Size of A
  tmp <- dim(A)
  n <- tmp[1]
  
  # Main loop
  x0 <- matrix(nrow=n,ncol=1)
  b0 <- matrix(nrow=n,ncol=1)
  Dx0 <- matrix(nrow=n,ncol=1)
  Db0 <- matrix(nrow=n,ncol=1)
  normx <- 1e+9
  normDx <- 0
  if (!is.null(iseed)) set.seed(iseed)
  for (i in 1:ncyc) {
    b <- matrix(stats::runif(n,min=0,max=bmax),ncol=1)
    Db <- matrix(stats::runif(n,min=0,max=Dbmax),ncol=1)
    x <- solve(A,b)
    Dx <- solve(A,Db)
    if (norm(x,"F") < normx) {
      x0 <- x
      b0 <- b
      normx <- norm(x,"F")
    }
    if (norm(Dx,"F") > normDx) {
      Dx0 <- Dx
      Db0 <- Db
      normDx <- norm(Dx,"F")
    }
  }
  
  # Value and max
  Ainv <- solve(A)
  ulim <- norm(Ainv,"F")*norm(A,"F")*norm(Db0,"F")/norm(b0,"F")
  relerr <- norm(Dx0,"F")/norm(x0,"F")
  msg <- paste0("Relative error: ",relerr," < ",ulim,": upper limit.\n")
  message(msg)
  msg <- paste0("Ratio: ",relerr/ulim,"\n")
  message(msg)
  
  return(list(b=b0,Db=Db0))
}

###--------------------------------------------------------------
### Auxiliary functions (not for users
###--------------------------------------------------------------
normJ <- function(matx,type="I") {
  # Flag for larger size
  flag <- TRUE
  
  # Number of columns
  tmp <- dim(matx)
  nr <- tmp[1]
  nc <- tmp[2]
  n <- nr
  if (nc > nr) flag <- FALSE
  if (!flag) n <- nc
  
  nrms <- c()
  for (i in 1:n) {
    if (flag) A <- matrix(matx[i,],ncol=1)
    if (!flag) A <- matrix(matx[,i],ncol=1)
    nrms <- c(nrms,norm(A,type=type))
  }
  
  return(nrms)
}