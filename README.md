
<!-- README.md is generated from README.Rmd. Please edit that file -->

*comphy* is an R package accompanying the book This package accompanies
the book  
[*Computational Physics with
R*](https://iopscience.iop.org/book/mono/978-0-7503-2632-2). It collects
all the functions described and used throughout the book, providing
ready-to-use implementations of numerical algorithms in computational
physics.

The package includes routines for:

- Numerical differentiation and integration

- Solvers for ordinary differential equations

- Sturm–Liouville eigenvalue problems

- Monte Carlo methods

- Other algorithms relevant to physics and applied mathematics

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("jfoadi/comphy")
```

## Example

This is a simple example where a numerical derivative of sin(x) is
computed using comphy:

``` r
library(comphy)

# Grid of known values
x <- seq(0,pi,length.out=11)

# Grid of values at which derivative must be calculated
# Same as x
x0 <- x

# Function known at tabulated points
f <- sin(x)+cos(x)

# Obtain a numerical derivative at those points
# using central differences (discard first and last value)
df <- deriv_reg(x0,x,f)

# Print derivative values
print(df)
#>  [1]         NA  0.6315304  0.2176105 -0.2176105 -0.6315304 -0.9836316
#>  [7] -1.2394482 -1.3739389 -1.3739389 -1.2394482         NA
```

The two `NA`’s are there because the derivative is calculated using a
centred difference.

## Documentation

All exported functions include documentation and usage examples. For a
comprehensive discussion and context, see the book [*Computational
Physics with
R*](https://iopscience.iop.org/book/mono/978-0-7503-2632-2).

## License

GPL (\>= 2)
