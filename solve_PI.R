# ===============================================================
# INTERIOR POINT METHOD PACKAGE
# ===============================================================
# This script implements a basic interior point algorithm for
# solving constrained optimization problems with inequality
# constraints of the form:
#
#   A x <= b
#
# Two objective functions are supported:
#   - Linear:     f(x) = a^T x + b0
#   - Quadratic:  f(x) = x1^2 + x2^2
#
# The algorithm uses a logarithmic barrier approach.
# ===============================================================


# ---------------------------------------------------------------
# FUNCTION 1 : INITIALIZATION OF THE INTERIOR POINT ALGORITHM
# ---------------------------------------------------------------

#' Initialize the interior point method
#'
#' @description
#' Checks feasibility of the initial point and initializes
#' the barrier parameter and Lagrange multipliers.
#'
#' @param x0 Initial strictly interior point (vector).
#' @param mu Initial barrier parameter (positive scalar).
#' @param A Constraint matrix.
#' @param b Constraint vector.
#'
#' @return A list containing the initial point, barrier parameter
#' and Lagrange multipliers.
#'
interior_init <- function(x0, mu, A, b) {

  # Check strict positivity (required by logarithmic barrier)
  if(any(x0 <= 0))
    stop("Initial point x0 must be strictly positive.")

  # Check strict feasibility of inequality constraints
  if(any(A %*% x0 >= b))
    stop("Initial point does not satisfy A x < b.")

  # Number of inequality constraints
  m <- nrow(A)

  # Initialize Lagrange multipliers (theoretical KKT completeness)
  lambda0 <- rep(1, m)

  # Return initialized variables
  return(list(
    x = x0,       # primal variable
    mu = mu,      # barrier parameter
    lambda = lambda0
  ))
}


# ---------------------------------------------------------------
# FUNCTION 2 : LOGARITHMIC BARRIER FUNCTION
# ---------------------------------------------------------------

#' Logarithmic barrier function for inequality constraints
#'
#' @description
#' Implements the logarithmic barrier associated with constraints
#' of the form A x <= b.
#'
#' @param x Current point.
#' @param A Constraint matrix.
#' @param b Constraint vector.
#'
#' @return Value of the barrier function.
#'
phi_barriere <- function(x, A, b) {

  # If constraints are violated, barrier is infinite
  if(any(A %*% x >= b))
    return(Inf)

  # Log-barrier penalizing proximity to constraint boundaries
  return(-sum(log(b - A %*% x)))
}


# ---------------------------------------------------------------
# FUNCTION 3 : GRADIENT OF THE BARRIER FUNCTION
# ---------------------------------------------------------------

#' Gradient of the logarithmic barrier
#'
#' @description
#' Computes the gradient of the barrier function:
#'   - sum log(b - A x)
#'
#' @param x Current point.
#' @param A Constraint matrix.
#' @param b Constraint vector.
#'
#' @return Gradient vector of the barrier function.
#'
grad_barriere <- function(x, A, b) {

  # Computes reciprocal slack variables
  d <- 1 / (b - A %*% x)

  # Gradient: Aᵀ * (b - A x)^(-1)
  return(t(A) %*% d)
}


# ---------------------------------------------------------------
# FUNCTION 4 : STEP SIZE SELECTION (BACKTRACKING)
# ---------------------------------------------------------------

#' Step size selection by backtracking
#'
#' @description
#' Finds a step size that keeps the iterate strictly feasible.
#'
#' @param x Current point.
#' @param direction Descent direction.
#' @param A Constraint matrix.
#' @param b Constraint vector.
#' @param alpha Reduction factor (0 < alpha < 1).
#'
#' @return Step size.
#'
choose_step <- function(x, direction, A, b, alpha = 0.5) {

  # Initial step length
  t <- 1

  # Reduce step until feasibility is preserved
  while(any(A %*% (x + t * direction) >= b)) {
    t <- t * alpha
  }

  return(t)
}


# ---------------------------------------------------------------
# FUNCTION 5 : BARRIER PARAMETER UPDATE
# ---------------------------------------------------------------
#' Update the barrier parameter
#'
#' @description
#' Reduces the barrier parameter to progressively approach
#' the constrained optimum.
#'
#' @param mu Current barrier parameter.
#' @param factor Reduction factor (0 < factor < 1).
#'
#' @return Updated barrier parameter.
#'
update_mu <- function(mu, factor = 0.3) {

  # Reduce influence of barrier over iterations
  return(mu * factor)
}


# ---------------------------------------------------------------
# FUNCTION 6 : LINEAR OPTIMIZATION PROBLEM
# ---------------------------------------------------------------

#' Solve a linear optimization problem using interior point method
#'
#' @description
#' Solves:
#'   min/max f(x) = a^T x + b0
#' subject to:
#'   A x <= b
#'
#' @param a Coefficient vector of the linear objective.
#' @param b0 Constant term (unused).
#' @param A Constraint matrix.
#' @param b Constraint vector.
#' @param x0 Initial strictly interior point.
#' @param mu Initial barrier parameter.
#' @param tol Convergence tolerance.
#' @param max_iter Maximum number of iterations.
#' @param mode "min" or "max".
#'
#' @return Approximate optimal solution vector.
#'
solve_linear_PI <- function(a, b0, A, b, x0, mu = 1,
                            tol = 1e-6, max_iter = 50,
                            mode = c("min", "max")) {

  # Select optimization mode
  mode <- match.arg(mode)

  # Convert maximization into minimization
  sign <- if (mode == "min") 1 else -1

  # Initialization
  init <- interior_init(x0, mu, A, b)
  x <- init$x
  mu <- init$mu

  for(k in 1:max_iter) {

    # Gradient of linear objective
    grad_f <- a

    # Gradient of barrier
    grad_phi <- grad_barriere(x, A, b)

    # Penalized objective gradient
    grad_total <- sign * grad_f + mu * grad_phi

    # Steepest descent direction
    direction <- -grad_total

    # Feasible step size
    step <- choose_step(x, direction, A, b)

    # Update variable
    x_new <- x + step * direction

    # Stopping criterion
    if(sum(abs(x_new - x)) < tol) break

    # Update iterate and barrier
    x <- x_new
    mu <- update_mu(mu)
  }

  return(x)
}


# ---------------------------------------------------------------
# FUNCTION 7 : QUADRATIC OPTIMIZATION PROBLEM
# ---------------------------------------------------------------

#' Solve a quadratic optimization problem using interior point method
#'
#' @description
#' Solves:
#'   min/max f(x) = x1^2 + x2^2
#' subject to:
#'   A x <= b
#'
#' @param A Constraint matrix.
#' @param b Constraint vector.
#' @param x0 Initial strictly interior point.
#' @param mu Initial barrier parameter.
#' @param tol Convergence tolerance.
#' @param max_iter Maximum number of iterations.
#' @param mode "min" or "max".
#'
#' @return Approximate optimal solution vector.
#'
solve_quadratic_PI <- function(A, b, x0, mu = 1,
                               tol = 1e-6, max_iter = 50,
                               mode = c("min", "max")) {

  # Optimization mode
  mode <- match.arg(mode)
  sign <- if (mode == "min") 1 else -1

  # Initialization
  init <- interior_init(x0, mu, A, b)
  x <- init$x
  mu <- init$mu

  for(k in 1:max_iter) {

    # Gradient of quadratic objective
    grad_f <- 2 * x

    # Gradient of barrier
    grad_phi <- grad_barriere(x, A, b)

    # Total gradient
    grad_total <- sign * grad_f + mu * grad_phi

    # Descent direction
    direction <- -grad_total

    # Step size selection
    step <- choose_step(x, direction, A, b)

    # Update
    x_new <- x + step * direction

    # Convergence check
    if(sum(abs(x_new - x)) < tol) break

    x <- x_new
    mu <- update_mu(mu)
  }

  return(x)
}


# ---------------------------------------------------------------
# FUNCTION 8 : MAIN INTERFACE FUNCTION
# ---------------------------------------------------------------

#' Solve an optimization problem using the interior point method
#'
#' @description
#' High-level interface that dispatches to either the linear or
#' quadratic interior point solver.
#'
#' @param type Problem type: "linear" or "quadratic".
#' @param a Linear objective coefficients.
#' @param b0 Constant term (optional).
#' @param A Constraint matrix.
#' @param b Constraint vector.
#' @param x0 Initial strictly interior point.
#' @param mu Initial barrier parameter.
#' @param mode Optimization mode: "min" or "max".
#'
#' @return Approximate optimal solution vector.
#'
solve_PI <- function(type = c("linear", "quadratic"),
                     a = NULL, b0 = NULL,
                     A, b, x0, mu = 1,
                     mode = c("min", "max")) {

  type <- match.arg(type)
  mode <- match.arg(mode)

  if(type == "linear") {
    return(solve_linear_PI(a, b0, A, b, x0, mu, mode = mode))
  }

  if(type == "quadratic") {
    return(solve_quadratic_PI(A, b, x0, mu, mode = mode))
  }
}
