### ---- TEST CAS LINÉAIRE -----

#Coefficients de f(x) = 3x1 + 2x2
a <- c(6, 4)

b0 <- 0   # constante (pas utilisée ici)

# Matrice A pour Ax <= b
A <- matrix(c(3,9,
              4,5,
              2,1,
             -1,0,
              0,-1), nrow=5, byrow=TRUE)

# Vecteur b
b <- c(81, 55, 20, 0, 0)


# Point initial
x0 <- c(1, 1)

# Résolution
solution_linear <- solve_PI(
  type = "linear",
  a = a,
  b0 = b0,
  A = A,
  b = b,
  x0 = x0,
  mu = 1,
  mode = "max"
)

solution_linear


 ### ---- TEST CAS QUADRATIQUE -----

# Matrice A
A <- matrix(c(1, 1,
              -1, 0,
              0, -1), nrow=3,, byrow=TRUE)

# Vecteur b
b <- c(5, 0, 0)

# Point initial strictement intérieur
x0 <- c(1, 1)

# Résolution
solution_quadratic <- solve_PI(
  type = "quadratic",
  A = A,
  b = b,
  x0 = x0,
  mu = 1,
  mode = "min"
)

solution_quadratic
