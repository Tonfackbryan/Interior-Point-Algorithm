# Interior-Point-Algorithm(R Package)
R package implementing the interior point method for solving constrained linear and quadratic optimization problems using logarithmic barrier function.

## Overview

This project implements a basic Interior Point Method (IPM) for solving
constrained optimization problems with inequality constraints.

The algorithm is designed for educational and academic purposes and
illustrates how interior point techniques can be implemented from scratch
in R using a logarithmic barrier approach.

The method ensures that all iterates remain strictly inside the feasible
region throughout the optimization process.

---

## Features

- Supports constrained optimization problems with inequality constraints
- Uses a logarithmic barrier method
- Handles both minimization and maximization problems
- Implements linear and quadratic objective functions
- Simple and readable implementation
- Designed for learning and experimentation

---

## Supported Optimization Problems

### Linear Optimization

The package can solve linear optimization problems where the objective
function is a linear combination of the decision variables, subject to
inequality constraints.

### Quadratic Optimization

The package also supports quadratic optimization problems with a quadratic
objective function, subject to inequality constraints.

---

## Algorithm Description

The interior point method transforms a constrained optimization problem
into a sequence of unconstrained problems by introducing a barrier term.

At each iteration:
- The gradient of the objective function is computed
- A barrier term penalizes approaches to the constraint boundaries
- A descent direction is calculated
- A step size is selected to preserve feasibility
- The barrier parameter is gradually reduced

This process continues until convergence.

---

## Installation

Clone the repository from GitHub:

```bash
git clone https://github.com/your-username/interior-point-method.git
