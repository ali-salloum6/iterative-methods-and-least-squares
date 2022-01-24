# Iterative methods and least squares (linear algebra)
Date of completion: 12 Apr 2021.

This project is coursework for the linear algebra course at Innopolis University. It is based on a previous project which you can check [here](https://github.com/ali-salloum6/Linear-systems-calculator "GitHub repository").

You can find the description of each problem in its directory in addition to an application of the least-squares with a [chart](Least%20squares%20approximation/Chart/report%20(code%20and%20chart).pdf) of a quadratic regression of given data after applying least-squares approximation using Gnuplot with C++.

The project consists of three main parts:
- Jacobi method.
- Seidel method.
- Least-squares approximation.

#### Jacobi method:
In numerical linear algebra, the Jacobi method is an iterative algorithm for determining the solutions of a strictly diagonally dominant system of linear equations. Each diagonal element is solved for, and an approximate value is plugged in. The process is then iterated until it converges. [(Wikipedia)](https://en.wikipedia.org/wiki/Jacobi_method)

#### Seidel method:

In numerical linear algebra, the Gauss-Seidel method, also known as the Liebmann method or the method of successive displacement, is an iterative method used to solve a system of linear equations. It is named after the German mathematicians Carl Friedrich Gauss and Philipp Ludwig von Seidel, and is similar to the Jacobi method. Though it can be applied to any matrix with non-zero elements on the diagonals, convergence is only guaranteed if the matrix is either strictly diagonally dominant, or symmetric and positive definite.[(Wikipedia)](https://en.wikipedia.org/wiki/Gauss–Seidel_method)

#### Least squares:

The method of least squares is a standard approach in regression analysis to approximate the solution of overdetermined systems (sets of equations in which there are more equations than unknowns) by minimizing the sum of the squares of the residuals (a residual being: the difference between an observed value, and the fitted value provided by a model) made in the results of each individual equation.[(Wikipedia)](https://en.wikipedia.org/wiki/Least_squares)

