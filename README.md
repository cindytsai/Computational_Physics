# Computational Physics

Course _Computational Physics_ from NTU.</br>
I will use c language.</br>
Most of the content will be my note and some external links I think it's useful.

## Contents
* ![Assignments](https://github.com/cindytsai/Computational_Physics#assignments)

|     Folder    |                                                                  Tags                                                                  |
|:-------------:|:--------------------------------------------------------------------------------------------------------------------------------------:|
| Problem Set 1 |                                                   Machine Precision, Differentiation                                                   |
| Problem Set 2 |                                                               Integration                                                              |
| Problem Set 3 |                                                  Random Walk, Monte Carlo Integration                                                  |
| Problem Set 4 |                                                             Estimate Error                                                             |
| Problem Set 5 |                                                             2D Ising Model                                                             |
| Problem Set 6 |                                     LU-Decomposition, Inverse Matrix, Conjugate Gradient Algorithm                                     |
| Problem Set 7 |                                        Solve ODE, Physical Pendulum, Real Time Graphics(OpenGL)                                        |
| Problem Set 8 | Solve PDE, Vibrating String, Heat Diffusion, Real Time Graphics(OpenGL), Poisson Equation|
| Problem Set 9 |                             1D Real Scalar Field, 1D Scalar Field with &#955;&#966;<sup>4</sup> Interaction                            |

* ![Lectures](https://github.com/cindytsai/Computational_Physics#lectures)

## Assignments

Problem set questions are named with `psX_2019.pdf`, my reports are named with `Problem Set X.pdf`.

### Problem Set 1
* Machine Precision
* Differentiation
  * Richardson Extapolation

### Problem Set 2
* Integration
  * Trapezoidal Method
  * Simpson Method
  * 5-points Formula
  * Romberg Integration

### Problem Set 3
* Random walk on a 2D-lattice
* Monte Carlo Integration in 10-dim
  * Simple Sampling
  * Rejection Method
  * Metropolis Algorithm

### Problem Set 4
* Estimate the error in Monte Carlo Integration
  * Integrated autocorrelation time
  * Binning Method
  * Jackknife Method
  * Binning with the Jackknife Method

### Problem Set 5
* Single Cluster Algorithm
  * Prove that it satisfies detailed balance.
* 2D Ising Model
  * Using Methods
    * Metropolis Algorithm
    * Heat Bath
    * Single Cluster Algorithm
  * Measure the expectation values
    * Energy density 
    * Specific heat
    * Magnetization density
    * Magnetic susceptibility
  * Estimate the error of the expectation values
    * Jackknife Method
  * Measure the integrated autocorrelation time
    * Compute dynamical critical exponent
  * Find critical temperature

### Problem Set 6
* LU-Decomposition with Pivoting
* Inverse of a matrix using LU-Decomposition
* Conjugate Gradient Algorithm
  * Prove that conjugate gradient algorithm also works for A = D<sup>&#8224;</sup>D.
* Solve complex linear system with conjugate gradient algorithm

### Problem Set 7
* Solve ODE
  * Euler Method
  * Modified Euler Method
  * Improved Euler Method
  * 3rd order Runge-Kutta Method
* 3rd order and 4th order Runge-Kutta Method
  * Prove the weight coefficient in 3rd order Runge-Kutta Method
  * Prove the weight coefficient in 4th order Runge-Kutta Method
  _These are wrong..._
* Physical Pendulum
  * Solve 2nd order differential equation of this system.
    * Adaptive 4th order Runge-Kutta Method
  * Animate its motion in real time graphics.

### Problem Set 8
* Vibrating String with gravity
  * Solve PDE of this system with fixed end.
  * Animate its motion with real time graphics
  * Initial Condition
    * Sinusoidal Wave
    * Plucked String
    * Gaussian Wave Packet
* Heat Diffusion
  * Solve for the temperature distribution on a circular plate.
  * Animate the heat diffusion with real time graphics.
* Poisson Equation
  * Jacobi Method
  * Gauss Seidel Method
  * Successive Over-Relaxation
  * Conjugate Gradient Algorithm with Double Precision
  * Conjugate Gradient Algorithm with Mixed Precision

### Problem Set 9
* 1D Real Scalar Field
  * Methropolis Algorithm
  * Hybrid Monte Carlo Method
* 1D Scalar Field with &#955;&#966;<sup>4</sup> Interaction

## Lecture
> TODO: Match the lecture topic to specific problem sets.

### Introduction
* Data Representation (IEEE Standard)</br>
How computer store number.
  * Integer
  * Single precision real number
  * Double precision real number
* Machine Precision </br>
Smallest number epsilon that the computer cannot tell their difference. There will be different epsilon for different data type.
  * [machine_precision.c](/Introduction/machine_precision.c)

### Differentiation
* First order
  * Use Taylor expansion to expand f(x) with different neighborhood, so that it can canceled out the higher order term.
* Second order
  * Three point formula
* Richardson Extrapolation </br>
The iteration steps relaiton can be used at first or second order.
  * [first_order_diff.c](/Differentiation/first_order_diff.c)

### Integration
* We can cut f(x) with fixed interval on x-axis, then use different degrees of polynomial to find the weight of f(x). So we can sum all pieces of area.</br>
With this, we can organized it as a rule or formula.
  * Trapezoidal Rule (Piecewise Linear)
  * Simposon's Rule (Quadratic)
  * 4-points formula
  * 5-points formula (Quartic)
* Romberg Integration </br>
The concepts are the same as Richardson Extrapolation, we try to cancel out the error. So the more iteration we made, the more accurate the approximation is.
