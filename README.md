# Computational Physics
Course _Computational Physics_ from NTU.</br>
I will use c language.</br>
Most of the content will be my note and some external links I think it's useful.

## Introduction
* Data Representation (IEEE Standard)</br>
How computer store number.
  * Integer
  * Single precision real number
  * Double precision real number
* Machine Precision </br>
Smallest number epsilon that the computer cannot tell their difference. There will be different epsilon for different data type.
  * [machine_precision.c](/Introduction/machine_precision.c)

## Differentiation
* First order
  * Use Taylor expansion to expand f(x) with different neighborhood, so that it can canceled out the higher order term.
* Second order
  * Three point formula
* Richardson Extrapolation </br>
The iteration steps relaiton can be used at first or second order.
  * [first_order_diff.c](/Differentiation/first_order_diff.c)

## Integration
* We can cut f(x) with fixed interval on x-axis, then use different degrees of polynomial to find the weight of f(x). So we can sum all pieces of area.</br>
With this, we can organized it as a rule or formula.
  * Trapezoidal Rule (Piecewise Linear)
  * Simposon's Rule (Quadratic)
  * 4-points formula
  * 5-points formula (Quartic)
* Romberg Integration </br>
The concepts are the same as Richardson Extrapolation, we try to cancel out the error. So the more iteration we made, the more accurate the approximation is.