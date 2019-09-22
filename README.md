# Computational Physics
Course _Computational Physics_ from NTU.</br>
I will use c language.</br>
Most of the content will be my note and some external links I think it's useful.

## Introduction
* Data Representation (IEEE Standard)
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