# Java Finite Element Library 
[![tests](https://github.com/Tovermodus/JFEM/actions/workflows/tests.yml/badge.svg)](https://github.com/Tovermodus/JFEM/actions/workflows/tests.yml)

JFEM is a Java library for the computational solution of partial differential equations using finite elements.
It has the functionality to approximate the solution to many partial differential equations with a range of boundary conditions, including but not limited to:
* the Poisson equation
* the heat equation
* the linear elasticity equation
* the simplified Stokes equations
* the time dependent Stokes equations
* the Darcy equation

We can solve these by a range of different finite elements, including:
* Qk-elements
* Discontinuous Qk-elements
* Raviart-Thomas elements
* Taylor-Hood elements
* Mixed elements

The linear system that results from the discretization can be solved by the following methods:
* LU-Decomposition
* (Preconditioned) CG method
* (Preconditioned) GMRes method
* (Preconditioned) BiCGStab method

## Getting started
This library is strongly object oriented. This makes the implementation of new functionality straightforward and also gives a straightforward interpretation to the final programs. To solve a partial differential equation, such as

<img src="http://www.sciweavers.org/tex2img.php?eq=%20-%5CDelta%20u%20%2B%205u%20%3D%20f&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt=" -\Delta u + 5u = f" width="111" height="19" />

with zero Dirichlet boundary conditions in weak form, which would be

<img src="http://www.sciweavers.org/tex2img.php?eq=%5Cint_%5COmega%20%5Cnabla%20u%5Ccdot%5Cnabla%20v%5C%2C%20dx%20%2B%20%5Cint_%5COmega5%20uv%5C%2Cdx%20%3D%20%5Cint_%5COmega%20fv%5C%2Cdx%5Cqquad%5Cforall%20v%5Cin%20H_0%5E1%28%5COmega%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\int_\Omega \nabla u\cdot\nabla v\, dx + \int_\Omega5 uv\,dx = \int_\Omega fv\,dx\qquad\forall v\in H_0^1(\Omega)" width="425" height="46" />

we first choose a finite element space such as the Qk-space of scalar valued functions. To do this in code we create an instance of the `ContinuousTPFESpace` (different spaces have different names). On this instance we then call `assembleCells()`, `assembleFunctions(<polynomial degree>)`, `initializeSystemMatrix()` and `initializeRhs()` in this order. Now we need to translate the equation into code. We do this by creating three Integral objects, one for each integral. Two of the integrals for this equation go into the system matrix and one into the right hand side, so we create two instances of `TPCellIntegral` once with the argument `TPCellIntegral.GRAD_GRAD` to represent the first integral and once with the arguments `5.0,TPCellIntegral.VALUE_VALUE` for the second integral. For the right hand side, we create an instance of `TPRightHandSideIntegral` to which we pass the function f which could for example be an anonymous class extending `ScalarFunction` and the argumment `TPRightHandSideIntegral.VALUE`. We then call the function `evaluateCellIntegrals` on our `ContinuousTPFESpace` object, to which we pass two lists , first a list containing the first two integrals, then a list containing only the third integral. If we also had integrals that are evaluated on faces we would also call `evaluateFaceIntegrals`.
Finally we obtain the solution vector by calling `space.getSystemMatrix().solve(space.getRhs())` by LU decomposition.

For more information on how to do this and further examples, see the Examples Module.

## Example Results
1. A very compressible ball that is moving to the top left is deformed by a fluid that flows in from the left: <img src=ball_in_water.gif/>
2. Two heat sources heat a plate that is cooled on all four sides:
