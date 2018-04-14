project: ddeabm
project_dir: ./src
output_dir: ./doc
project_github: https://github.com/jacobwilliams/ddeabm
summary: Modern Fortran Implementation of the DDEABM Adams-Bashforth-Moulton ODE Solver
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         private
         protected
source: true
graph: true
media_dir: media
exclude: pyplot_module.f90
         test.f90
         zeroin_test.f90
         performance_test.f90
         ddeabm_test_vec.f90
         ddeabm_test.f90
         ddeabm_stepsize_test.f90
         ddeabm_example.f90
         ddeabm_fixed_step_test.f90
extra_mods: pyplot_module:https://github.com/jacobwilliams/pyplot-fortran
            iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

Brief description
---------------

This is a modern object-oriented Fortran implementation of the DDEABM Adams-Bashforth-Moulton ODE solver. The original Fortran 77 code was obtained from the [SLATEC library](http://www.netlib.org/slatec/src/). It has been extensively refactored.

DDEABM uses the Adams-Bashforth-Moulton predictor-corrector formulas of orders 1 through 12 to integrate a system of first order ordinary differential equations of the form `dx/dt = f(t,x)`. Also included is an event-location capability, where the equations can be integrated until a specified function `g(t,x) = 0`.

This project is hosted on [GitHub](https://github.com/jacobwilliams/ddeabm).

## Keywords

Adams-Bashforth-Moulton Method, DEPAC, Initial Value Problems,
ODE, Ordinary Differential Equations, Predictor-Corrector

## References

1. L. F. Shampine, M. K. Gordon, "Solving ordinary differential equations with ODE, STEP, and INTRP",  Report SLA-73-1060, Sandia Laboratories, 1973.
2. L. F. Shampine, M. K. Gordon, "Computer solution of ordinary differential equations, the initial value problem", W. H. Freeman and Company, 1975.
3. L. F. Shampine, H. A. Watts, "DEPAC - Design of a user oriented package of ode solvers", Report SAND79-2374, Sandia Laboratories, 1979.
4. H. A. Watts, "A smoother interpolant for DE/STEP, INTRP and DEABM: II", Report SAND84-0293, Sandia Laboratories, 1984.
5. R. P. Brent, "[An algorithm with guaranteed convergence for finding a zero of a function](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)", The Computer Journal, Vol 14, No. 4., 1971.
6. R. P. Brent, "[Algorithms for minimization without derivatives](http://maths-people.anu.edu.au/~brent/pub/pub011.html)", Prentice-Hall, Inc., 1973.

# License

The ddeabm source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/ddeabm/blob/master/LICENSE) (BSD-style).  The original DDEABM Fortran 77 code is [public domain](http://www.netlib.org/slatec/guide).
