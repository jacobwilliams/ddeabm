project: ddeabm
project_dir: ./src
output_dir: ./doc
project_github: https://github.com/jacobwilliams/ddeabm
summary: Modern Fortran implementation of the DDEABM Adams-Bashforth algorithm
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
display: private
display: protected
source: true
graph: true

Brief description
---------------

The goal of this project is to produce a modern object-oriented Fortran implementation of the DDEABM Adams-Bashforth ODE solver.

The original Fortran 77 code is public domain, and was obtained from the [SLATEC library](http://www.netlib.org/slatec/src/). It has been extensively refactored.

*This is a work in progress.*

## References

1.  L. F. Shampine, M. K. Gordon, "Solving ordinary differential equations with ODE, STEP, and INTERP",  Report SLA-73-1060, Sandia Laboratories, 1973.
2.  L. F. Shampine, M. K. Gordon, "Computer solution of ordinary differential equations, the initial value problem", W. H. Freeman and Company, 1975.
3. L. F. Shampine, H. A. Watts, "DEPAC - Design of a user oriented package of ode solvers", Report SAND79-2374, Sandia Laboratories, 1979.
4.  H. A. Watts, "A smoother interpolant for DE/STEP, INTRP and DEABM: II", Report SAND84-0293, Sandia Laboratories, 1984.
