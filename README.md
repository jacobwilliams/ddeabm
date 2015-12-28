ddeabm
======

# Status

[![Build Status](https://img.shields.io/travis/jacobwilliams/bspline-fortran/master.svg?style=plastic)](https://travis-ci.org/jacobwilliams/bspline-fortran)

# Description

This is a modern object-oriented Fortran implementation of the DDEABM Adams-Bashforth-Moulton ODE solver. The original Fortran 77 code was obtained from the [SLATEC library](http://www.netlib.org/slatec/src/). It has been extensively refactored.

DDEABM uses the Adams-Bashforth-Moulton predictor-corrector formulas of orders 1 through 12 to integrate a system of first order ordinary differential equations of the form `dx/dt = f(t,x)`. Also included is an event-location capability, where the equations can be integrated until a specified function `g(t,x) = 0`.

This project is hosted on [GitHub](https://github.com/jacobwilliams/ddeabm).

# Examples

The `ddeabm_module` provides a thread-safe and object-oriented interface to the DDEABM method. Some example use cases are presented below:

## Basic integration

This example shows how to integrate a conic orbit (6 state equations) around the Earth from an initial time `t0` to a final time `tf`:

```Fortran
program ddeabm_example

use ddeabm_module
use kind_module

implicit none

real(wp),parameter :: mu = 398600.436233_wp !! Earth gravitational parameter (km^3/s^2)
integer,parameter  :: n = 6                 !! number of state variables

type(ddeabm_class) :: s
real(wp),dimension(n) :: x0,x
real(wp) :: t0,tf,t
integer :: idid

call s%initialize(n,maxnum=10000,df=twobody,rtol=[1.0e-12_wp],atol=[1.0e-12_wp])

!initial conditions:
x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
        1.0_wp,2.0_wp,3.0_wp]
t0 = 0.0_wp       !initial time (sec)
tf = 1000.0_wp    !final time (sec)

write(*,'(A/,*(F15.6/))') 'Initial time:',t0
write(*,'(A/,*(F15.6/))') 'Initial state:',x0
t = t0
x = x0
call s%integrate(t,x,tf,idid=idid)
write(*,'(A/,*(F15.6/))') 'Final time:',t
write(*,'(A/,*(F15.6/))') 'Initial time:',x

contains

    subroutine twobody(me,t,x,xdot)

        !! derivative routine for two-body orbit propagation

        implicit none

        class(ddeabm_class),intent(inout) :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot

        real(wp),dimension(3) :: r,v,a_grav
        real(wp) :: rmag

        r = x(1:3)
        v = x(4:6)
        rmag = norm2(r)
        a_grav = -mu/rmag**3 * r ! acceleration due to gravity

        xdot(1:3) = v
        xdot(4:6) = a_grav

    end subroutine twobody

end program ddeabm_example

```
It produces the following output:

```
Initial time:
       0.000000

Initial state:
   10000.000000
   10000.000000
   10000.000000
       1.000000
       2.000000
       3.000000

Final time:
    1000.000000

Initial time:
   10667.963305
   11658.055962
   12648.148619
       0.377639
       1.350074
       2.322509
```

## Reporting of intermediate points

The intermediate integration points can also be reported to a user-defined procedure.  For the above example, the following subroutine could be defined:

```Fortran
subroutine twobody_report(me,t,x)

    !! report function - write time,state to console

    implicit none

    class(ddeabm_class),intent(inout)    :: me
    real(wp),intent(in)                  :: t
    real(wp),dimension(:),intent(in)     :: x

    write(*,'(*(F15.6,1X))') t,x

end subroutine twobody_report
```
Which can be added to the class on initialization:

```Fortran
call s%initialize(n,maxnum=10000,df=twobody,&
                  rtol=[1.0e-12_wp],atol=[1.0e-12_wp],&
                  report=twobody_report)

```
This function is then called at each time step
if the equations are integrated using the `integration_mode=2` option like so:
```Fortran
call s%integrate(t,x,tf,idid=idid,integration_mode=2)
```

## Event location

A user-defined event function `g(t,x)` can also be defined in order to stop the integration at a specified event (i.e., when `g(t,x)=0`). In the above example, say it is desired that the integration stop when `z = x(3) = 12,000 km`.  The event function for this would be:

```Fortran
subroutine twobody_event(me,t,x,g)

    !! event function for z = 12,000 km

    implicit none

    class(ddeabm_with_event_class),intent(inout) :: me
    real(wp),intent(in)                          :: t
    real(wp),dimension(:),intent(in)             :: x
    real(wp),intent(out)                         :: g

    g = 12000.0_wp - x(3)

end subroutine twobody_event
```
For event finding, the `ddeabm_with_event_class` type is used (which is an extension of the main `ddeabm_class`).  For example:

```Fortran
type(ddeabm_with_event_class) :: s
...
call s%initialize_event(n,maxnum=10000,df=twobody,&
                        rtol=[1.0e-12_wp],atol=[1.0e-12_wp],&
                        g=twobody_event,root_tol=1.0e-12_wp)
...
call s%integrate_to_event(t,x,tf,idid=idid,gval=gval)
```
In this case, `root_tol` is the tolerance for the event location, and `gval` is the value of the event function at the final time (note that the integration will stop when `g(t,x)=0` or at `t=tf`, whichever occurs first).

# Documentation

The latest API documentation can be found [here](http://jacobwilliams.github.io/ddeabm/). This was generated from the source code using [FORD](https://github.com/cmacmackin/ford) (note that the included `build.sh` script will also generate these files).

# License

The ddeabm source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/ddeabm/blob/master/LICENSE) (BSD-style).  The original DDEABM Fortran 77 code is [public domain](http://www.netlib.org/slatec/guide).

# References

1. L. F. Shampine, M. K. Gordon, "Solving ordinary differential equations with ODE, STEP, and INTRP",  Report SLA-73-1060, Sandia Laboratories, 1973.
2. L. F. Shampine, M. K. Gordon, "Computer solution of ordinary differential equations, the initial value problem", W. H. Freeman and Company, 1975.
3. L. F. Shampine, H. A. Watts, "DEPAC - Design of a user oriented package of ode solvers", Report SAND79-2374, Sandia Laboratories, 1979.
4. H. A. Watts, "A smoother interpolant for DE/STEP, INTRP and DEABM: II", Report SAND84-0293, Sandia Laboratories, 1984.
5. R. P. Brent, "[An algorithm with guaranteed convergence for finding a zero of a function](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)", The Computer Journal, Vol 14, No. 4., 1971.
6. R. P. Brent, "[Algorithms for minimization without derivatives](http://maths-people.anu.edu.au/~brent/pub/pub011.html)", Prentice-Hall, Inc., 1973.
