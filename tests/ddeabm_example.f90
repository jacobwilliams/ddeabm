program ddeabm_example

use ddeabm_module
use ddeabm_kinds

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

write(*,'(A/,*(F15.6/))') 'Initial time: ',t0
write(*,'(A/,*(F15.6/))') 'Initial state:',x0
t = t0
x = x0
call s%integrate(t,x,tf,idid=idid)
write(*,'(A/,*(F15.6/))') 'Final time: ',t
write(*,'(A/,*(F15.6/))') 'Final state:',x

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
