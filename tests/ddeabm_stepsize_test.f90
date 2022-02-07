!*****************************************************************************************
!> author: Jacob Williams
!  date: 1/1/2018
!
!  Unit test for [[ddeabm_class]].
!  Tests the different options for specifying the initial step size.

    program ddeabm_stepsize_test

    use ddeabm_module, wp => ddeabm_rk

    implicit none

    type,extends(ddeabm_with_event_class) :: spacecraft
        !! spacecraft propagation type.
        !! extends the [[ddeabm_class]] to include data used in the deriv routine
        real(wp) :: mu = 0.0_wp      !! central body gravitational parameter (km3/s2)
        integer :: fevals = 0        !! number of function evaluations
        logical :: first = .true.    !! first point is being exported
    end type spacecraft

    integer,parameter  :: n = 6            !! number of state variables
    real(wp),parameter :: tol = 1.0e-12_wp !! event location tolerance
    integer,parameter :: integration_mode = 1  !! don't report points
    real(wp),parameter :: initial_step_size = 1.56e-4_wp !! for `initial_step_mode` 3

    type(spacecraft) :: s
    real(wp),dimension(n) :: x0,xf,x02,x
    real(wp) :: t0,tf,dt,gf,tf_actual,t,gval,z_target
    integer :: idid,initial_step_mode

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' ddeabm_test'
    write(*,*) '---------------'
    write(*,*) ''

    !***************************************************************************

    do initial_step_mode = 1, 3  ! test each mode

        !constructor (main body is Earth):
        call s%initialize(n,maxnum=10000,df=twobody,&
                            rtol=[1.0e-12_wp],atol=[1.0e-12_wp],&
                            report=twobody_report,&
                            initial_step_mode=initial_step_mode,&
                            initial_step_size=initial_step_size)
        s%mu = 398600.436233_wp  !earth

        !initial conditions:
        x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
                1.0_wp,2.0_wp,3.0_wp]
        t0 = 0.0_wp       !initial time (sec)
        tf = 1000.0_wp    !final time (sec)
        s%fevals = 0

        if (initial_step_mode==1) then
            write(*,'(A/,*(F15.6/))') 'Initial time:',t0
            write(*,'(A/,*(F15.6/))') 'Initial state:',x0
        end if
        s%fevals = 0
        s%first = .true.
        t = t0
        x = x0
        call s%first_call()  ! restart integration
        call s%integrate(t,x,tf,idid=idid,integration_mode=integration_mode)
        xf = x
        write(*,*) ''
        write(*,*) '------------------------------------------------------'
        write(*,'(A/,*(I5/))')    'initial_step_mode: ',initial_step_mode
        write(*,'(A/,*(I5/))')    'idid: ',idid
        write(*,'(A/,*(F15.6/))') 'Final time:',t
        write(*,'(A/,*(F15.6/))') 'Final state:',xf
        write(*,'(A,I5)') 'Function evaluations:', s%fevals
        write(*,*) ''

    end do

    contains
!*****************************************************************************************

    !*********************************************************
        subroutine twobody(me,t,x,xdot)

        !! derivative routine for two-body orbit propagation

        implicit none

        class(ddeabm_class),intent(inout) :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot

        real(wp),dimension(3) :: r,v,a_grav
        real(wp) :: rmag

        select type (me)
        class is (spacecraft)

            r = x(1:3)
            v = x(4:6)
            rmag = norm2(r)
            a_grav = -me%mu/rmag**3 * r !acceleration due to gravity

            xdot(1:3) = v
            xdot(4:6) = a_grav

            me%fevals = me%fevals + 1

        end select

        end subroutine twobody
    !*********************************************************

    !*********************************************************
        subroutine twobody_report(me,t,x)

        !! report function - write time,state to console

        implicit none

        class(ddeabm_class),intent(inout)    :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(:),intent(in)     :: x

        select type (me)
        class is (spacecraft)
            if (me%first) then  !print header
                write(*,*) ''
                write(*,'(*(A15,1X))')  'time (sec)','x (km)','y (km)','z (km)',&
                                        'vx (km/s)','vy (km/s)','vz (km/s)'
                me%first = .false.
            end if
        end select

        write(*,'(*(F15.6,1X))') t,x

        end subroutine twobody_report
    !*********************************************************

    end program ddeabm_stepsize_test
!*****************************************************************************************
