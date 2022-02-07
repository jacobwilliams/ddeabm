!*****************************************************************************************
!> author: Jacob Williams
!  date: 4/7/2018
!
!  Unit test for [[ddeabm_class]].
!  Test fixed-step size reporting.

    program ddeabm_fixed_step_test

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
    real(wp),parameter :: dt = 100.0_wp    !! output step size (sec)

    type(spacecraft) :: s
    real(wp),dimension(n) :: x0,xf,x02,x
    real(wp) :: t0,tf,gf,tf_actual,t,gval,z_target
    integer :: idid

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' ddeabm_fixed_step_test'
    write(*,*) '---------------'
    write(*,*) ''

    !***************************************************************************

    write(*,*) ''
    write(*,*) '-----------------------'
    write(*,*) 'forward/backward test:'
    write(*,*) '-----------------------'
    write(*,*) ''

    !constructor (main body is Earth):

    call s%initialize(n,maxnum=10000,df=twobody,rtol=[1.0e-12_wp],atol=[1.0e-12_wp],&
                        report=twobody_report)
    s%mu = 398600.436233_wp  !earth

    !initial conditions:
    x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
            1.0_wp,2.0_wp,3.0_wp]
    t0 = 0.0_wp       !initial time (sec)
    tf = 1002.1_wp    !final time (sec)
    s%fevals = 0

    write(*,'(A/,*(F15.6/))') 'Initial time:',t0
    write(*,'(A/,*(F15.6/))') 'Initial state:',x0
    s%fevals = 0
    s%first = .true.
    t = t0
    x = x0
    call s%first_call()
    call s%integrate(t,x,tf,idid=idid,integration_mode=2,tstep=dt)    !forward (report points)
    xf = x
    write(*,*) ''
    write(*,'(A/,*(I5/))')    'idid: ',idid
    write(*,'(A/,*(F15.6/))') 'Final time:',t
    write(*,'(A/,*(F15.6/))') 'Final state:',xf
    write(*,'(A,I5)') 'Function evaluations:', s%fevals
    write(*,*) ''

    t = tf
    x = xf
    s%fevals = 0
    call s%first_call()  !restarting the integration
    call s%integrate(t,x,t0,idid=idid,integration_mode=2,tstep=dt)  !backwards
    x02 = x

    write(*,'(A/,*(I5/))')    'idid: ',idid
    write(*,'(A/,*(F15.6/))') 'Initial state:',x02
    write(*,*) ''
    write(*,'(A/,*(E20.12/))') 'Error:',x02-x0
    write(*,'(A,I5)') 'Function evaluations:', s%fevals
    write(*,*) ''

    write(*,*) ''
    write(*,*) '-----------------------'
    write(*,*) 'continue integration:'
    write(*,*) '-----------------------'
    write(*,*) ''

    t = t0
    x = x02
    tf = t0 - 100.0_wp
    s%fevals = 0
    call s%first_call()  !restarting the integration (go from t -> tf)
    call s%integrate(t,x,tf,idid=idid)  !backwards
    xf = x
    write(*,'(A/,*(I5/))')    'idid: ',idid
    write(*,'(A/,*(F15.6/))') 'Final time:',t
    write(*,'(A/,*(F15.6/))') 'Final state:',xf
    write(*,'(A,I5)') 'Function evaluations:', s%fevals
    write(*,*) ''

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

    !*********************************************************
        subroutine twobody_event(me,t,x,g)

        !! event function for z = 12,000 km

        implicit none

        class(ddeabm_with_event_class),intent(inout) :: me
        real(wp),intent(in)                          :: t
        real(wp),dimension(:),intent(in)             :: x
        real(wp),intent(out)                         :: g

        g = z_target - x(3)  !z coordinate

        end subroutine twobody_event
    !*********************************************************

    end program ddeabm_fixed_step_test
!*****************************************************************************************
