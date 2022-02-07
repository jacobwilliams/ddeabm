!*****************************************************************************************
!> author: Jacob Williams
!  date: 12/16/2015
!
!  Unit test for [[ddeabm_class]].
!  Integrate a two-body orbit around the Earth.
!  Also tests integration to an event.

    program ddeabm_test

    use ddeabm_module, wp => ddeabm_rk

    implicit none

    integer,parameter  :: n = 6            !! number of state variables
    real(wp),parameter :: tol = 1.0e-12_wp !! event location tolerance
    real(wp),parameter :: rtol = 1.0e-12_wp !! integrator tolerances
    real(wp),parameter :: atol = 1.0e-12_wp !! integrator tolerances
    integer,parameter :: maxnum = 10000  !! max func evals
    real(wp),parameter :: mu = 398600.436233_wp !! MU for Earth

    real(wp),dimension(6),parameter :: x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !! initial state [r,v] (km,km/s)
                                             1.0_wp,2.0_wp,3.0_wp]
    real(wp),parameter :: t0 = 0.0_wp       !! initial time (sec)


    type,extends(ddeabm_with_event_class) :: spacecraft
        !! spacecraft propagation type.
        !! extends the [[ddeabm_class]] to include data used in the deriv routine
        real(wp) :: mu = mu          !! central body gravitational parameter (km3/s2)
        integer :: fevals = 0        !! number of function evaluations
        logical :: first = .true.    !! first point is being exported
    end type spacecraft

    type(spacecraft) :: s, s2
    real(wp),dimension(n) :: xf,x02,x
    real(wp) :: tf,dt,gf,tf_actual,t,gval,z_target
    integer :: idid
    integer :: icase

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' ddeabm_test'
    write(*,*) '---------------'
    write(*,*) ''

    ! !***************************************************************************

    write(*,*) ''
    write(*,*) '-----------------------'
    write(*,*) 'forward/backward test:'
    write(*,*) '-----------------------'
    write(*,*) ''

    !constructor (main body is Earth):

    call s%initialize(n,maxnum=maxnum,df=twobody,rtol=[rtol],atol=[atol],&
                        report=twobody_report)

    !initial conditions:
    tf = 1000.0_wp    !final time (sec)
    s%fevals = 0

    write(*,'(A/,*(F15.6/))') 'Initial time:',t0
    write(*,'(A/,*(F15.6/))') 'Initial state:',x0
    s%fevals = 0
    s%first = .true.
    t = t0
    x = x0
    call s%first_call()
    call s%integrate(t,x,tf,idid=idid,integration_mode=2)    !forward (report points)
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
    call s%integrate(t,x,t0,idid=idid)  !backwards
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

    ! ! !***************************************************************************

    !integration to event test (integrate until z-coordinate is 12,000 km)

    write(*,*) ''
    write(*,*) '-----------------------'
    write(*,*) 'integration to event [default steps]:'
    write(*,*) '-----------------------'
    write(*,*) ''

    call s%destroy()
    call s%initialize_event(n,maxnum=maxnum,df=twobody,rtol=[rtol],atol=[atol],&
                                g=twobody_event,root_tol=tol,report=twobody_report)
    write(*,'(A/,*(F15.6/))') 'Initial time:',t0
    write(*,'(A/,*(F15.6/))') 'Initial state:',x0
    s%fevals = 0
    s%first = .true.
    t = t0
    x = x0
    tf = 1000.0_wp    !max final time (sec)
    z_target = 12000.0_wp
    call s%first_call()
    call s%integrate_to_event(t,x,tf,idid=idid,gval=gval,integration_mode=2)
    xf = x
    write(*,*) ''
    write(*,'(A/,*(I5/))')    'idid: ',idid
    write(*,'(A/,*(F30.16/))') 'gval: ',gval
    write(*,'(A/,*(F15.6/))') 'Final time:',t
    write(*,'(A/,*(F15.6/))') 'Final state:',xf
    write(*,'(A,I5)') 'Function evaluations:', s%fevals
    write(*,*) ''

    ! write(*,*) ''
    ! write(*,*) '-----------------------'
    ! write(*,*) 'integration to event [fixed steps]:'
    ! write(*,*) '-----------------------'
    ! write(*,*) ''

    ! call s%destroy()
    ! call s%initialize_event(n,maxnum=maxnum,df=twobody,rtol=[rtol],atol=[atol],&
    !                             g=twobody_event,root_tol=tol,report=twobody_report)
    ! write(*,'(A/,*(F15.6/))') 'Initial time:',t0
    ! write(*,'(A/,*(F15.6/))') 'Initial state:',x0
    ! s%fevals = 0
    ! s%first = .true.
    ! t = t0
    ! x = x0
    ! tf = 1000.0_wp    !max final time (sec)
    ! z_target = 12000.0_wp
    ! call s%first_call()
    ! call s%integrate_to_event(t,x,tf,idid=idid,gval=gval,integration_mode=2,tstep=50.0_wp)
    ! xf = x
    ! write(*,*) ''
    ! write(*,'(A/,*(I5/))')    'idid: ',idid
    ! write(*,'(A/,*(F30.16/))') 'gval: ',gval
    ! write(*,'(A/,*(F15.6/))') 'Final time:',t
    ! write(*,'(A/,*(F15.6/))') 'Final state:',xf
    ! write(*,'(A,I5)') 'Function evaluations:', s%fevals
    ! write(*,*) ''

    ! write(*,*) '======================='
    ! write(*,*) 'integrate normally to root and compare ...'
    ! write(*,*) '======================='
    ! call s%first_call()  !have to restart the integration after a root finding
    ! tf = t
    ! t = t0
    ! x = x0
    ! call s%integrate(t,x,tf,idid=idid,integration_mode=2)
    ! call twobody_event(s,t,x,gval)
    ! write(*,*) ''
    ! write(*,'(A/,*(I5/))')    'idid: ',idid
    ! write(*,'(A/,*(F30.16/))') 'gval: ',gval
    ! write(*,'(A/,*(F15.6/))') 'Final time:',t
    ! write(*,'(A/,*(F15.6/))') 'Final state:',x
    ! write(*,*) ''

    ! write(*,*) '======================='
    ! write(*,*) 'integrate with fixed step size to root and compare ...'
    ! write(*,*) '======================='
    ! call s%first_call()  !have to restart the integration after a root finding
    ! tf = t
    ! t = t0
    ! x = x0
    ! call s%integrate(t,x,tf,idid=idid,integration_mode=2,tstep=50.0_wp)
    ! call twobody_event(s,t,x,gval)
    ! write(*,*) ''
    ! write(*,'(A/,*(I5/))')    'idid: ',idid
    ! write(*,'(A/,*(F30.16/))') 'gval: ',gval
    ! write(*,'(A/,*(F15.6/))') 'Final time:',t
    ! write(*,'(A/,*(F15.6/))') 'Final state:',x
    ! write(*,*) ''

    ! stop
    ! ! ....................

    ! !--------------------------------------------------------------------------------------

    ! Now, continue integration (until z-coordinate is 13,000 km)

    write(*,*) '======================='
    write(*,*) 'continue integration...'
    write(*,*) '======================='

    s%fevals = 0
    s%first = .true.
    x = xf
    tf = 2000.0_wp    !max final time (sec)
    z_target = 13000.0_wp !change the target value

    !...new ... try to continue without restarting ...
    call s%integrate_to_event(t,x,tf,idid=idid,gval=gval,integration_mode=2,tstep=50.0_wp,continue=.true.) ! test dense output here
    !call s%integrate_to_event(t,x,tf,idid=idid,gval=gval,integration_mode=2,continue=.true.) ! default points

    !...original... restart integration ...
    ! call s%first_call()  !have to restart the integration after a root finding
    ! call s%integrate_to_event(t,x,tf,idid=idid,gval=gval,integration_mode=2,tstep=50.0_wp) ! test dense output here
    ! !call s%integrate_to_event(t,x,tf,idid=idid,gval=gval,integration_mode=2) ! default points

    xf = x
    write(*,*) ''
    write(*,'(A/,*(I5/))')    'idid: ',idid
    write(*,'(A/,*(F30.16/))') 'gval: ',gval
    write(*,'(A/,*(F15.6/))') 'Final time:',t
    write(*,'(A/,*(F15.6/))') 'Final state:',xf
    write(*,'(A,I5)') 'Function evaluations:', s%fevals
    write(*,*) ''

    write(*,*) '======================='
    write(*,*) 'integrate normally to root and compare ...'
    write(*,*) '======================='

    call s2%initialize(n,maxnum=maxnum,df=twobody,rtol=[rtol],atol=[atol])
    call s2%first_call()  !have to restart the integration after a root finding
    s2%fevals = 0
    tf = t
    t = t0
    x = x0
    call s2%integrate(t,x,tf,idid=idid,integration_mode=1)
    call twobody_event(s2,t,x,gval)
    write(*,*) ''
    write(*,'(A/,*(I5/))')     'idid: ',idid
    write(*,'(A/,*(F30.16/))') 'gval: ',gval
    write(*,'(A/,*(F15.6/))')  'Initial time:',t0
    write(*,'(A/,*(F15.6/))')  'Final time:',t
    write(*,'(A/,*(F15.6/))')  'Final state:',x
    write(*,'(A,I5)') 'Function evaluations:', s2%fevals
    write(*,*) ''

    ! write(*,*) '======================='
    ! write(*,*) 'integrate with fixed step size to root and compare ...'
    ! write(*,*) '======================='
    ! call s2%initialize(n,maxnum=maxnum,df=twobody,rtol=[rtol],atol=[atol])
    ! call s2%first_call()  !have to restart the integration after a root finding
    ! s2%fevals = 0
    ! tf = t
    ! t = t0
    ! x = x0
    ! call s2%integrate(t,x,tf,idid=idid,integration_mode=1,tstep=50.0_wp)
    ! call twobody_event(s2,t,x,gval)
    ! write(*,*) ''
    ! write(*,'(A/,*(I5/))')     'idid: ',idid
    ! write(*,'(A/,*(F30.16/))') 'gval: ',gval
    ! write(*,'(A/,*(F15.6/))')  'Initial time:',t0
    ! write(*,'(A/,*(F15.6/))')  'Final time:',t
    ! write(*,'(A/,*(F15.6/))')  'Final state:',x
    ! write(*,'(A,I5)') 'Function evaluations:', s2%fevals
    ! write(*,*) ''

    ! !****************************************************************************************

    ! !
    ! ! try continuing without root finding. This accurate for both default and fixed steps
    ! !

    ! do icase = 1, 2  !first with default steps, next with fixed steps

    !     write(*,*) ''
    !     write(*,*) '========================================================='
    !     if (icase==1) then
    !         write(*,*) ' default steps'
    !     else
    !         write(*,*) ' fixed 50 sec step'
    !     end if
    !     write(*,*) '========================================================='
    !     write(*,*) ''

    !     write(*,*) ''
    !     write(*,*) '====================================='
    !     write(*,*) ' integrate from t=0 to t=300'

    !     call s%initialize(n,maxnum=maxnum,df=twobody,rtol=[rtol],atol=[atol],&
    !                         report=twobody_report)

    !     !initial conditions:
    !     tf = 300.0_wp    !final time (sec)
    !     s%fevals = 0
    !     s%fevals = 0
    !     s%first = .true.
    !     t = t0
    !     x = x0
    !     call s%first_call()
    !     if (icase==1) then
    !         call s%integrate(t,x,tf,idid=idid,integration_mode=1)
    !     else
    !         call s%integrate(t,x,tf,idid=idid,integration_mode=1,tstep=50.0_wp)
    !     end if
    !     xf = x
    !     write(*,*) ''
    !     write(*,'(A/,*(I5/))')     'idid: ',idid
    !     write(*,'(A/,*(F15.6/))')  'Initial time:',t0
    !     write(*,'(A/,*(F15.6/))')  'Final time:',t
    !     write(*,'(A/,*(F15.6/))')  'Final state:',x
    !     write(*,'(A,I5)') 'Function evaluations:', s%fevals
    !     write(*,*) ''

    !     write(*,*) ' continue to t=600'
    !     t = tf
    !     tf = 600.0_wp
    !     if (icase==1) then
    !         call s%integrate(t,x,tf,idid=idid,integration_mode=1)
    !     else
    !         call s%integrate(t,x,tf,idid=idid,integration_mode=1,tstep=50.0_wp)
    !     end if
    !     write(*,*) ''
    !     write(*,'(A/,*(I5/))')     'idid: ',idid
    !     write(*,'(A/,*(F15.6/))')  'Final time:',t
    !     write(*,'(A/,*(F15.6/))')  'Final state:',x
    !     write(*,'(A,I5)') 'Function evaluations:', s%fevals
    !     write(*,*) ''

    !     write(*,*) ''
    !     write(*,*) '====================================='
    !     write(*,*) ' integrate all at once from t=0 to t=600'

    !     call s%initialize(n,maxnum=maxnum,df=twobody,rtol=[rtol],atol=[atol],&
    !                         report=twobody_report)
    !     !initial conditions:
    !     tf = 600.0_wp    !final time (sec)
    !     s%fevals = 0
    !     s%fevals = 0
    !     s%first = .true.
    !     t = t0
    !     x = x0
    !     call s%first_call()
    !     if (icase==1) then
    !         call s%integrate(t,x,tf,idid=idid,integration_mode=1)
    !     else
    !         call s%integrate(t,x,tf,idid=idid,integration_mode=1,tstep=50.0_wp)
    !     end if
    !     xf = x
    !     write(*,*) ''
    !     write(*,'(A/,*(I5/))')     'idid: ',idid
    !     write(*,'(A/,*(F15.6/))')  'Initial time:',t0
    !     write(*,'(A/,*(F15.6/))')  'Final time:',t
    !     write(*,'(A/,*(F15.6/))')  'Final state:',x
    !     write(*,'(A,I5)') 'Function evaluations:', s%fevals
    !     write(*,*) ''

    !     write(*,*) '====================================='

    ! end do

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

    end program ddeabm_test
!*****************************************************************************************
