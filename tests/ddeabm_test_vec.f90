!*****************************************************************************************
!> author: Jacob Williams
!  date: 12/17/2017
!
!  Unit test for [[ddeabm_class]].
!  Tests multiple event functions.

    program ddeabm_test_vec

    use ddeabm_kinds
    use ddeabm_module

    implicit none

    type,extends(ddeabm_with_event_class_vec) :: spacecraft
        !! spacecraft propagation type.
        !! extends the [[ddeabm_class]] to include data used in the deriv routine
        real(wp) :: mu = 0.0_wp      !! central body gravitational parameter (km3/s2)
        integer :: fevals = 0        !! number of function evaluations
        logical :: first = .true.    !! first point is being exported
    end type spacecraft

    integer,parameter  :: n = 6            !! number of state variables
    real(wp),parameter :: tol = 1.0e-12_wp !! event location tolerance
    integer,parameter  :: ng = 3           !! number of event functions

    type(spacecraft) :: s
    real(wp),dimension(n) :: x0,xf,x02,x
    real(wp) :: t0,tf,dt,gf,tf_actual,t
    integer :: idid
    real(wp),dimension(ng) :: gval
    real(wp),dimension(3) :: r_target

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' ddeabm_test'
    write(*,*) '---------------'
    write(*,*) ''

    !integration to event test (integrate until z-coordinate is 12,000 km)

    write(*,*) ''
    write(*,*) '-----------------------'
    write(*,*) 'integration to event:'
    write(*,*) '-----------------------'
    write(*,*) ''

    call s%initialize_event(n,maxnum=10000,df=twobody,rtol=[1.0e-12_wp],atol=[1.0e-12_wp],&
                                ng=ng,g=twobody_event,root_tol=[1.0e-12_wp],report=twobody_report)
    s%mu = 398600.436233_wp  !earth
    s%fevals = 0

    !initial conditions:
    x0 = [10000.0_wp,10000.0_wp,10000.0_wp,&   !initial state [r,v] (km,km/s)
            1.0_wp,2.0_wp,3.0_wp]
    t0 = 0.0_wp       !initial time (sec)

    write(*,'(A/,*(F15.6/))') 'Initial time:',t0
    write(*,'(A/,*(F15.6/))') 'Initial state:',x0
    s%fevals = 0
    s%first = .true.
    t = t0
    x = x0
    tf = 2000.0_wp    !max final time (sec)
    r_target = [0.0_wp, 11000.0_wp, 12000.0_wp]  ! three event functions (x,y,z)
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

    ! Now, continue integration (until next event)

    write(*,*) '======================='
    write(*,*) 'continue integration...'
    write(*,*) '======================='

    s%mu = 398600.436233_wp  !earth
    s%fevals = 0
    s%first = .true.
    x = xf
    write(*,'(A/,*(F15.6/))') 'Initial time     :',t
    write(*,'(A/,*(F15.6/))') 'Max final time   :',tf
    write(*,'(A/,*(F15.6/))') 'Initial state    :',x
    !call s%first_call()  ! instead of restarting, we use the "continue" argument, for efficiency
    ! test fixed output step here:
    call s%integrate_to_event(t,x,tf,idid=idid,gval=gval,&
                                integration_mode=2,tstep=25.0_wp,continue=.true.) 
    xf = x
    write(*,*) ''
    write(*,'(A/,*(I5/))')    'idid: ',idid
    write(*,'(A/,*(F30.16/))') 'gval: ',gval
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
        subroutine twobody_event(me,t,x,g,ig)

        !! 3 event functions for [x,y,z] = r_target

        implicit none

        class(ddeabm_with_event_class_vec),intent(inout) :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(:),intent(in)     :: x
        real(wp),dimension(:),intent(out)    :: g
        integer,intent(in),optional          :: ig  !! the event function to compute

        logical,dimension(ng) :: compute
        integer :: i !! counter

        if (.not. present(ig)) then
            compute = .true. ! compute them all
        else
            compute = .false.
            compute(ig) = .true. ! only compute this one
        end if

        !where (compute) g = r_target - x   ! gfortran 11.1 is crashing here.
        do i = 1, ng 
            if (compute(i)) g(i) = r_target(i) - x(i)
        end do

        end subroutine twobody_event
    !*********************************************************

    end program ddeabm_test_vec
!*****************************************************************************************
