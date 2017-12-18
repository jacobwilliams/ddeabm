!*****************************************************************************************
!>
!  Test program for [[ddeabm_class]].
!  Generates a speed/accuracy plot.
!
!@note Requires [pyplot-fortran](https://github.com/jacobwilliams/pyplot-fortran).
!
!### Results
!  ![Plot of results](|media|/ddeabm_performance_test.png)

    program ddeabm_performance_test

    use ddeabm_module
    use ddeabm_kinds
    use pyplot_module
    use iso_fortran_env

    integer,parameter :: n = 6  !! number of state variables

    type,extends(ddeabm_class) :: spacecraft
        !! spacecraft propagation type.
        !! extends the [[ddeabm_class]] to include data used in the deriv routine
        real(wp) :: mu     = 0.0_wp  !! central body gravitational parameter ( \( km^3/s^2 \) )
        integer  :: fevals = 0       !! number of function evaluations
        logical  :: first  = .true.  !! first point is being exported
    end type spacecraft

    integer,parameter :: exp_min = 8   !! min exponent for tolerance
    integer,parameter :: exp_max = 13  !! max exponent for tolerance

    real(wp),dimension(n) :: x0, acc
    real(wp) :: tol
    integer :: i,j,num_func_evals
    type(pyplot) :: plt
    real(wp),dimension(exp_min:exp_max) :: err_vec_x,err_vec_y,err_vec_z
    real(wp),dimension(exp_min:exp_max) :: err_vec_vx,err_vec_vy,err_vec_vz
    real(wp),dimension(exp_min:exp_max) :: fevals_vec
    character(len=10) :: istr
    character(len=:),allocatable :: kind_str
    integer :: istat !! status code

    character(len=1),dimension(7) :: colors = ['b','g','r','c','m','k','y'] !! plot colors

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' ddeabm_performance_test'
    write(*,*) '---------------'
    write(*,*) ''

    select case(wp)
    case(real32)
        kind_str = '(real32)'
    case(real64)
        kind_str = '(real64)'
    case(real128)
        kind_str = '(real128)'
    case default
        error stop 'error: unknown real kind'
    end select

    !initialize the plot:
    call plt%initialize(grid=.true.,xlabel='Number of Digits of Accuracy',&
                        ylabel='Number of Function Evaluations',&
                        title='DDEABM Performance '//kind_str,legend=.true.)

    !5 cases with different initial states
    do j=1,5

        write(istr,'(I5)') j  !orbit case string

        !initial state [r,v] (km,km/s)
        x0 = [10000.0_wp*j,10000.0_wp*j,10000.0_wp*j,&
                1.0_wp,2.0_wp,3.0_wp]

        do i=exp_min,exp_max
            tol = 10.0_wp**(-i)
            call go(x0,tol,tol,num_func_evals,acc)
            err_vec_x(i)  = acc(1)
            err_vec_y(i)  = acc(2)
            err_vec_z(i)  = acc(3)
            err_vec_vx(i) = acc(4)
            err_vec_vy(i) = acc(5)
            err_vec_vz(i) = acc(6)
            fevals_vec(i) = num_func_evals
        end do

        !generate the plot:
        call plt%add_plot(dble(err_vec_x),dble(fevals_vec),label='Test '//trim(adjustl(istr))//' r',&
                            linestyle=colors(1+mod(j,size(colors)))//'o-',&
                            markersize=5,linewidth=2,istat=istat)
        call plt%add_plot(dble(err_vec_y),dble(fevals_vec),label='',&
                            linestyle=colors(1+mod(j,size(colors)))//'o-',&
                            markersize=5,linewidth=2,istat=istat)
        call plt%add_plot(dble(err_vec_z),dble(fevals_vec),label='',&
                            linestyle=colors(1+mod(j,size(colors)))//'o-',&
                            markersize=5,linewidth=2,istat=istat)
        call plt%add_plot(dble(err_vec_vx),dble(fevals_vec),label='Test '//trim(adjustl(istr))//' v',&
                            linestyle=colors(1+mod(j,size(colors)))//'o--',&
                            markersize=5,linewidth=2,istat=istat)
        call plt%add_plot(dble(err_vec_vy),dble(fevals_vec),label='',&
                            linestyle=colors(1+mod(j,size(colors)))//'o--',&
                            markersize=5,linewidth=2,istat=istat)
        call plt%add_plot(dble(err_vec_vz),dble(fevals_vec),label='',&
                            linestyle=colors(1+mod(j,size(colors)))//'o--',&
                            markersize=5,linewidth=2,istat=istat)

    end do

    call plt%savefig('ddeabm_performance_test.png',istat=istat)

    contains
!*****************************************************************************************

!***************************************************************************
    subroutine go(x0,rtol,atol,fevals,num_digits)

        !! run the test.

        implicit none

        real(wp),dimension(n),intent(in)  :: x0          !! initial state [r,v]
        real(wp),intent(in)               :: rtol        !! absolute tolerance
        real(wp),intent(in)               :: atol        !! relative tolerance
        integer,intent(out)               :: fevals      !! number of function evaluations for forward propagation
        real(wp),dimension(n),intent(out) :: num_digits  !! number of digits of accuracy [r,v]

        type(spacecraft) :: s
        real(wp),dimension(n) :: xf,x02,x,errvec
        real(wp) :: t0,tf,dt,gf,tf_actual,t,rerr,verr
        integer :: idid

        !constructor (main body is Earth):
        call s%initialize(n,maxnum=10000,df=twobody,rtol=[rtol],atol=[atol])
        s%mu = 398600.436233_wp  !earth
        t0   = 0.0_wp            !initial time (sec)
        tf   = 1000.0_wp         !final time (sec)

        s%fevals = 0
        s%first = .true.
        t = t0
        x = x0
        call s%first_call()
        call s%integrate(t,x,tf,idid=idid)    !forward
        xf = x
        fevals = s%fevals  !number of function evaluations

        t = tf
        x = xf
        s%fevals = 0
        call s%first_call()  !restarting the integration
        call s%integrate(t,x,t0,idid=idid)  !backwards
        x02 = x

        !number of digits of accuracy
        errvec = abs(x02 - x0) / abs(x0) !relative error
        where (errvec==0.0_wp)
            errvec = epsilon(1.0_wp) !small number
        end where
        num_digits = abs(log10(errvec))

        !write(*,*) x02(1), x0(1), x02(1)-x0(1), errvec(1), num_digits(1)

    end subroutine go
    !*********************************************************

    !*********************************************************
    subroutine twobody(me,t,x,xdot)

        !! derivative routine for two-body orbit propagation.

        implicit none

        class(ddeabm_class),intent(inout) :: me
        real(wp),intent(in)               :: t    !! time [\(s\)]
        real(wp),dimension(:),intent(in)  :: x    !! state vector [\(km\), \(km/s\)]
        real(wp),dimension(:),intent(out) :: xdot !! state vector derivative [\(km/s\), \(km/s^2\)]

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

    end program ddeabm_performance_test
!*****************************************************************************************
