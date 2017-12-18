!*****************************************************************************************
!> author: Jacob Williams
!  date: 12/22/2015
!
!  Unit tests for [[zeroin]] subroutine.

    program zeroin_test

    use ddeabm_kinds
    use root_module, only: zeroin

    implicit none

    real(wp) :: ax,bx,xzero,fzero
    integer :: iflag,fevals

    real(wp),parameter :: tol     = 1.0e-12_wp  !! tolerance for root location
    real(wp),parameter :: rad2deg = 180.0_wp / acos(-1.0_wp) !! radians to degrees

    write(*,*) ''
    write(*,*) '---------------'
    write(*,*) ' zeroin_test'
    write(*,*) '---------------'
    write(*,*) ''

    fevals = 0
    ax = 0.0_wp
    bx = 4.0_wp
    call zeroin(sin_func,ax,bx,tol,xzero,fzero,iflag)
    write(*,*) 'sin:', iflag, xzero*rad2deg, fzero, fevals

    fevals = 0
    ax = -1.0_wp
    bx = 4.0_wp
    call zeroin(cos_func,ax,bx,tol,xzero,fzero,iflag)
    write(*,*) 'cos:', iflag, xzero*rad2deg, fzero, fevals

    fevals = 0
    ax = 2.0_wp
    bx = 4.5_wp
    call zeroin(tan_func,ax,bx,tol,xzero,fzero,iflag)
    write(*,*) 'tan:', iflag, xzero*rad2deg, fzero, fevals

    fevals = 0
    ax = 0.5_wp
    bx = 2.0_wp
    call zeroin(log_func,ax,bx,tol,xzero,fzero,iflag)
    write(*,*) 'log:', iflag, xzero, fzero, fevals

    contains

        function sin_func(x) result(f)  !! \( \sin(x) \)
        implicit none
        real(wp),intent(in)  :: x
        real(wp)             :: f
        f = sin(x)
        fevals = fevals + 1
        end function sin_func

        function cos_func(x) result(f)  !! \( \cos(x) \)
        implicit none
        real(wp),intent(in)  :: x
        real(wp)             :: f
        f = cos(x)
        fevals = fevals + 1
        end function cos_func

        function tan_func(x) result(f)  !! \( \tan(x) \)
        implicit none
        real(wp),intent(in)  :: x
        real(wp)             :: f
        f = tan(x)
        fevals = fevals + 1
        end function tan_func

        function log_func(x) result(f)  !! \( \log(x) \)
        implicit none
        real(wp),intent(in)  :: x
        real(wp)             :: f
        f = log(x)
        fevals = fevals + 1
        end function log_func

    end program zeroin_test
!*****************************************************************************************
