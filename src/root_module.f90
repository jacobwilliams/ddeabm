!*****************************************************************************************
!> author: Jacob Williams
!  date: 12/22/2015
!
!  Root-finding routines. Used by [[ddeabm_module]].

    module root_module

    use kind_module

    implicit none

    private

    !real parameters:
    real(wp),parameter :: zero  = 0.0_wp
    real(wp),parameter :: one   = 1.0_wp
    real(wp),parameter :: two   = 2.0_wp
    real(wp),parameter :: three = 3.0_wp
    real(wp),parameter :: eps   = epsilon(one)    !! original code had d1mach(4)

    abstract interface
        function func(t) result(f)  !! interface to the [[zeroin]] function
        import :: wp
        implicit none
        real(wp),intent(in)  :: t
        real(wp)             :: f
        end function func
    end interface

    public :: zeroin

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Find a zero of the function `f(x)` in the interval `ax`,`bx`.
!
!  It is assumed that `f(ax)` and `f(bx) `have opposite signs
!  this is checked, and an error message is printed if this is not
!  satisfied. zeroin returns a zero `x` in the given interval
!  `ax`,`bx` to within a tolerance `4*macheps*abs(x)+tol`, where `macheps` is
!  the relative machine precision defined as the smallest representable
!  number such that `1.0_wp+macheps>1.0_wp`.
!
!#References
!  * R. P. Brent, "An algorithm with guaranteed convergence for
!    finding a zero of a function", The Computer Journal,
!    Vol 14, No. 4. (1971). [[link](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)]
!  * R. P. Brent, "Algorithms for minimization without derivatives",
!    prentice-hall, inc. (1973).
!
!# See also
!  [1] [zeroin.f](http://www.netlib.org/go/zeroin.f)

    subroutine zeroin(f,ax,bx,tol,xzero,fzero,iflag,fax,fbx)

    use iso_fortran_env, only: error_unit

    implicit none

    procedure(func)              :: f     !! the function
    real(wp),intent(in)          :: ax    !! left endpoint of initial interval
    real(wp),intent(in)          :: bx    !! right endpoint of initial interval
    real(wp),intent(in)          :: tol   !! desired length of the interval of uncertainty of the final result (>=0)
    real(wp),intent(out)         :: xzero !! abscissa approximating a zero of f in the interval ax,bx
    real(wp),intent(out)         :: fzero !! f(xzero)
    integer,intent(out)          :: iflag !! status flag (-1=error, 0=root found)
    real(wp),intent(in),optional :: fax   !! if f(ax) is already known, it can be input here
    real(wp),intent(in),optional :: fbx   !! if f(ax) is already known, it can be input here

    real(wp) :: a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s

    tol1 = eps+one

    a=ax
    b=bx

    if (present(fax)) then
        fa = fax
    else
        fa=f(a)
    end if
    if (present(fbx)) then
        fb = fbx
    else
        fb=f(b)
    end if

    !check that f(ax) and f(bx) have different signs
    if (fa == zero .or. fb == zero .or. fa * (fb/abs(fb)) <= zero) then

    iflag = 0

20  c=a
    fc=fa
    d=b-a
    e=d

30  if (abs(fc)<abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
    end if

40  tol1=two*eps*abs(b)+0.5_wp*tol
    xm = 0.5_wp*(c-b)
    if ((abs(xm)<=tol1).or.(fb==zero)) go to 150

    ! see if a bisection is forced

    if ((abs(e)>=tol1).and.(abs(fa)>abs(fb))) go to 50
    d=xm
    e=d
    go to 110
50  s=fb/fa
    if (a/=c) go to 60

    ! linear interpolation

    p=two*xm*s
    q=one-s
    go to 70

    ! inverse quadratic interpolation

60  q=fa/fc
    r=fb/fc
    p=s*(two*xm*q*(q-r)-(b-a)*(r-one))
    q=(q-one)*(r-one)*(s-one)
70  if (p<=zero) go to 80
    q=-q
    go to 90

80  p=-p
90  s=e
    e=d
    if (((two*p)>=(three*xm*q-abs(tol1*q))).or.(p>=&
        abs(0.5_wp*s*q))) go to 100
    d=p/q
    go to 110

100 d=xm
    e=d
110 a=b
    fa=fb
    if (abs(d)<=tol1) go to 120
    b=b+d
    go to 140

120 if (xm<=zero) go to 130
    b=b+tol1
    go to 140

130 b=b-tol1
140 fb=f(b)
    if ((fb*(fc/abs(fc)))>zero) go to 20
    go to 30

150 xzero = b
    fzero = fb

    else
        iflag = -1
        write(error_unit,'(A)')&
            'Error in zeroin: f(ax) and f(bx) do not have different signs.'
    end if

    end subroutine zeroin
!*****************************************************************************************

!*****************************************************************************************
    end module root_module
!*****************************************************************************************
