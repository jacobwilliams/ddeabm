!*****************************************************************************************
!>
!  Unit tests for [[ddeabm_module]].
!
!  Checks each part of the algorithm against known solutions:
!  the integrator (accuracy, convergence, direction, tolerances,
!  initial step modes, order changes, stiffness detection, step failures),
!  the starting step algorithm, dense output (interpolation),
!  the reporting modes, the event finding wrappers (bracketing functions,
!  tolerances, continuation), user-triggered errors, and input validation.
!
!  The program stops with a nonzero exit code if any check fails.
!
!### See also
!  * `ddeabm_regression_test.f90` -- tests for specific bugs that were fixed.

module ddeabm_unit_test_module

   use ddeabm_module, wp => ddeabm_rk
   use, intrinsic :: iso_fortran_env, only: real64

   implicit none

   real(wp), parameter :: pi = acos(-1.0_wp)

   ! problems for the derivative routine:
   integer, parameter :: p_decay = 1     !! `y' = -y`
   integer, parameter :: p_osc = 2       !! `x'' = -x`
   integer, parameter :: p_cos = 3       !! `y' = cos(t)`
   integer, parameter :: p_zero = 4      !! `y' = 0`
   integer, parameter :: p_const = 5     !! `y' = 1`
   integer, parameter :: p_t = 6         !! `y' = t`
   integer, parameter :: p_blowup = 7    !! `y' = y**2`
   integer, parameter :: p_stiff = 8     !! `y' = -1000*(y - cos(t))`
   integer, parameter :: p_vdp = 9       !! Van der Pol oscillator, `mu = 5`
   integer, parameter :: p_stop = 10     !! oscillator that calls `stop_integration` for `t > t_stop`

   integer :: iprob = p_osc                !! problem to solve
   real(wp) :: t_stop = huge(1.0_wp)       !! for `p_stop`
   integer :: n_feval = 0                  !! number of derivative evaluations
   integer :: n_reports = 0                !! number of reported points
   real(wp), dimension(100000) :: t_reports = 0.0_wp  !! the reported times
   real(wp), dimension(100000) :: x_reports = 0.0_wp  !! the first state component at the reported times
   integer :: n_failed = 0                 !! number of failed checks
   integer :: n_checks = 0                 !! number of checks

   ! tolerances in the tests are specified for real64. for a lower precision
   ! real kind (REAL32), they are scaled to the same fraction of the available
   ! digits: x -> x**pw, where pw = log(eps)/log(eps_real64).
   real(wp), parameter :: pw = merge(log(epsilon(1.0_wp))/log(real(epsilon(1.0_real64), wp)), 1.0_wp, &
                                     epsilon(1.0_wp) > real(epsilon(1.0_real64), wp))

contains

   elemental function sc(x) result(xs)
      !! scale a tolerance specified for real64 to the real kind used (see `pw`)
      real(wp), intent(in) :: x
      real(wp) :: xs
      xs = x**pw
   end function sc

   subroutine deriv(me, t, x, xdot)
      !! derivative routine for all the test problems
      class(ddeabm_class), intent(inout) :: me
      real(wp), intent(in) :: t
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: xdot
      n_feval = n_feval + 1
      select case (iprob)
      case (p_decay); xdot = -x
      case (p_osc); xdot = [x(2), -x(1)]
      case (p_cos); xdot = cos(t)
      case (p_zero); xdot = 0.0_wp
      case (p_const); xdot = 1.0_wp
      case (p_t); xdot = t
      case (p_blowup); xdot = x**2
      case (p_stiff); xdot = -1000.0_wp*(x - cos(t))
      case (p_vdp); xdot = [x(2), 5.0_wp*(1.0_wp - x(1)**2)*x(2) - x(1)]
      case (p_stop)
         if (t > t_stop) then
            call me%stop_integration()
            xdot = 0.0_wp
         else
            xdot = [x(2), -x(1)]
         end if
      end select
   end subroutine deriv

   subroutine report(me, t, x)
      !! records the reported points
      class(ddeabm_class), intent(inout) :: me
      real(wp), intent(in) :: t
      real(wp), dimension(:), intent(in) :: x
      n_reports = n_reports + 1
      t_reports(n_reports) = t
      x_reports(n_reports) = x(1)
   end subroutine report

   subroutine g_x1(me, t, x, g)
      !! scalar event: `x(1) = 0`
      class(ddeabm_with_event_class), intent(inout) :: me
      real(wp), intent(in) :: t
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: g
      g = x(1)
   end subroutine g_x1



   subroutine g_pulse(me, t, x, g)
      !! scalar event with two close roots at `t = 1` and `t = 1.001`
      !! (negative only between them). It is scaled so the slope at the
      !! roots is about +/- 1000, for the function value tolerance.
      class(ddeabm_with_event_class), intent(inout) :: me
      real(wp), intent(in) :: t
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: g
      g = 1.0e6_wp*(t - 1.0_wp)*(t - 1.001_wp)
   end subroutine g_pulse

   logical function bracket_decreasing(me, t1, t2, x1, x2, g1, g2) result(bracketed)
      !! only accept roots where `g` is decreasing
      class(ddeabm_with_event_class), intent(inout) :: me
      real(wp), intent(in) :: t1, t2, g1, g2
      real(wp), dimension(:), intent(in) :: x1, x2
      bracketed = g1 > 0.0_wp
   end function bracket_decreasing

   subroutine g_vec(me, t, x, g, ig)
      !! vector events: `x(1) = 0` and `x(2) = 0`
      class(ddeabm_with_event_class_vec), intent(inout) :: me
      real(wp), intent(in) :: t
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: g
      integer, intent(in), optional :: ig
      if (present(ig)) then
         g = 0.0_wp
         g(ig) = x(ig)    ! only compute the one requested
      else
         g = x(1:2)
      end if
   end subroutine g_vec

   logical function bracket_vec_decreasing(me, ig, t1, t2, x1, x2, g1, g2) result(bracketed)
      !! only accept roots where `g(ig)` is decreasing
      class(ddeabm_with_event_class_vec), intent(inout) :: me
      integer, intent(in), optional :: ig
      real(wp), intent(in) :: t1, t2, g1, g2
      real(wp), dimension(:), intent(in) :: x1, x2
      bracketed = g1 > 0.0_wp
   end function bracket_vec_decreasing

   subroutine check(name, ok)
      !! record the result of a check
      character(len=*), intent(in) :: name
      logical, intent(in) :: ok
      n_checks = n_checks + 1
      if (ok) then
         write (*, '(A)') '  PASS: '//name
      else
         write (*, '(A)') '  FAIL: '//name
         n_failed = n_failed + 1
      end if
   end subroutine check

   logical function success(idid)
      !! if `idid` indicates a successful integration to `tout`
      integer, intent(in) :: idid
      success = idid == 2 .or. idid == 3
   end function success

end module ddeabm_unit_test_module

program ddeabm_unit_test

   use ddeabm_unit_test_module

   implicit none

   type(ddeabm_class) :: s, s2
   type(ddeabm_with_event_class) :: se
   type(ddeabm_with_event_class_vec) :: sv
   real(wp) :: t, y1(1), y(2), yi(2), gval, gv(2), tout, err, t_root, x_ref(2)
   real(wp), dimension(:), allocatable :: t_ref, xr_ref   !! reference reported points
   real(wp) :: hmax  !! maximum step size
   integer, dimension(ddeabm_max_order_limit) :: nf  !! number of evaluations for each maximum order
   integer :: n1
   real(wp) :: tolk, errk
   real(wp), dimension(3) :: errs
   real(wp), dimension(3), parameter :: tols = [1.0e-6_wp**pw, 1.0e-9_wp**pw, 1.0e-12_wp**pw]
   real(wp), dimension(2), parameter :: y0 = [1.0_wp, 0.0_wp]
   real(wp), parameter :: tol = 1.0e-11_wp**pw
   real(wp), parameter :: etol = 1.0e-8_wp**pw
   integer :: idid, i, k, istat
   logical :: ok

   write (*, '(A)') ''
   write (*, '(A)') '---------------------------'
   write (*, '(A)') ' ddeabm_unit_test'
   write (*, '(A)') '---------------------------'

   !===================================================================
   write (*, '(/A)') 'Integrator accuracy against analytic solutions'

   iprob = p_decay
   call s%initialize(1, 10000, deriv, [tol], [tol])
   call s%first_call()
   t = 0.0_wp; y1 = 1.0_wp
   call s%integrate(t, y1, 5.0_wp, idid=idid)
   call check('y''=-y: y(5) = exp(-5)', success(idid) .and. t == 5.0_wp .and. &
              abs(y1(1) - exp(-5.0_wp)) < etol)

   iprob = p_osc
   call s%initialize(2, 10000, deriv, [tol], [tol])
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 10.0_wp, idid=idid)
   call check('x''''=-x: x(10) = cos(10), v(10) = -sin(10)', success(idid) .and. &
              all(abs(y - [cos(10.0_wp), -sin(10.0_wp)]) < etol))

   iprob = p_cos
   call s%initialize(1, 10000, deriv, [tol], [tol])
   call s%first_call()
   t = 0.0_wp; y1 = 0.0_wp
   call s%integrate(t, y1, 3.0_wp, idid=idid)
   call check('y''=cos(t) (nonautonomous): y(3) = sin(3)', success(idid) .and. &
              abs(y1(1) - sin(3.0_wp)) < etol)

   !===================================================================
   write (*, '(/A)') 'Convergence: error decreases with the tolerance'
   iprob = p_osc
   do i = 1, size(tols)
      call s%initialize(2, 100000, deriv, [tols(i)], [tols(i)])
      call s%first_call()
      t = 0.0_wp; y = y0
      call s%integrate(t, y, 20.0_wp, idid=idid)
      errs(i) = maxval(abs(y - [cos(20.0_wp), -sin(20.0_wp)]))
   end do
   write (*, '(A,3(1X,E10.3))') '  errors for tol = 1e-6, 1e-9, 1e-12:', errs
   call check('errors decrease monotonically', errs(2) < errs(1) .and. errs(3) < errs(2))
   call check('error is within 100 times the tolerance', all(errs < 100.0_wp*tols))

   !===================================================================
   write (*, '(/A)') 'Direction of integration and continuation calls'
   iprob = p_osc
   call s%initialize(2, 10000, deriv, [tol], [tol])
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, -5.0_wp, idid=idid)
   call check('backward integration: x(-5) = cos(-5)', success(idid) .and. t == -5.0_wp .and. &
              abs(y(1) - cos(-5.0_wp)) < etol)

   call s%first_call()
   t = 0.0_wp; y = y0
   ok = .true.
   do i = 1, 10   ! continuation calls (interval mode), without restarting
      tout = real(i, wp)
      call s%integrate(t, y, tout, idid=idid)
      ok = ok .and. success(idid) .and. t == tout .and. abs(y(1) - cos(tout)) < etol
   end do
   call check('successive tout values without restarting', ok)

   call s%first_call()   ! restart with the same object gives the same answer
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 10.0_wp, idid=idid)
   call s2%initialize(2, 10000, deriv, [tol], [tol])
   call s2%first_call()
   t = 0.0_wp; yi = y0
   call s2%integrate(t, yi, 10.0_wp, idid=idid)
   call check('restart (first_call) is equivalent to a new object', all(y == yi))

   !===================================================================
   write (*, '(/A)') 'Tolerances'
   iprob = p_osc
   call s%initialize(2, 10000, deriv, [tol, tol], [tol, tol])   ! vector tolerances
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 10.0_wp, idid=idid)
   call check('vector tolerances: x(10) = cos(10)', success(idid) .and. abs(y(1) - cos(10.0_wp)) < etol)

   call s%initialize(2, 10000, deriv, [tol], [tol, tol])   ! size-1 rtol with vector atol
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 10.0_wp, idid=idid)
   call check('size-1 rtol with vector atol: x(10) = cos(10)', success(idid) .and. &
              abs(y(1) - cos(10.0_wp)) < etol)

   call s%initialize(2, 10000, deriv, [sc(1.0e-6_wp)], [0.0_wp])   ! pure relative error
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 1.0_wp, idid=idid)
   call check('atol=0 with a zero solution component: idid=-3', idid == -3 .and. t == 0.0_wp)

   call s%initialize(2, 10000, deriv, [sc(1.0e-6_wp)], [0.0_wp])
   call s%first_call()
   t = 0.5_wp; y = [cos(0.5_wp), -sin(0.5_wp)]
   call s%integrate(t, y, 1.0_wp, idid=idid)
   call check('atol=0 (pure relative error) with nonzero solution', success(idid) .and. &
              abs(y(1) - cos(1.0_wp)) < sc(1.0e-5_wp))

   !===================================================================
   write (*, '(/A)') 'Initial step size modes'
   iprob = p_osc
   do k = 1, 3
      if (k == 3) then
         call s%initialize(2, 10000, deriv, [tol], [tol], initial_step_mode=k, initial_step_size=sc(1.0e-3_wp))
      else
         call s%initialize(2, 10000, deriv, [tol], [tol], initial_step_mode=k)
      end if
      call s%first_call()
      t = 0.0_wp; y = y0
      call s%integrate(t, y, 10.0_wp, idid=idid)
      call check('initial_step_mode='//achar(iachar('0') + k)//': x(10) = cos(10)', &
                 success(idid) .and. abs(y(1) - cos(10.0_wp)) < etol)
   end do

   !===================================================================
   write (*, '(/A)') 'Starting step algorithm (dhstrt) special cases'
   ! y'=0: zero first and second derivative bounds
   iprob = p_zero
   call s%initialize(1, 10000, deriv, [tol], [tol])
   call s%first_call()
   t = 0.0_wp; y1 = 3.0_wp
   call s%integrate(t, y1, 10.0_wp, idid=idid)
   call check('y''=0: y stays constant (to roundoff)', success(idid) .and. &
              abs(y1(1) - 3.0_wp) <= 10.0_wp*spacing(3.0_wp))
   ! y'=1: zero second derivative bound, nonzero first derivative bound
   iprob = p_const
   call s%first_call()
   t = 0.0_wp; y1 = 0.0_wp
   call s%integrate(t, y1, 10.0_wp, idid=idid)
   call check('y''=1: y(10) = 10', success(idid) .and. abs(y1(1) - 10.0_wp) < etol)
   ! y'=t, y(0)=0: zero initial derivative (null perturbation vector)
   iprob = p_t
   call s%first_call()
   t = 0.0_wp; y1 = 0.0_wp
   call s%integrate(t, y1, 2.0_wp, idid=idid)
   call check('y''=t from y''(0)=0: y(2) = 2', success(idid) .and. abs(y1(1) - 2.0_wp) < etol)

   ! a user initial step that is too small for the machine precision:
   iprob = p_osc
   call s%initialize(2, 10000, deriv, [tol], [tol], initial_step_mode=3, initial_step_size=epsilon(1.0_wp)**2)
   call s%first_call()
   t = 1.0_wp; y = [cos(1.0_wp), -sin(1.0_wp)]
   call s%integrate(t, y, 2.0_wp, idid=idid)
   call check('initial step too small: idid=-2', idid == -2 .and. t == 1.0_wp)
   call s%resume_integration()   ! the step size has been increased
   call s%integrate(t, y, 2.0_wp, idid=idid)
   call check('resumed: x(2) = cos(2)', success(idid) .and. abs(y(1) - cos(2.0_wp)) < etol)

   !===================================================================
   write (*, '(/A)') 'Order changes (Van der Pol, mu=5)'
   ! the solution has rapid transitions, so the order is raised and lowered.
   ! compare with a tighter tolerance solution.
   iprob = p_vdp
   call s%initialize(2, 100000, deriv, [sc(1.0e-13_wp)], [sc(1.0e-13_wp)])
   call s%first_call()
   t = 0.0_wp; x_ref = [2.0_wp, 0.0_wp]
   call s%integrate(t, x_ref, 20.0_wp, idid=idid)
   call s%initialize(2, 100000, deriv, [sc(1.0e-8_wp)], [sc(1.0e-8_wp)])
   call s%first_call()
   t = 0.0_wp; y = [2.0_wp, 0.0_wp]
   call s%integrate(t, y, 20.0_wp, idid=idid)
   err = maxval(abs(y - x_ref))
   write (*, '(A,E10.3)') '  difference from reference solution: ', err
   call check('agrees with the reference solution', success(idid) .and. err < sc(1.0e-5_wp))
   ! interpolated output (each point is interpolated within a step, after order changes):
   call s%initialize(2, 100000, deriv, [sc(1.0e-13_wp)], [sc(1.0e-13_wp)], report=report)
   call s%first_call()
   t = 0.0_wp; y = [2.0_wp, 0.0_wp]; n_reports = 0
   call s%integrate(t, y, 20.0_wp, idid=idid, integration_mode=2, tstep=0.01_wp)
   t_ref = t_reports(1:n_reports)
   xr_ref = x_reports(1:n_reports)
   call s%initialize(2, 100000, deriv, [sc(1.0e-9_wp)], [sc(1.0e-9_wp)], report=report)
   call s%first_call()
   t = 0.0_wp; y = [2.0_wp, 0.0_wp]; n_reports = 0
   call s%integrate(t, y, 20.0_wp, idid=idid, integration_mode=2, tstep=0.01_wp)
   ok = n_reports == size(t_ref)
   if (ok) ok = all(t_reports(1:n_reports) == t_ref) .and. &
                maxval(abs(x_reports(1:n_reports) - xr_ref)) < 1000.0_wp*sc(1.0e-9_wp)
   call check('interpolated output at 2001 points agrees with the reference', success(idid) .and. ok)

   !===================================================================
   write (*, '(/A)') 'Maximum order'
   iprob = p_osc
   ! the default is 12, and it gives identical results to max_order=12:
   call s%initialize(2, 100000, deriv, [tol], [tol])
   call s%first_call()
   t = 0.0_wp; y = y0; n_feval = 0
   call s%integrate(t, y, 10.0_wp, idid=idid)
   n1 = n_feval; yi = y
   call s%initialize(2, 100000, deriv, [tol], [tol], max_order=12)
   call s%first_call()
   t = 0.0_wp; y = y0; n_feval = 0
   call s%integrate(t, y, 10.0_wp, idid=idid)
   call check('default (12): identical results to max_order=12', ddeabm_default_max_order == 12 .and. &
              n_feval == n1 .and. all(y == yi))

   ! every allowed maximum order gives an accurate solution:
   tolk = sc(1.0e-10_wp)
   ok = .true.
   do k = 1, ddeabm_max_order_limit
      call s%initialize(2, 10000000, deriv, [tolk], [tolk], max_order=k)
      call s%first_call()
      t = 0.0_wp; y = y0; n_feval = 0
      call s%integrate(t, y, 1.0_wp, idid=idid)
      nf(k) = n_feval
      errk = maxval(abs(y - [cos(1.0_wp), -sin(1.0_wp)]))
      ok = ok .and. success(idid) .and. errk < 1000.0_wp*tolk
   end do
   write (*, '(A,5(1X,I0))') '  evaluations for max_order = 1, 2, 4, 8, 12:', nf([1, 2, 4, 8, 12])
   call check('max_order = 1 to ddeabm_max_order_limit: accurate solution', ok)
   ! (on this short interval, orders above 8 are not used)
   call check('lower maximum orders need more evaluations (1 > 2 > 4 > 8 >= 12)', &
              nf(1) > nf(2) .and. nf(2) > nf(4) .and. nf(4) > nf(8) .and. nf(8) >= nf(12))

   ! for higher precision and tight tolerances, a maximum order above 12 reduces the work:
   if (precision(1.0_wp) >= 30) then
      tolk = 1.0e-30_wp
      do k = 12, 16, 4
         call s%initialize(2, 10000000, deriv, [tolk], [tolk], max_order=k)
         call s%first_call()
         t = 0.0_wp; y = y0; n_feval = 0
         call s%integrate(t, y, 20.0_wp, idid=idid)
         nf(k) = n_feval
         errk = maxval(abs(y - [cos(20.0_wp), -sin(20.0_wp)]))
         call check('tol=1e-30, max_order='//achar(iachar('0') + k/10)//achar(iachar('0') + mod(k, 10))// &
                    ': accurate solution', success(idid) .and. errk < 1.0e-27_wp)
      end do
      write (*, '(A,2(1X,I0))') '  evaluations for max_order = 12, 16:', nf(12), nf(16)
      call check('tol=1e-30: max_order=16 needs fewer evaluations than 12', nf(16) < nf(12))
   end if

   ! stiffness detection is not done for max_order <= 4:
   iprob = p_stiff
   call s%initialize(1, 500, deriv, [sc(1.0e-3_wp)], [sc(1.0e-3_wp)], max_order=4)
   call s%first_call()
   t = 0.0_wp; y1 = 1.0_wp
   ok = .true.
   do i = 1, 1000
      call s%integrate(t, y1, 10.0_wp, idid=idid)
      if (idid == -4) ok = .false.
      if (idid /= -1 .and. idid /= -4) exit
      call s%resume_integration()
   end do
   call check('max_order=4, stiff problem: stiffness is not reported', ok .and. success(idid) .and. &
              abs(y1(1) - cos(10.0_wp)) < sc(1.0e-2_wp))

   ! event finding and interpolation with the lowest and highest maximum orders:
   iprob = p_osc
   tolk = sc(1.0e-8_wp)
   call se%initialize_event(2, 10000000, deriv, [tolk], [tolk], g_x1, sc(1.0e-12_wp), max_order=1)
   call se%first_call()
   t = 0.0_wp; y = y0
   call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval, tstep=0.3_wp)
   call check('max_order=1, scalar events: root at pi/2', idid == 1000 .and. &
              abs(t - pi/2.0_wp) < 1000.0_wp*tolk)
   call sv%initialize_event(2, 10000000, deriv, [tol], [tol], 2, g_vec, [sc(1.0e-12_wp)], &
                            max_order=ddeabm_max_order_limit)
   call sv%first_call()
   t = 0.1_wp; y = [cos(t), -sin(t)]
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv)
   call sv%interpolate(t - 0.01_wp, yi, istat)
   call check('max_order=20, vector events: root at pi/2, interpolation', idid == 1001 .and. &
              abs(t - pi/2.0_wp) < etol .and. istat == 0 .and. &
              abs(yi(1) - cos(t - 0.01_wp)) < etol)

   !===================================================================
   write (*, '(/A)') 'Stiffness detection and work limit'
   ! (the work limit, idid=-1, can be reached before the stiffness is detected)
   iprob = p_stiff
   call s%initialize(1, 500, deriv, [sc(1.0e-3_wp)], [sc(1.0e-3_wp)])
   call s%first_call()
   t = 0.0_wp; y1 = 1.0_wp
   call s%integrate(t, y1, 10.0_wp, idid=idid)
   call check('stiff problem: work limit reached (idid=-1 or -4)', idid == -1 .or. idid == -4)
   ok = .false.   ! if stiffness was detected
   do i = 1, 100
      if (idid == -4) ok = .true.
      if (idid /= -1 .and. idid /= -4) exit
      call s%resume_integration()
      call s%integrate(t, y1, 10.0_wp, idid=idid)
   end do
   call check('stiff problem: stiffness detected (idid=-4)', ok)
   call check('stiff problem: resumed to completion, y(10) = cos(10)', success(idid) .and. &
              abs(y1(1) - cos(10.0_wp)) < sc(1.0e-2_wp))

   !===================================================================
   write (*, '(/A)') 'Step failures near a singularity (y''=y^2, blows up at t=1)'
   iprob = p_blowup
   call s%initialize(1, 100000, deriv, [sc(1.0e-10_wp)], [sc(1.0e-10_wp)])
   call s%first_call()
   t = 0.0_wp; y1 = 1.0_wp
   call s%integrate(t, y1, 2.0_wp, idid=idid)
   call check('integration stops at the singularity (idid<0)', idid < 0 .and. abs(t - 1.0_wp) < sc(1.0e-6_wp))
   call check('solution is very large', y1(1) > 0.1_wp/sc(1.0e-6_wp))

   !===================================================================
   write (*, '(/A)') 'tout = t on the first call'
   iprob = p_osc
   call s%initialize(2, 10000, deriv, [tol], [tol])
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 0.0_wp, idid=idid)
   call check('idid=2 and no change to y', idid == 2 .and. t == 0.0_wp .and. all(y == y0))

   !===================================================================
   write (*, '(/A)') 'tstop'
   iprob = p_osc
   call s%initialize(2, 10000, deriv, [tol], [tol])
   ok = .true.
   do k = 1, 20   ! (reaching tstop may also use the extrapolation branch)
      call s%first_call()
      t = 0.0_wp; y = y0; n_feval = 0
      tout = 0.1_wp*k
      call s%integrate(t, y, tout, tstop=tout, idid=idid)
      ok = ok .and. success(idid) .and. t == tout .and. abs(y(1) - cos(tout)) < etol
   end do
   call check('integrate to tout = tstop', ok)
   call s%initialize(2, 10000, deriv, [tol], [tol], report=report)
   ok = .true.
   do k = 1, 20   ! intermediate output: steps end exactly at (or very near) tstop
      call s%first_call()
      t = 0.0_wp; y = y0
      tout = 0.37_wp*k
      call s%integrate(t, y, tout, tstop=tout, idid=idid, integration_mode=2)
      ok = ok .and. success(idid) .and. t == tout .and. abs(y(1) - cos(tout)) < etol
   end do
   call check('intermediate output to tout = tstop', ok)
   ! moving tstop: when within roundoff of tstop, the code extrapolates
   call s%first_call()
   t = 0.0_wp; y = y0
   tout = 1.0_wp - 2.0_wp*spacing(1.0_wp)
   call s%integrate(t, y, tout, tstop=tout, idid=idid)
   call s%integrate(t, y, 1.0_wp, tstop=1.0_wp, idid=idid)
   call check('tstop moved by roundoff: extrapolated to tout', idid == 3 .and. t == 1.0_wp .and. &
              abs(y(1) - cos(1.0_wp)) < etol)
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 2.0_wp, tstop=1.0_wp, idid=idid)
   call check('tout beyond tstop: idid=-33', idid == -33)

   !===================================================================
   write (*, '(/A)') 'Maximum step size'
   iprob = p_osc
   hmax = 0.05_wp
   ! without a maximum (the default), steps are larger than hmax:
   call s%initialize(2, 10000, deriv, [tol], [tol], report=report)
   call s%first_call()
   t = 0.0_wp; y = y0; n_reports = 0
   call s%integrate(t, y, 5.0_wp, idid=idid, integration_mode=2)
   call check('no maximum: some steps are larger than hmax', success(idid) .and. &
              maxval(t_reports(2:n_reports) - t_reports(1:n_reports - 1)) > hmax)
   ! with a maximum:
   do k = 1, 3   ! each initial step mode (mode 3 with a too-large initial step)
      if (k == 3) then
         call s%initialize(2, 10000, deriv, [tol], [tol], report=report, initial_step_mode=k, &
                           initial_step_size=1.0_wp, max_step_size=hmax)
      else
         call s%initialize(2, 10000, deriv, [tol], [tol], report=report, initial_step_mode=k, &
                           max_step_size=hmax)
      end if
      call s%first_call()
      t = 0.0_wp; y = y0; n_reports = 0
      call s%integrate(t, y, 5.0_wp, idid=idid, integration_mode=2)
      ok = success(idid) .and. n_reports > 1
      if (ok) ok = maxval(t_reports(2:n_reports) - t_reports(1:n_reports - 1)) <= hmax*(1.0_wp + sc(1.0e-12_wp))
      call check('initial_step_mode='//achar(iachar('0') + k)//': all steps <= max_step_size', ok)
      call check('initial_step_mode='//achar(iachar('0') + k)//': x(5) = cos(5)', &
                 abs(y(1) - cos(5.0_wp)) < etol)
   end do
   ! backward:
   call s%first_call()
   t = 0.0_wp; y = y0; n_reports = 0
   call s%integrate(t, y, -5.0_wp, idid=idid, integration_mode=2)
   ok = success(idid) .and. n_reports > 1
   if (ok) ok = maxval(t_reports(1:n_reports - 1) - t_reports(2:n_reports)) <= hmax*(1.0_wp + sc(1.0e-12_wp))
   call check('backward: all steps <= max_step_size, x(-5) = cos(-5)', ok .and. &
              abs(y(1) - cos(-5.0_wp)) < etol)

   ! events: two close roots within one step are not detected (no sign change),
   ! unless the step size is limited to less than their separation:
   call se%initialize_event(2, 100000, deriv, [tol], [tol], g_pulse, sc(1.0e-12_wp))
   call se%first_call()
   t = 0.0_wp; y = y0
   call se%integrate_to_event(t, y, 2.0_wp, idid=idid, gval=gval)
   call check('events, no maximum: close roots are stepped over', success(idid) .and. t == 2.0_wp)
   call se%initialize_event(2, 100000, deriv, [tol], [tol], g_pulse, sc(1.0e-12_wp), &
                            max_step_size=4.0e-4_wp)
   call se%first_call()
   t = 0.0_wp; y = y0
   call se%integrate_to_event(t, y, 2.0_wp, idid=idid, gval=gval)
   call check('events, max_step_size < root separation: first root at t=1', &
              idid == 1000 .and. abs(t - 1.0_wp) < etol)
   call se%integrate_to_event(t, y, 2.0_wp, idid=idid, gval=gval, continue=.true.)
   call check('events, max_step_size < root separation: second root at t=1.001', &
              idid == 1000 .and. abs(t - 1.001_wp) < etol)

   call sv%initialize_event(2, 10000, deriv, [tol], [tol], 2, g_vec, [sc(1.0e-12_wp)], &
                            report=report, max_step_size=hmax)
   call sv%first_call()
   t = 0.1_wp; y = [cos(t), -sin(t)]; n_reports = 0
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, integration_mode=2)
   ok = idid == 1001 .and. abs(t - pi/2.0_wp) < etol .and. n_reports > 1
   if (ok) ok = maxval(t_reports(2:n_reports) - t_reports(1:n_reports - 1)) <= hmax*(1.0_wp + sc(1.0e-12_wp))
   call check('vector events: all steps <= max_step_size, root at pi/2', ok)

   !===================================================================
   write (*, '(/A)') 'Dense output (interpolate)'
   iprob = p_osc
   call s%initialize(2, 10000, deriv, [tol], [tol], report=report)
   call s%first_call()
   t = 0.0_wp; y = y0; n_reports = 0
   call s%integrate(t, y, 10.0_wp, idid=idid, integration_mode=2)  ! one step at a time
   ! interpolate within the last step:
   ok = n_reports > 2
   if (ok) then
      do i = 0, 10
         tout = t_reports(n_reports - 1) + (t_reports(n_reports) - t_reports(n_reports - 1))*i/10.0_wp
         call s%interpolate(tout, yi, istat)
         ok = ok .and. istat == 0 .and. all(abs(yi - [cos(tout), -sin(tout)]) < etol)
      end do
   end if
   call check('interpolated states in the last step are correct', ok)
   call s2%initialize(2, 10000, deriv, [tol], [tol])
   call s2%interpolate(1.0_wp, yi)   ! before any step (prints an error message)
   call check('interpolate before any step (no istat): returns zero state', all(yi == 0.0_wp))

   !===================================================================
   write (*, '(/A)') 'Reporting (integration_mode=2)'
   iprob = p_osc
   call s%initialize(2, 10000, deriv, [tol], [tol], report=report)
   call s%first_call()
   t = 0.0_wp; y = y0; n_reports = 0
   call s%integrate(t, y, 10.0_wp, idid=idid, integration_mode=2)
   ok = n_reports > 2 .and. t_reports(1) == 0.0_wp .and. t_reports(n_reports) == 10.0_wp
   ok = ok .and. all(t_reports(2:n_reports) > t_reports(1:n_reports - 1))
   ok = ok .and. all(abs(x_reports(1:n_reports) - cos(t_reports(1:n_reports))) < etol)
   call check('default steps: increasing times from t0 to tout, correct states', success(idid) .and. ok)

   call s%first_call()
   t = 0.0_wp; y = y0; n_reports = 0
   call s%integrate(t, y, 10.0_wp, idid=idid, integration_mode=2, tstep=0.75_wp)
   ok = n_reports == 15 .and. t_reports(n_reports) == 10.0_wp
   do i = 1, n_reports - 1
      ok = ok .and. abs(t_reports(i) - 0.75_wp*(i - 1)) < sc(1.0e-12_wp)
   end do
   ok = ok .and. all(abs(x_reports(1:n_reports) - cos(t_reports(1:n_reports))) < etol)
   call check('fixed step: points at t0 + k*tstep and tout, correct states', success(idid) .and. ok)

   call s%first_call()
   t = 0.0_wp; y = y0; n_reports = 0
   call s%integrate(t, y, -3.0_wp, idid=idid, integration_mode=2, tstep=1.0_wp)
   ok = n_reports == 4 .and. all(abs(t_reports(1:4) - [0.0_wp, -1.0_wp, -2.0_wp, -3.0_wp]) < sc(1.0e-12_wp))
   call check('fixed step, backward: points at 0,-1,-2,-3', success(idid) .and. ok)

   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 1.0_wp, idid=idid, integration_mode=3)
   call check('invalid integration_mode: idid=-33', idid == -33)
   call s2%initialize(2, 10000, deriv, [tol], [tol])   ! no report function
   call s2%first_call()
   t = 0.0_wp; y = y0
   call s2%integrate(t, y, 1.0_wp, idid=idid, integration_mode=2)
   call check('integration_mode=2 without a report function: idid=-33', idid == -33)

   !===================================================================
   write (*, '(/A)') 'User-triggered stop (stop_integration)'
   iprob = p_stop
   t_stop = 0.5_wp
   call s%initialize(2, 10000, deriv, [tol], [tol], report=report)
   do k = 1, 2
      call s%first_call()
      t = 0.0_wp; y = y0
      call s%integrate(t, y, 1.0_wp, idid=idid, integration_mode=k)
      call check('integrate (mode '//achar(iachar('0') + k)//'): idid=-1000, t <= t_stop', &
                 idid == -1000 .and. t <= t_stop .and. abs(y(1) - cos(t)) < etol)
   end do
   call se%initialize_event(2, 10000, deriv, [tol], [tol], g_x1, sc(1.0e-12_wp))
   call se%first_call()
   t = 0.0_wp; y = y0
   call se%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gval)
   call check('integrate_to_event: idid=-1000, last good point returned', &
              idid == -1000 .and. t <= t_stop .and. abs(y(1) - cos(t)) < etol .and. gval == y(1))
   call sv%initialize_event(2, 10000, deriv, [tol], [tol], 2, g_vec, [sc(1.0e-12_wp)])
   call sv%first_call()
   t = 0.0_wp; y = [cos(0.1_wp), -sin(0.1_wp)]
   call sv%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gv)
   call check('vector integrate_to_event: idid=-1000, last good point returned', &
              idid == -1000 .and. t <= t_stop .and. abs(y(1) - cos(t + 0.1_wp)) < etol .and. all(gv == y))
   t_stop = huge(1.0_wp)

   !===================================================================
   write (*, '(/A)') 'Input validation in ddes'
   iprob = p_osc
   call s%initialize(2, 10000, deriv, [-sc(1.0e-6_wp)], [tol])
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 1.0_wp, idid=idid)
   call check('negative rtol: idid=-33', idid == -33)
   call s%initialize(2, 10000, deriv, [tol], [-sc(1.0e-6_wp)])
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 1.0_wp, idid=idid)
   call check('negative atol: idid=-33', idid == -33)

   call s%initialize(2, 10000, deriv, [tol], [tol])
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 1.0_wp, idid=idid)
   call s%integrate(t, y, 1.0_wp, idid=idid)
   call check('continuation call with t = tout: idid=-33', idid == -33)

   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 1.0_wp, idid=idid)
   t = 0.5_wp
   call s%integrate(t, y, 2.0_wp, idid=idid)
   call check('continuation call with a changed t: idid=-33', idid == -33)

   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 1.0_wp, idid=idid)
   call s%integrate(t, y, -1.0_wp, idid=idid)
   call check('change of direction without restarting: idid=-33', idid == -33)
   call s%integrate(t, y, 2.0_wp, idid=idid)
   call check('call again after invalid input (info(1)=-1): idid=-33', idid == -33)
   ok = .true.
   do i = 1, 6   ! repeated calls without fixing the problem (infinite loop detection)
      call s%integrate(t, y, 2.0_wp, idid=idid)
      ok = ok .and. idid == -33
   end do
   call check('repeated invalid calls keep returning idid=-33', ok)

   call s%initialize(0, 10000, deriv, [tol], [tol])
   call s%first_call()
   t = 0.0_wp
   call s%integrate(t, y(1:0), 1.0_wp, idid=idid)
   call check('neq=0: idid=-33', idid == -33)

   call s%destroy()   ! the derivative routine is no longer associated
   t = 0.0_wp
   call s%integrate(t, y(1:0), 1.0_wp, idid=idid)
   call check('after destroy: idid=-33', idid == -33)

   !===================================================================
   write (*, '(/A)') 'Event finding (scalar)'
   iprob = p_osc
   call se%initialize_event(2, 10000, deriv, [tol], [tol], g_x1, sc(1.0e-12_wp), report=report)
   do k = 1, 2   ! default steps and fixed step
      call se%first_call()
      t = 0.0_wp; y = y0
      if (k == 1) then
         call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval)
      else
         call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval, tstep=0.4_wp)
      end if
      call check(trim(merge('default steps', 'fixed step   ', k == 1))//': root at pi/2', &
                 idid == 1000 .and. abs(t - pi/2.0_wp) < etol .and. abs(gval) < etol .and. &
                 abs(y(2) + 1.0_wp) < etol)
   end do

   call se%first_call()   ! backward
   t = 0.0_wp; y = y0
   call se%integrate_to_event(t, y, -10.0_wp, idid=idid, gval=gval)
   call check('backward: root at -pi/2', idid == 1000 .and. abs(t + pi/2.0_wp) < etol)

   call se%first_call()   ! root at the initial time is ignored
   t = pi/2.0_wp; y = [cos(t), -sin(t)]
   call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval)
   call check('root at the initial time is ignored: next root at 3pi/2', &
              idid == 1000 .and. abs(t - 1.5_wp*pi) < etol)

   call se%first_call()   ! with tstop and vector tolerances
   call se%initialize_event(2, 10000, deriv, [tol, tol], [tol, tol], g_x1, sc(1.0e-12_wp))
   call se%first_call()
   t = 0.0_wp; y = y0
   call se%integrate_to_event(t, y, 1.0_wp, tstop=1.0_wp, idid=idid, gval=gval)
   call check('tstop and vector tolerances, no root: stops at tmax', &
              success(idid) .and. t == 1.0_wp .and. abs(gval - cos(1.0_wp)) < etol)

   call se%initialize_event(2, 10000, deriv, [tol], [tol], g_x1, sc(1.0e-12_wp), &
                            bracket=bracket_decreasing, report=report)
   call se%first_call()   ! bracket function: only decreasing crossings
   t = 0.0_wp; y = y0
   call se%integrate_to_event(t, y, 20.0_wp, idid=idid, gval=gval)
   t_root = t
   call se%integrate_to_event(t, y, 20.0_wp, idid=idid, gval=gval, continue=.true.)
   call check('bracket function: roots at pi/2 then 5pi/2 (3pi/2 rejected)', &
              abs(t_root - pi/2.0_wp) < etol .and. idid == 1000 .and. abs(t - 2.5_wp*pi) < etol)

   call se%initialize_event(2, 10000, deriv, [tol], [tol], g_x1, sc(1.0e-12_wp), report=report)
   call se%first_call()   ! continuing in the opposite direction is not allowed
   t = 0.0_wp; y = y0
   call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval)
   call se%integrate_to_event(t, y, 0.0_wp, idid=idid, gval=gval, continue=.true.)
   call check('continue in the opposite direction: idid=-33', idid == -33)

   call se%first_call()   ! continue with tmax inside the rest of the last step (default steps)
   t = 0.0_wp; y = y0
   call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval)
   t_root = t
   tout = t_root + sc(1.0e-6_wp)
   call se%integrate_to_event(t, y, tout, idid=idid, gval=gval, continue=.true.)
   call check('continue with tmax in the last step: stops at tmax', &
              idid == 3 .and. t == tout .and. abs(y(1) - cos(tout)) < etol)

   call se%first_call()   ! continue with a smaller tstep: points in the rest of the step are interpolated
   t = 0.0_wp; y = y0; n_reports = 0
   call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval, integration_mode=2, tstep=3.0_wp)
   t_root = t
   n_reports = 0
   call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval, integration_mode=2, &
                              tstep=0.25_wp, continue=.true.)
   ok = idid == 1000 .and. abs(t - 1.5_wp*pi) < etol .and. n_reports > 4
   do i = 1, 4   ! (3-pi/2)/0.25 = 5.7 points in the rest of the first step
      ok = ok .and. abs(t_reports(i) - (t_root + 0.25_wp*i)) < sc(1.0e-12_wp) .and. &
           abs(x_reports(i) - cos(t_reports(i))) < etol
   end do
   ok = ok .and. all(t_reports(2:n_reports) > t_reports(1:n_reports - 1))
   call check('continue with smaller tstep: grid from the root, correct states, next root', ok)

   ! the root accuracy does not depend on tstep (the root is located within
   ! an integration step, where the interpolating polynomial is valid):
   ok = .true.
   do k = 1, 8
      call se%first_call()
      t = 0.0_wp; y = y0
      call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval, tstep=0.5_wp*k)
      ok = ok .and. idid == 1000 .and. abs(t - pi/2.0_wp) < sc(1.0e-10_wp) .and. &
           abs(y(1) - cos(t)) < sc(1.0e-10_wp)
   end do
   call check('fixed step: root accuracy does not depend on tstep (0.5 to 4)', ok)

   call se%initialize_event(2, 10000, deriv, [tol], [tol], g_x1, sc(1.0e-12_wp))
   call se%first_call()
   t = 0.0_wp; y = y0
   call se%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gval, integration_mode=3)
   call check('invalid integration_mode: idid=-33', idid == -33)
   call se%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gval, tstep=0.0_wp)
   call check('tstep=0: idid=-33', idid == -33)
   t = 1.8_wp; y = [cos(t), -sin(t)]
   call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval, tstep=0.01_wp*spacing(t))
   call check('tstep lost in roundoff: idid=-33', idid == -33 .and. t == 1.8_wp)

   call se%first_call()   ! when continuing, tmax must not be beyond tstop
   t = 0.0_wp; y = y0
   call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval)
   t_root = t
   call se%integrate_to_event(t, y, t_root + 0.002_wp, tstop=t_root + 0.001_wp, idid=idid, gval=gval, &
                              continue=.true.)
   call check('continue with tmax beyond tstop: idid=-33', idid == -33)
   call se%integrate_to_event(t, y, t_root + 0.001_wp, tstop=t_root + 0.001_wp, idid=idid, gval=gval, &
                              continue=.true.)
   call check('continue with tmax = tstop: stops at tmax', success(idid) .and. &
              t == t_root + 0.001_wp .and. abs(y(1) - cos(t)) < etol)
   call se%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gval, integration_mode=2)
   call check('integration_mode=2 without a report function: idid=-33', idid == -33)

   !===================================================================
   write (*, '(/A)') 'Event finding (vector)'
   iprob = p_osc
   call sv%initialize_event(2, 10000, deriv, [tol], [tol], 2, g_vec, [sc(1.0e-12_wp), sc(1.0e-12_wp)], &
                            report=report)
   call sv%first_call()   ! earliest root of all the functions, in order
   t = 0.1_wp; y = [cos(t), -sin(t)]
   ok = .true.
   do i = 1, 4   ! roots: x(1) at pi/2, x(2) at pi, x(1) at 3pi/2, x(2) at 2pi
      call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, continue=i > 1)
      ok = ok .and. idid == 1000 + 2 - mod(i, 2) .and. abs(t - i*pi/2.0_wp) < etol .and. &
           abs(gv(idid - 1000)) < etol
   end do
   call check('roots of both functions found in order (pi/2, pi, 3pi/2, 2pi)', ok)

   call sv%first_call()   ! backward
   t = -0.1_wp; y = [cos(t), -sin(t)]
   call sv%integrate_to_event(t, y, -10.0_wp, idid=idid, gval=gv, integration_mode=2)
   call check('backward: first root is x(1) at -pi/2', idid == 1001 .and. abs(t + pi/2.0_wp) < etol)

   call sv%first_call()   ! interrupted task (maxnum)
   call sv%initialize_event(2, 5, deriv, [tol], [tol], 2, g_vec, [sc(1.0e-12_wp)])
   call sv%first_call()
   t = 0.1_wp; y = [cos(t), -sin(t)]
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv)
   call check('interrupted (idid=-1): gval = g(t,y)', idid == -1 .and. all(gv == y))
   do i = 1, 1000
      call sv%resume_integration()
      call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv)
      if (idid /= -1) exit
   end do
   call check('resumed: first root is x(1) at pi/2', idid == 1001 .and. abs(t - pi/2.0_wp) < etol)

   call sv%initialize_event(2, 10000, deriv, [tol, tol], [tol, tol], 2, g_vec, [sc(1.0e-12_wp)], &
                            bracket=bracket_vec_decreasing)
   call sv%first_call()   ! bracket function, vector tolerances, tstop
   t = 0.1_wp; y = [cos(t), -sin(t)]
   call sv%integrate_to_event(t, y, 10.0_wp, tstop=10.0_wp, idid=idid, gval=gv)
   t_root = t
   call sv%integrate_to_event(t, y, 10.0_wp, tstop=10.0_wp, idid=idid, gval=gv, continue=.true.)
   ! x(1)=cos decreases through 0 at pi/2; x(2)=-sin increases through 0 at pi (rejected),
   ! decreases through 0 at 2pi
   call check('bracket function: x(1) at pi/2, then x(2) at 2pi', &
              abs(t_root - pi/2.0_wp) < etol .and. idid == 1002 .and. abs(t - 2.0_wp*pi) < etol)

   call sv%initialize_event(2, 10000, deriv, [tol], [tol], 2, g_vec, [sc(1.0e-12_wp), sc(1.0e-12_wp), 1.0_wp])
   call sv%first_call()   ! wrong size root_tol: the first element is used
   t = 0.1_wp; y = [cos(t), -sin(t)]
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv)
   call check('wrong size root_tol: first element used, root found', &
              idid == 1001 .and. abs(t - pi/2.0_wp) < etol)

   call sv%first_call()   ! continue with tmax inside the rest of the last step
   t = 0.1_wp; y = [cos(t), -sin(t)]
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv)
   tout = t + sc(1.0e-6_wp)
   call sv%integrate_to_event(t, y, tout, idid=idid, gval=gv, continue=.true.)
   call check('continue with tmax in the last step: stops at tmax', &
              idid == 3 .and. t == tout .and. abs(y(1) - cos(tout)) < etol)
   call sv%integrate_to_event(t, y, 0.0_wp, idid=idid, gval=gv, continue=.true.)
   call check('continue in the opposite direction: idid=-33', idid == -33)

   call sv%first_call()
   t = 0.0_wp; y = y0
   call sv%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gv, integration_mode=3)
   call check('invalid integration_mode: idid=-33', idid == -33)
   call sv%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gv, tstep=0.0_wp)
   call check('tstep=0: idid=-33', idid == -33)
   t = 1.8_wp; y = [cos(t), -sin(t)]
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, tstep=0.01_wp*spacing(t))
   call check('tstep lost in roundoff: idid=-33', idid == -33 .and. t == 1.8_wp)

   call sv%first_call()   ! when continuing, tmax must not be beyond tstop
   t = 0.1_wp; y = [cos(t), -sin(t)]
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv)
   t_root = t
   call sv%integrate_to_event(t, y, t_root + 0.002_wp, tstop=t_root + 0.001_wp, idid=idid, gval=gv, &
                              continue=.true.)
   call check('continue with tmax beyond tstop: idid=-33', idid == -33)
   call sv%integrate_to_event(t, y, t_root + 0.001_wp, tstop=t_root + 0.001_wp, idid=idid, gval=gv, &
                              continue=.true.)
   call check('continue with tmax = tstop: stops at tmax', success(idid) .and. &
              t == t_root + 0.001_wp .and. abs(y(1) - cos(t)) < etol)
   call sv%initialize_event(2, 10000, deriv, [tol], [tol], 2, g_vec, [sc(1.0e-12_wp)])
   call sv%first_call()
   call sv%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gv, integration_mode=2)
   call check('integration_mode=2 without a report function: idid=-33', idid == -33)

   call sv%destroy()   ! uninitialized: no event function
   call sv%integrate_to_event(t, y(1:0), 1.0_wp, idid=idid, gval=gv(1:0))
   call check('no event function associated: idid=-33', idid == -33)
   call se%destroy()
   call se%integrate_to_event(t, y(1:0), 1.0_wp, idid=idid, gval=gval)
   call check('scalar: no event function associated: idid=-33', idid == -33)

   !===================================================================
   write (*, '(A)') ''
   if (n_failed > 0) then
      write (*, '(A,I0,A,I0,A)') 'FAILED: ', n_failed, ' of ', n_checks, ' checks failed'
      error stop 1
   else
      write (*, '(A,I0,A)') 'All ', n_checks, ' checks passed.'
   end if

end program ddeabm_unit_test
!*****************************************************************************************
