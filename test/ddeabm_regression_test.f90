!*****************************************************************************************
!>
!  Regression tests for bugs fixed in [[ddeabm_module]].
!  All tests use the harmonic oscillator
!  `x'' = -x`, `x(0) = 1`, `x'(0) = 0`, with the solution `x = cos(t)`.
!
!  The program stops with a nonzero exit code if any check fails.

module ddeabm_regression_test_module

   use ddeabm_module, wp => ddeabm_rk
   use, intrinsic :: iso_fortran_env, only: real64

   implicit none

   real(wp), parameter :: pi = acos(-1.0_wp)
   real(wp), parameter :: pi2 = pi/2.0_wp

   real(wp) :: tmax_seen = -huge(1.0_wp)   !! largest `t` passed to the derivative routine
   real(wp) :: t_event2 = 0.0_wp           !! root of the second vector event function
   integer :: n_reports = 0                !! number of reported points
   real(wp), dimension(10000) :: t_reports = 0.0_wp  !! the reported times
   integer :: n_failed = 0                 !! number of failed checks

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

   subroutine osc(me, t, x, xdot)
      !! derivative routine: harmonic oscillator
      class(ddeabm_class), intent(inout) :: me
      real(wp), intent(in) :: t
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: xdot
      tmax_seen = max(tmax_seen, t)
      xdot = [x(2), -x(1)]
   end subroutine osc

   subroutine report(me, t, x)
      !! records the reported times
      class(ddeabm_class), intent(inout) :: me
      real(wp), intent(in) :: t
      real(wp), dimension(:), intent(in) :: x
      n_reports = n_reports + 1
      t_reports(n_reports) = t
   end subroutine report

   subroutine g_scalar(me, t, x, g)
      !! scalar event: `x = 0` (roots at `pi/2 + k*pi`)
      class(ddeabm_with_event_class), intent(inout) :: me
      real(wp), intent(in) :: t
      real(wp), dimension(:), intent(in) :: x
      real(wp), intent(out) :: g
      g = x(1)
   end subroutine g_scalar

   subroutine g_vec(me, t, x, g, ig)
      !! vector events:
      !!  1. `x = 0` (roots at `pi/2 + k*pi`)
      !!  2. `t = t_event2`
      class(ddeabm_with_event_class_vec), intent(inout) :: me
      real(wp), intent(in) :: t
      real(wp), dimension(:), intent(in) :: x
      real(wp), dimension(:), intent(out) :: g
      integer, intent(in), optional :: ig
      g(1) = x(1)
      g(2) = t - t_event2
   end subroutine g_vec

   subroutine check(name, ok)
      !! record the result of a check
      character(len=*), intent(in) :: name
      logical, intent(in) :: ok
      if (ok) then
         write (*, '(A)') '  PASS: '//name
      else
         write (*, '(A)') '  FAIL: '//name
         n_failed = n_failed + 1
      end if
   end subroutine check

end module ddeabm_regression_test_module

program ddeabm_regression_test

   use ddeabm_regression_test_module

   implicit none

   type(ddeabm_class) :: s
   type(ddeabm_with_event_class) :: se
   type(ddeabm_with_event_class_vec) :: sv
   real(wp) :: t, y(2), yi(2), gval, gv(2), err1, err2, t_root
   integer :: idid, istat, i, k
   logical :: ok

   real(wp), parameter :: tol = 1.0e-10_wp**pw    !! integration tolerance
   real(wp), parameter :: etol = 1.0e-8_wp**pw    !! tolerance for checking results
   real(wp), dimension(2), parameter :: y0 = [1.0_wp, 0.0_wp]

   write (*, '(A)') ''
   write (*, '(A)') '---------------------------'
   write (*, '(A)') ' ddeabm_regression_test'
   write (*, '(A)') '---------------------------'

   !-------------------------------------------------------------------
   write (*, '(/A)') '1. stale idid value in the caller does not affect a new call'
   call s%initialize(2, 1000, osc, [tol], [tol])
   do k = 1, 2
      call s%first_call()
      t = 0.0_wp; y = y0
      idid = merge(-33, -2, k == 1)   ! stale value from a previous call
      call s%integrate(t, y, 1.0_wp, idid=idid)
      call check('stale idid='//trim(merge('-33', '-2 ', k == 1))//' -> success, t=1, y=cos(1)', &
                 (idid == 2 .or. idid == 3) .and. t == 1.0_wp .and. abs(y(1) - cos(1.0_wp)) < etol)
   end do

   !-------------------------------------------------------------------
   write (*, '(/A)') '2. vector events: a function that is zero at t0 does not mask the others'
   call sv%initialize_event(2, 1000, osc, [tol], [tol], 2, g_vec, [sc(1.0e-12_wp)])
   t_event2 = 0.0_wp   ! g(2) is zero at the initial time
   call sv%first_call()
   t = 0.0_wp; y = y0
   ! the first interval [0,2.5] contains one root of g(1):
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, tstep=2.5_wp)
   call check('root of g(1) at pi/2 found in the first interval', &
              idid == 1001 .and. abs(t - pi2) < etol)

   !-------------------------------------------------------------------
   write (*, '(/A)') '3. initial_step_mode=3 respects tout and tstop'
   call s%initialize(2, 1000, osc, [tol], [tol], initial_step_mode=3, initial_step_size=5.0_wp)
   call s%first_call()
   t = 0.0_wp; y = y0; tmax_seen = -huge(1.0_wp)
   call s%integrate(t, y, 1.0_wp, tstop=1.0_wp, idid=idid)
   call check('df not evaluated past tstop', tmax_seen <= 1.0_wp)
   call check('solution at t=1 is correct', (idid == 2 .or. idid == 3) .and. abs(y(1) - cos(1.0_wp)) < etol)

   !-------------------------------------------------------------------
   write (*, '(/A)') '4. gval is set when no root is found'
   call se%initialize_event(2, 1000, osc, [tol], [tol], g_scalar, sc(1.0e-12_wp))
   call se%first_call()
   t = 0.0_wp; y = y0; gval = -999.0_wp
   call se%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gval)
   call check('scalar (default steps): gval = cos(1)', &
              idid == 3 .and. abs(gval - cos(1.0_wp)) < etol)
   call se%first_call()
   t = 0.0_wp; y = y0; gval = -999.0_wp
   call se%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gval, tstep=0.3_wp)
   call check('scalar (fixed step): gval = cos(1)', abs(gval - cos(1.0_wp)) < etol)

   t_event2 = 100.0_wp
   call sv%first_call()
   t = 0.0_wp; y = y0; gv = -999.0_wp
   call sv%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gv)
   call check('vector: gval = [cos(1), 1-100]', &
              abs(gv(1) - cos(1.0_wp)) < etol .and. abs(gv(2) - (1.0_wp - 100.0_wp)) < etol)

   call se%initialize_event(2, 5, osc, [tol], [tol], g_scalar, sc(1.0e-12_wp))
   call se%first_call()
   t = 0.0_wp; y = y0; gval = -999.0_wp
   call se%integrate_to_event(t, y, 1.0_wp, idid=idid, gval=gval, tstep=1.0_wp)
   call check('scalar, interrupted (idid=-1): gval = g(t,y)', idid == -1 .and. gval == y(1))

   !-------------------------------------------------------------------
   write (*, '(/A)') '5. interrupted tasks can be continued with resume_integration()'
   call s%initialize(2, 5, osc, [tol], [tol])
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 20.0_wp, idid=idid)
   call check('maxnum exceeded: idid=-1', idid == -1)
   do i = 1, 1000
      call s%resume_integration()
      call s%integrate(t, y, 20.0_wp, idid=idid)
      if (idid /= -1) exit
   end do
   call check('after resuming: success, solution at t=20 is correct', &
              (idid == 2 .or. idid == 3) .and. t == 20.0_wp .and. abs(y(1) - cos(20.0_wp)) < sc(1.0e-7_wp))

   ! (start where no solution component is zero, since a pure
   !  relative error test is not possible there: idid=-3)
   call s%initialize(2, 1000, osc, [0.0_wp], [0.0_wp])
   call s%first_call()
   t = 0.5_wp; y = [cos(0.5_wp), -sin(0.5_wp)]
   call s%integrate(t, y, 1.5_wp, idid=idid)
   call check('rtol=atol=0: idid=-2', idid == -2)
   do i = 1, 10   ! (the code may also relax the tolerances again: idid=-2)
      call s%resume_integration()
      call s%integrate(t, y, 1.5_wp, idid=idid)
      if (idid /= -2) exit
   end do
   call check('after resuming: success, solution at t=1.5 is accurate', &
              (idid == 2 .or. idid == 3) .and. t == 1.5_wp .and. abs(y(1) - cos(1.5_wp)) < sc(1.0e-12_wp))

   call s%first_call()   ! restarting resets the tolerances to the user values
   t = 0.5_wp; y = [cos(0.5_wp), -sin(0.5_wp)]
   call s%integrate(t, y, 1.5_wp, idid=idid)
   call check('restart with first_call: idid=-2 again', idid == -2)

   !-------------------------------------------------------------------
   write (*, '(/A)') '6. continue: the rest of the step after a root is checked'
   call sv%initialize_event(2, 1000, osc, [tol], [tol], 2, g_vec, [sc(1.0e-12_wp)], report=report)
   do k = 1, 2   ! fixed step and default steps
      if (k == 1) then
         t_event2 = pi2 + 0.1_wp   ! both roots within one fixed step [0,2]
      else
         t_event2 = pi2 + 1.0e-4_wp  ! both roots within one default step
      end if
      call sv%first_call()
      t = 0.0_wp; y = y0; n_reports = 0
      if (k == 1) then
         call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, integration_mode=2, tstep=2.0_wp)
      else
         call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, integration_mode=2)
      end if
      call check(trim(merge('fixed  ', 'default', k == 1))//': first root: g(1) at pi/2', &
                 idid == 1001 .and. abs(t - pi2) < etol)
      t_root = t
      if (k == 1) then
         call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, integration_mode=2, &
                                    tstep=2.0_wp, continue=.true.)
      else
         call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, integration_mode=2, &
                                    continue=.true.)
      end if
      call check(trim(merge('fixed  ', 'default', k == 1))//': continue: g(2) at t_event2', &
                 idid == 1002 .and. abs(t - t_event2) < etol .and. abs(y(1) - cos(t)) < etol)
      if (k == 1) then
         call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, integration_mode=2, &
                                    tstep=2.0_wp, continue=.true.)
      else
         call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, integration_mode=2, &
                                    continue=.true.)
      end if
      call check(trim(merge('fixed  ', 'default', k == 1))//': continue: g(1) at 3pi/2', &
                 idid == 1001 .and. abs(t - 3.0_wp*pi2) < etol)
      ok = all(t_reports(2:n_reports) > t_reports(1:n_reports - 1))
      call check(trim(merge('fixed  ', 'default', k == 1))//': reported times are increasing', ok)
   end do

   ! tmax within the remainder of the last step:
   t_event2 = pi2 + 0.35_wp
   call sv%first_call()
   t = 0.0_wp; y = y0
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, tstep=2.0_wp)
   call check('tmax in last step: first root g(1) at pi/2', idid == 1001 .and. abs(t - pi2) < etol)
   call sv%integrate_to_event(t, y, 1.8_wp, idid=idid, gval=gv, tstep=2.0_wp, continue=.true.)
   call check('tmax in last step: stops at tmax=1.8', &
              idid == 3 .and. t == 1.8_wp .and. abs(y(1) - cos(1.8_wp)) < etol .and. &
              abs(gv(1) - y(1)) < etol)
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, tstep=2.0_wp, continue=.true.)
   call check('tmax in last step: then g(2) at t_event2', idid == 1002 .and. abs(t - t_event2) < etol)
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv, tstep=2.0_wp, continue=.true.)
   call check('tmax in last step: then g(1) at 3pi/2', idid == 1001 .and. abs(t - 3.0_wp*pi2) < etol)

   ! scalar, fixed step: the reporting grid continues from the root:
   call se%initialize_event(2, 1000, osc, [tol], [tol], g_scalar, sc(1.0e-12_wp), report=report)
   call se%first_call()
   t = 0.0_wp; y = y0; n_reports = 0
   call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval, integration_mode=2, tstep=0.5_wp)
   t_root = t
   n_reports = 0
   call se%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gval, integration_mode=2, &
                              tstep=0.5_wp, continue=.true.)
   call check('scalar: roots at pi/2 then 3pi/2', abs(t_root - pi2) < etol .and. &
              idid == 1000 .and. abs(t - 3.0_wp*pi2) < etol)
   call check('scalar: first report after continue is at root+tstep', &
              abs(t_reports(1) - (t_root + 0.5_wp)) < sc(1.0e-12_wp))

   !-------------------------------------------------------------------
   write (*, '(/A)') '7. vector events: the state returned at a root is at the root time'
   t_event2 = 100.0_wp
   call sv%first_call()
   t = 0.5_wp; y = [cos(0.5_wp), -sin(0.5_wp)]
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv)
   call sv%interpolate(t, yi)
   call check('y = interpolate(t) and y(1) = cos(t)', &
              idid == 1001 .and. all(y == yi) .and. abs(y(1) - cos(t)) < etol .and. gv(1) == y(1))

   !-------------------------------------------------------------------
   write (*, '(/A)') '8. interpolate before any step'
   call s%initialize(2, 1000, osc, [tol], [tol])
   call s%interpolate(1.0_wp, yi, istat)
   call check('istat = -1 and no crash', istat == -1 .and. all(yi == 0.0_wp))
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 1.0_wp, idid=idid)
   call s%interpolate(0.9_wp, yi, istat)
   call check('after a step: istat = 0 and correct', istat == 0 .and. abs(yi(1) - cos(0.9_wp)) < etol)

   !-------------------------------------------------------------------
   write (*, '(/A)') '9. tolerance and gval sizes'
   ! loose rtol for all components:
   call s%initialize(2, 1000, osc, [sc(1.0e-1_wp)], [sc(1.0e-12_wp)])
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 10.0_wp, idid=idid)
   err1 = abs(y(1) - cos(10.0_wp))
   ! vector rtol with a size-1 atol (previously only rtol(1) was used):
   call s%initialize(2, 1000, osc, [sc(1.0e-1_wp), sc(1.0e-12_wp)], [sc(1.0e-12_wp)])
   call s%first_call()
   t = 0.0_wp; y = y0
   call s%integrate(t, y, 10.0_wp, idid=idid)
   err2 = abs(y(1) - cos(10.0_wp))
   write (*, '(A,2(1X,E12.4))') '  errors (scalar rtol=0.1, vector rtol=[0.1,1e-12]):', err1, err2
   call check('vector rtol with size-1 atol is used', (idid == 2 .or. idid == 3) .and. err2 < err1*sc(1.0e-3_wp))

   call sv%initialize_event(2, 1000, osc, [tol], [tol], 2, g_vec, [sc(1.0e-12_wp)])
   call sv%first_call()
   t = 0.0_wp; y = y0
   call sv%integrate_to_event(t, y, 10.0_wp, idid=idid, gval=gv(1:1))
   call check('gval of the wrong size: idid=-33', idid == -33)

   !-------------------------------------------------------------------
   write (*, '(A)') ''
   if (n_failed > 0) then
      write (*, '(A,I0,A)') 'FAILED: ', n_failed, ' check(s) failed'
      error stop 1
   else
      write (*, '(A)') 'All checks passed.'
   end if

end program ddeabm_regression_test
!*****************************************************************************************
