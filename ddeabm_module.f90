!*****************************************************************************************    
    module ddeabm_module
!*****************************************************************************************
!****h* SLATEC/ddeabm_module
!
!  NAME
!    ddeabm_module
!
!  DESCRIPTION
!    Modern Fortran implementation of the DDEABM Adams-Bashforth algorithm.
!
!  SEE ALSO
!	http://www.netlib.org/slatec/src/
!
!  HISTORY
!	Jacob Williams : July 2014 : Created module from the SLATEC Fortran 77 code.
!
!  LICENSE
!	The original SLATEC code is a public domain work of the US Government.
!
!	The modifications are Copyright (c) 2014, Jacob Williams.
!	See https://github.com/jacobwilliams/ddeabm/blob/master/LICENSE for license terms.
!
!*****************************************************************************************

    use, intrinsic :: iso_fortran_env, wp=>real64    !double precision
    
	implicit none

    private
		    
    !parameters:    
    real(wp),parameter :: d1mach2 = huge(1.0_wp)                        ! the largest magnitude
    real(wp),parameter :: d1mach4 = radix(1.0_wp)**(1-digits(1.0_wp))   ! the largest relative spacing

	type,public :: ddeabm_class
	
		private
		
		! the number of (first order) differential equations to be integrated (>=0)
		integer :: neq = 0		
				
		!  the expense of solving the problem is monitored by counting the
		!  number of  steps attempted. when this exceeds  maxnum, the counter
		!  is reset to zero and the user is informed about possible excessive
		!  work.
		integer :: maxnum = 500		
		
		!  to define the system of first order differential equations
		!  which is to be solved.  for the given values of x and the
		!  vector  u(*)=(u(1),u(2),...,u(neq)) , the subroutine must
		!  evaluate the neq components of the system of differential
		!  equations  du/dx=df(x,u)  and store the derivatives in the
		!  array uprime(*), that is,  uprime(i) = * du(i)/dx *  for
		!  equations i=1,...,neq.
		procedure(func),pointer :: df => null()

		!work arrays:
		real(wp),dimension(:),allocatable	:: rwork
		integer 							:: lrw	= 0
		integer,dimension(:),allocatable	:: iwork
		integer 							:: liw	= 0
		
		!tolerances
		logical :: scalar_tols = .true.
		real(wp),dimension(:),allocatable :: rtol, atol	!the user input tols
		real(wp),dimension(:),allocatable :: rtol_tmp, atol_tmp	!the tols used internally
		
		!info array:
		integer,dimension(4) :: info = 0
		
	contains
	
		procedure,public :: initialize 	=> ddeabm_initialize
		procedure,public :: integrate 	=> ddeabm_wrapper
		procedure,public :: destroy 	=> destroy_ddeabm
		procedure,public :: first_call 	=> ddeabm_new_problem
		
		!support routines:					
		procedure :: ddeabm
		procedure :: ddes
		procedure :: dhstrt
		procedure :: dsteps
		
	end type ddeabm_class
	
	abstract interface
		subroutine func(me,t,x,xdot)
		import :: wp,ddeabm_class
		implicit none
		class(ddeabm_class),intent(inout) :: me
		real(wp),intent(in) :: t
		real(wp),dimension(:),intent(in) :: x
		real(wp),dimension(:),intent(out) :: xdot
		end subroutine func
	end interface
    
    contains
!*****************************************************************************************    
    
!*****************************************************************************************    
    subroutine ddeabm_new_problem(me)
!***************************************************************************************** 
!
!  Call this to indicate that a new problem is being solved (see ddeabm documentation).
!
!***************************************************************************************** 
    
    implicit none
    
  	class(ddeabm_class),intent(inout)	:: me	
    
    me%info(1) = 0
    
!*****************************************************************************************    
   end subroutine ddeabm_new_problem
!*****************************************************************************************    
   
!*****************************************************************************************    
    subroutine ddeabm_initialize(me,neq,maxnum,df,rtol,atol)
!*****************************************************************************************   
!
!  Initialize the class, and set the variables that cannot be changed during a problem. 
!
!*****************************************************************************************   
    
    implicit none
    
  	class(ddeabm_class),intent(inout)	:: me	
  	integer,intent(in)					:: neq
  	integer,intent(in) 					:: maxnum
  	procedure(func) 					:: df
  	real(wp),dimension(:),intent(in)	:: rtol
  	real(wp),dimension(:),intent(in)	:: atol
  	
  	logical :: vector_tols
  	
  	!initialize the class:
  	
  	call me%destroy()
  	
  	!number of equations:
  	
  	me%neq = neq
  	
  	!maximum number of steps:
  	
  	me%maxnum = maxnum
  	
  	!set the derivative routine pointer:
  	
  	me%df => df
  	
  	!allocate the work arrays:

 	me%lrw	= 130+21*neq
 	me%liw	= 51
 	
 	allocate(me%rwork(me%lrw))
 	me%rwork = 0.0_wp
 	
 	allocate(me%iwork(me%liw))
    me%iwork = 0
    
    !tolerances:   
	! [for now, we are considering these unchangeable, although they don't have to be]

	vector_tols = size(rtol)==neq .and. size(atol)==neq .and. neq>1
    me%scalar_tols = .not. vector_tols
    					
    if (me%scalar_tols) then
    	allocate(me%rtol(1)) ; allocate(me%rtol_tmp(1))
    	allocate(me%atol(1)) ; allocate(me%atol_tmp(1))
    	me%rtol = rtol(1)
    	me%atol = atol(1)
    else
    	allocate(me%rtol(size(rtol))) ; allocate(me%rtol_tmp(size(rtol)))
    	allocate(me%atol(size(atol))) ; allocate(me%atol_tmp(size(atol)))
    	me%rtol = rtol
    	me%atol = atol
    end if
    
!*****************************************************************************************    
    end subroutine ddeabm_initialize
!*****************************************************************************************    
        
!*****************************************************************************************    
    subroutine destroy_ddeabm(me)
!*****************************************************************************************    
    
    implicit none
    
   	class(ddeabm_class),intent(out)	:: me	
   	
!*****************************************************************************************    
    end subroutine destroy_ddeabm
!*****************************************************************************************    
            
!*****************************************************************************************    
    subroutine ddeabm_wrapper(me,t,y,tout,tstop,intermediate_steps,idid)
!*****************************************************************************************    
    
	implicit none

	class(ddeabm_class),intent(inout) 		:: me	
	real(wp),intent(inout) 					:: t
	real(wp),dimension(me%neq),intent(inout):: y
	real(wp),intent(in)						:: tout
	real(wp),intent(in),optional			:: tstop			  !not used if not present
	logical,intent(in),optional				:: intermediate_steps	 !false if not present
	integer,intent(out)						:: idid
	
	!set info array:
		
	!info(1) is set when ddeabm_new_problem is called
	
	!info(2)
	if (me%scalar_tols) then
		me%info(2) = 0
	else
		me%info(2) = 1
	end if
	
	!info(3)
	if (present(intermediate_steps)) then
		if (intermediate_steps) then
			me%info(3) = 1
		else
			me%info(3) = 0
		end if
	else
		me%info(3) = 0
	end if
			
	!info(4)
	if (present(tstop)) then
		me%info(4) = 1
		me%rwork(1) = tstop
	else
		me%info(4) = 0
	end if
	
	!make a copy of the tols, since the routine might change them:
	me%rtol_tmp = me%rtol
	me%atol_tmp = me%atol
	
   	!call the lower-level routine:
	call me%ddeabm( neq		= me%neq,&
					t		= t,&
					y		= y,&
					tout	= tout,&
					info	= me%info,&
					rtol 	= me%rtol_tmp,&
					atol	= me%atol_tmp,&
					idid	= idid,&
					rwork	= me%rwork,&
					lrw		= me%lrw,&
					iwork	= me%iwork,&
					liw		= me%liw)
	
	!Note: currently not using the recommended tols if idid=-2
					  
!*****************************************************************************************    
    end subroutine ddeabm_wrapper
!*****************************************************************************************    
    
!*****************************************************************************************    
	subroutine ddeabm (me,neq,t,y,tout,info,rtol,atol,idid,rwork,lrw,iwork,liw)
!*****************************************************************************************
!****f* ddeabm_module/ddeabm
!
!  NAME
!    ddeabm
!   
!  DESCRIPTION
!   
!	solve an initial value problem in ordinary differential
!	equations using an adams-bashforth method.
!
!   subroutine ddeabm uses the adams-bashforth-moulton
!   predictor-corrector formulas of orders one through twelve to
!   integrate a system of neq first order ordinary differential
!   equations of the form
!                         du/dx = df(x,u)
!   when the vector y(*) of initial values for u(*) at x=t is given.
!   the subroutine integrates from t to tout. it is easy to continue the
!   integration to get results at additional tout.  this is the interval
!   mode of operation.  it is also easy for the routine to return with
!   the solution at each intermediate step on the way to tout.  this is
!   the intermediate-output mode of operation.
!
! **********************************************************************
! * description of the arguments to ddeabm (an overview) *
! **********************************************************************
!
!   the parameters are
!
!      df -- this is the name of a subroutine which you provide to
!             define the differential equations.
!
!      neq -- this is the number of (first order) differential
!             equations to be integrated.
!
!      t -- this is a double precision value of the independent
!           variable.
!
!      y(*) -- this double precision array contains the solution
!              components at t.
!
!      tout -- this is a double precision point at which a solution is
!              desired.
!
!      info(*) -- the basic task of the code is to integrate the
!             differential equations from t to tout and return an
!             answer at tout.  info(*) is an integer array which is used
!             to communicate exactly how you want this task to be
!             carried out.
!
!      rtol, atol -- these double precision quantities represent
!                    relative and absolute error tolerances which you
!                    provide to indicate how accurately you wish the
!                    solution to be computed.  you may choose them to be
!                    both scalars or else both vectors.
!
!      idid -- this scalar quantity is an indicator reporting what
!             the code did.  you must monitor this integer variable to
!             decide what action to take next.
!
!      rwork(*), lrw -- rwork(*) is a double precision work array of
!             length lrw which provides the code with needed storage
!             space.
!
!      iwork(*), liw -- iwork(*) is an integer work array of length liw
!             which provides the code with needed storage space and an
!             across call flag.
!
!  quantities which are used as input items are
!             neq, t, y(*), tout, info(*),
!             rtol, atol, rwork(1), lrw and liw.
!
!  quantities which may be altered by the code are
!             t, y(*), info(1), rtol, atol,
!             idid, rwork(*) and iwork(*).
!
! **********************************************************************
! * input -- what to do on the first call to ddeabm *
! **********************************************************************
!
!   the first call of the code is defined to be the start of each new
!   problem.  read through the descriptions of all the following items,
!   provide sufficient storage space for designated arrays, set
!   appropriate variables for the initialization of the problem, and
!   give information about how you want the problem to be solved.
!
!
!      df -- provide a subroutine of the form
!                               df(x,u,uprime)
!             to define the system of first order differential equations
!             which is to be solved.  for the given values of x and the
!             vector  u(*)=(u(1),u(2),...,u(neq)) , the subroutine must
!             evaluate the neq components of the system of differential
!             equations  du/dx=df(x,u)  and store the derivatives in the
!             array uprime(*), that is,  uprime(i) = * du(i)/dx *  for
!             equations i=1,...,neq.
!
!             subroutine df must not alter x or u(*).  you must declare
!             the name df in an external statement in your program that
!             calls ddeabm.  you must dimension u and uprime in df.
!
!      neq -- set it to the number of differential equations.
!             (neq >= 1)
!
!      t -- set it to the initial point of the integration.
!             you must use a program variable for t because the code
!             changes its value.
!
!      y(*) -- set this vector to the initial values of the neq solution
!             components at the initial point.  you must dimension y at
!             least neq in your calling program.
!
!      tout -- set it to the first point at which a solution
!             is desired.  you can take tout = t, in which case the code
!             will evaluate the derivative of the solution at t and
!             return. integration either forward in t  (tout > t)  or
!             backward in t  (tout < t)  is permitted.
!
!             the code advances the solution from t to tout using
!             step sizes which are automatically selected so as to
!             achieve the desired accuracy.  if you wish, the code will
!             return with the solution and its derivative following
!             each intermediate step (intermediate-output mode) so that
!             you can monitor them, but you still must provide tout in
!             accord with the basic aim of the code.
!
!             the first step taken by the code is a critical one
!             because it must reflect how fast the solution changes near
!             the initial point.  the code automatically selects an
!             initial step size which is practically always suitable for
!             the problem. by using the fact that the code will not step
!             past tout in the first step, you could, if necessary,
!             restrict the length of the initial step size.
!
!             for some problems it may not be permissible to integrate
!             past a point tstop because a discontinuity occurs there
!             or the solution or its derivative is not defined beyond
!             tstop.  when you have declared a tstop point (see info(4)
!             and rwork(1)), you have told the code not to integrate
!             past tstop.  in this case any tout beyond tstop is invalid
!             input.
!
!      info(*) -- use the info array to give the code more details about
!             how you want your problem solved.  this array should be
!             dimensioned of length 4.  you must respond to all of
!             the following items which are arranged as questions.  the
!             simplest use of the code corresponds to answering all
!             questions as yes ,i.e. setting all entries of info to 0.
!
!        info(1) -- this parameter enables the code to initialize
!               itself.  you must set it to indicate the start of every
!               new problem.
!
!            **** is this the first call for this problem ...
!                  yes -- set info(1) = 0
!                   no -- not applicable here.
!                         see below for continuation calls.  ****
!
!        info(2) -- how much accuracy you want of your solution
!               is specified by the error tolerances rtol and atol.
!               the simplest use is to take them both to be scalars.
!               to obtain more flexibility, they can both be vectors.
!               the code must be told your choice.
!
!            **** are both error tolerances rtol, atol scalars ...
!                  yes -- set info(2) = 0
!                         and input scalars for both rtol and atol
!                   no -- set info(2) = 1
!                         and input arrays for both rtol and atol ****
!
!        info(3) -- the code integrates from t in the direction
!               of tout by steps.  if you wish, it will return the
!               computed solution and derivative at the next
!               intermediate step (the intermediate-output mode) or
!               tout, whichever comes first.  this is a good way to
!               proceed if you want to see the behavior of the solution.
!               if you must have solutions at a great many specific
!               tout points, this code will compute them efficiently.
!
!            **** do you want the solution only at
!                 tout (and not at the next intermediate step) ...
!                  yes -- set info(3) = 0
!                   no -- set info(3) = 1 ****
!
!        info(4) -- to handle solutions at a great many specific
!               values tout efficiently, this code may integrate past
!               tout and interpolate to obtain the result at tout.
!               sometimes it is not possible to integrate beyond some
!               point tstop because the equation changes there or it is
!               not defined past tstop.  then you must tell the code
!               not to go past.
!
!            **** can the integration be carried out without any
!                 restrictions on the independent variable t ...
!                  yes -- set info(4)=0
!                   no -- set info(4)=1
!                         and define the stopping point tstop by
!                         setting rwork(1)=tstop ****
!
!      rtol, atol -- you must assign relative (rtol) and absolute (atol)
!             error tolerances to tell the code how accurately you want
!             the solution to be computed.  they must be defined as
!             program variables because the code may change them.  you
!             have two choices --
!                  both rtol and atol are scalars. (info(2)=0)
!                  both rtol and atol are vectors. (info(2)=1)
!             in either case all components must be non-negative.
!
!             the tolerances are used by the code in a local error test
!             at each step which requires roughly that
!                     abs(local error) <= rtol*abs(y)+atol
!             for each vector component.
!             (more specifically, a euclidean norm is used to measure
!             the size of vectors, and the error test uses the magnitude
!             of the solution at the beginning of the step.)
!
!             the true (global) error is the difference between the true
!             solution of the initial value problem and the computed
!             approximation.  practically all present day codes,
!             including this one, control the local error at each step
!             and do not even attempt to control the global error
!             directly.  roughly speaking, they produce a solution y(t)
!             which satisfies the differential equations with a
!             residual r(t),    dy(t)/dt = df(t,y(t)) + r(t)   ,
!             and, almost always, r(t) is bounded by the error
!             tolerances.  usually, but not always, the true accuracy of
!             the computed y is comparable to the error tolerances. this
!             code will usually, but not always, deliver a more accurate
!             solution if you reduce the tolerances and integrate again.
!             by comparing two such solutions you can get a fairly
!             reliable idea of the true error in the solution at the
!             bigger tolerances.
!
!             setting atol=0.d0 results in a pure relative error test on
!             that component. setting rtol=0. results in a pure absolute
!             error test on that component.  a mixed test with non-zero
!             rtol and atol corresponds roughly to a relative error
!             test when the solution component is much bigger than atol
!             and to an absolute error test when the solution component
!             is smaller than the threshold atol.
!
!             proper selection of the absolute error control parameters
!             atol  requires you to have some idea of the scale of the
!             solution components.  to acquire this information may mean
!             that you will have to solve the problem more than once. in
!             the absence of scale information, you should ask for some
!             relative accuracy in all the components (by setting  rtol
!             values non-zero) and perhaps impose extremely small
!             absolute error tolerances to protect against the danger of
!             a solution component becoming zero.
!
!             the code will not attempt to compute a solution at an
!             accuracy unreasonable for the machine being used.  it will
!             advise you if you ask for too much accuracy and inform
!             you as to the maximum accuracy it believes possible.
!
!      rwork(*) -- dimension this double precision work array of length
!             lrw in your calling program.
!
!      rwork(1) -- if you have set info(4)=0, you can ignore this
!             optional input parameter.  otherwise you must define a
!             stopping point tstop by setting   rwork(1) = tstop.
!             (for some problems it may not be permissible to integrate
!             past a point tstop because a discontinuity occurs there
!             or the solution or its derivative is not defined beyond
!             tstop.)
!
!      lrw -- set it to the declared length of the rwork array.
!             you must have  lrw >= 130+21*neq
!
!      iwork(*) -- dimension this integer work array of length liw in
!             your calling program.
!
!      liw -- set it to the declared length of the iwork array.
!             you must have  liw >= 51
!
! **********************************************************************
! * output -- after any return from ddeabm *
! **********************************************************************
!
!   the principal aim of the code is to return a computed solution at
!   tout, although it is also possible to obtain intermediate results
!   along the way.  to find out whether the code achieved its goal
!   or if the integration process was interrupted before the task was
!   completed, you must check the idid parameter.
!
!
!      t -- the solution was successfully advanced to the
!             output value of t.
!
!      y(*) -- contains the computed solution approximation at t.
!             you may also be interested in the approximate derivative
!             of the solution at t.  it is contained in
!             rwork(21),...,rwork(20+neq).
!
!      idid -- reports what the code did
!
!                         *** task completed ***
!                   reported by positive values of idid
!
!             idid = 1 -- a step was successfully taken in the
!                       intermediate-output mode.  the code has not
!                       yet reached tout.
!
!             idid = 2 -- the integration to tout was successfully
!                       completed (t=tout) by stepping exactly to tout.
!
!             idid = 3 -- the integration to tout was successfully
!                       completed (t=tout) by stepping past tout.
!                       y(*) is obtained by interpolation.
!
!                         *** task interrupted ***
!                   reported by negative values of idid
!
!             idid = -1 -- a large amount of work has been expended.
!                       (500 steps attempted)
!
!             idid = -2 -- the error tolerances are too stringent.
!
!             idid = -3 -- the local error test cannot be satisfied
!                       because you specified a zero component in atol
!                       and the corresponding computed solution
!                       component is zero.  thus, a pure relative error
!                       test is impossible for this component.
!
!             idid = -4 -- the problem appears to be stiff.
!
!             idid = -5,-6,-7,..,-32  -- not applicable for this code
!                       but used by other members of depac or possible
!                       future extensions.
!
!                         *** task terminated ***
!                   reported by the value of idid=-33
!
!             idid = -33 -- the code has encountered trouble from which
!                       it cannot recover.  a message is printed
!                       explaining the trouble and control is returned
!                       to the calling program. for example, this occurs
!                       when invalid input is detected.
!
!      rtol, atol -- these quantities remain unchanged except when
!             idid = -2.  in this case, the error tolerances have been
!             increased by the code to values which are estimated to be
!             appropriate for continuing the integration.  however, the
!             reported solution at t was obtained using the input values
!             of rtol and atol.
!
!      rwork, iwork -- contain information which is usually of no
!             interest to the user but necessary for subsequent calls.
!             however, you may find use for
!
!             rwork(11)--which contains the step size h to be
!                        attempted on the next step.
!
!             rwork(12)--if the tolerances have been increased by the
!                        code (idid = -2) , they were multiplied by the
!                        value in rwork(12).
!
!             rwork(13)--which contains the current value of the
!                        independent variable, i.e. the farthest point
!                        integration has reached. this will be different
!                        from t only when interpolation has been
!                        performed (idid=3).
!
!             rwork(20+i)--which contains the approximate derivative
!                        of the solution component y(i).  in ddeabm, it
!                        is obtained by calling subroutine df to
!                        evaluate the differential equation using t and
!                        y(*) when idid=1 or 2, and by interpolation
!                        when idid=3.
!
! **********************************************************************
! * input -- what to do to continue the integration *
! *             (calls after the first)             *
! **********************************************************************
!
!        this code is organized so that subsequent calls to continue the
!        integration involve little (if any) additional effort on your
!        part. you must monitor the idid parameter in order to determine
!        what to do next.
!
!        recalling that the principal task of the code is to integrate
!        from t to tout (the interval mode), usually all you will need
!        to do is specify a new tout upon reaching the current tout.
!
!        do not alter any quantity not specifically permitted below,
!        in particular do not alter neq, t, y(*), rwork(*), iwork(*) or
!        the differential equation in subroutine df. any such alteration
!        constitutes a new problem and must be treated as such, i.e.
!        you must start afresh.
!
!        you cannot change from vector to scalar error control or vice
!        versa (info(2)) but you can change the size of the entries of
!        rtol, atol.  increasing a tolerance makes the equation easier
!        to integrate.  decreasing a tolerance will make the equation
!        harder to integrate and should generally be avoided.
!
!        you can switch from the intermediate-output mode to the
!        interval mode (info(3)) or vice versa at any time.
!
!        if it has been necessary to prevent the integration from going
!        past a point tstop (info(4), rwork(1)), keep in mind that the
!        code will not integrate to any tout beyond the currently
!        specified tstop.  once tstop has been reached you must change
!        the value of tstop or set info(4)=0.  you may change info(4)
!        or tstop at any time but you must supply the value of tstop in
!        rwork(1) whenever you set info(4)=1.
!
!        the parameter info(1) is used by the code to indicate the
!        beginning of a new problem and to indicate whether integration
!        is to be continued.  you must input the value  info(1) = 0
!        when starting a new problem.  you must input the value
!        info(1) = 1  if you wish to continue after an interrupted task.
!        do not set  info(1) = 0  on a continuation call unless you
!        want the code to restart at the current t.
!
!                         *** following a completed task ***
!         if
!             idid = 1, call the code again to continue the integration
!                     another step in the direction of tout.
!
!             idid = 2 or 3, define a new tout and call the code again.
!                     tout must be different from t. you cannot change
!                     the direction of integration without restarting.
!
!                         *** following an interrupted task ***
!                     to show the code that you realize the task was
!                     interrupted and that you want to continue, you
!                     must take appropriate action and reset info(1) = 1
!         if
!             idid = -1, the code has attempted 500 steps.
!                     if you want to continue, set info(1) = 1 and
!                     call the code again. an additional 500 steps
!                     will be allowed.
!
!             idid = -2, the error tolerances rtol, atol have been
!                     increased to values the code estimates appropriate
!                     for continuing.  you may want to change them
!                     yourself.  if you are sure you want to continue
!                     with relaxed error tolerances, set info(1)=1 and
!                     call the code again.
!
!             idid = -3, a solution component is zero and you set the
!                     corresponding component of atol to zero.  if you
!                     are sure you want to continue, you must first
!                     alter the error criterion to use positive values
!                     for those components of atol corresponding to zero
!                     solution components, then set info(1)=1 and call
!                     the code again.
!
!             idid = -4, the problem appears to be stiff.  it is very
!                     inefficient to solve such problems with ddeabm.
!                     the code ddebdf in depac handles this task
!                     efficiently.  if you are absolutely sure you want
!                     to continue with ddeabm, set info(1)=1 and call
!                     the code again.
!
!                         *** following a terminated task ***
!         if
!             idid = -33, you cannot continue the solution of this
!                     problem.  an attempt to do so will result in your
!                     run being terminated.
!
! **********************************************************************
! *long description:
!
! ....   ddeabm is a variable order (one through twelve) adams code.
! ....   its complexity lies somewhere between that of dderkf and
! ....   ddebdf.  ddeabm is primarily designed to solve non-stiff and
! ....   mildly stiff differential equations when derivative evaluations
! ....   are expensive, high accuracy results are needed or answers at
! ....   many specific points are required. ddeabm attempts to discover
! ....   when it is not suitable for the task posed.
!
! *********************************************************************
!
!  AUTHORS
!
!	L. F. shampine 
!	H. A. Watts 
!	M. K. gordon
!	
!  SEE ALSO
!	l. f. shampine and h. a. watts, depac - design of a user
!		oriented package of ode solvers, report sand79-2374,
!		sandia laboratories, 1979.
!
!  HISTORY
!
!   820301  date written
!   890531  changed all specific intrinsics to generic.  (wrb)
!   890831  modified array declarations.  (wrb)
!   891006  cosmetic changes to prologue.  (wrb)
!   891024  changed references from dvnorm to dhvnrm.  (wrb)
!   891024  revision date from version 3.2
!   891214  prologue converted to version 4.0 format.  (bab)
!   900510  convert xerrwv calls to xermsg calls.  (rwc)
!   920501  reformatted the references section.  (wrb)
!	July, 2014 : Major refactoring into modern Fortran (jw)
!
!*****************************************************************************************    

	implicit none

	class(ddeabm_class),intent(inout)		:: me	
	integer,intent(in)						:: neq
	real(wp),intent(inout)					:: t
	real(wp),dimension(neq),intent(inout)	:: y
	real(wp),intent(in)						:: tout
	integer,dimension(4),intent(inout)		:: info
	real(wp),dimension(:),intent(inout)		:: rtol
	real(wp),dimension(:),intent(inout)		:: atol
	integer,intent(out)						:: idid
	real(wp),dimension(lrw),intent(inout)	:: rwork
	integer,intent(in)						:: lrw
	integer,dimension(liw),intent(inout)	:: iwork
	integer,intent(in)						:: liw

	integer ialpha, ibeta, idelsn, ifouru, ig, ihold,&
		ip, iphi, ipsi, isig, itold, itstar, itwou,&
		iv, iw, iwt, iyp, iypout, iyy, igi,ixold
	logical start,phase1,nornd,stiff,intout
	character(len=8) :: xern1
	character(len=16) :: xern3
	      
!
!     check for an apparent infinite loop
!
!***first executable statement  ddeabm
      if ( info(1) == 0 ) iwork(liw) = 0
      if (iwork(liw) >= 5) then
         if (t == rwork(21 + neq)) then
            write (xern3, '(1pe15.6)') t
            call xermsg ('slatec', 'ddeabm',&
               'an apparent infinite loop has been detected.$$' //&
               'you have made repeated calls at t = ' // xern3 //&
               ' and the integration has not advanced.  check the ' //&
               'way you have set parameters for the call to the ' //&
               'code, particularly info(1).', 13, 2)
            return
         end if
      end if
!
!     check lrw and liw for sufficient storage allocation
!
      idid=0
      if (lrw < 130+21*neq) then
         write (xern1, '(i8)') lrw
         call xermsg ('slatec', 'ddeabm', 'the length of the rwork ' //&
            'array must be at least 130 + 21*neq.$$' //&
            'you have called the code with lrw = ' // xern1, 1, 1)
         idid=-33
      end if

      if (liw < 51) then
         write (xern1, '(i8)') liw
         call xermsg ('slatec', 'ddeabm', 'the length of the iwork ' //&
            'array must be at least 51.$$you have called the code ' //&
            'with liw = ' // xern1, 2, 1)
         idid=-33
      end if
      
	!make sure the deriv function was set:
		if (.not. associated(me%df)) then
         call xermsg ('slatec', 'ddeabm', 'the derivative function DF '//&
         				' has not been associated.' // xern1, 0, 0)		
         idid=-33
		end if
    
!
!     compute the indices for the arrays to be stored in the work array
!
      iypout = 21
      itstar = neq + 21
      iyp = 1 + itstar
      iyy = neq + iyp
      iwt = neq + iyy
      ip = neq + iwt
      iphi = neq + ip
      ialpha = (neq*16) + iphi
      ibeta = 12 + ialpha
      ipsi = 12 + ibeta
      iv = 12 + ipsi
      iw = 12 + iv
      isig = 12 + iw
      ig = 13 + isig
      igi = 13 + ig
      ixold = 11 + igi
      ihold = 1 + ixold
      itold = 1 + ihold
      idelsn = 1 + itold
      itwou = 1 + idelsn
      ifouru = 1 + itwou

      rwork(itstar) = t
      
      if (info(1) /= 0) then      
		  start = iwork(21) /= (-1)
		  phase1 = iwork(22) /= (-1)
		  nornd = iwork(23) /= (-1)
		  stiff = iwork(24) /= (-1)
		  intout = iwork(25) /= (-1)
      end if

      call me%ddes(neq,t,y,tout,info,rtol,atol,idid,rwork(iypout),&
               rwork(iyp),rwork(iyy),rwork(iwt),rwork(ip),rwork(iphi),&
               rwork(ialpha),rwork(ibeta),rwork(ipsi),rwork(iv),&
               rwork(iw),rwork(isig),rwork(ig),rwork(igi),rwork(11),&
               rwork(12),rwork(13),rwork(ixold),rwork(ihold),&
               rwork(itold),rwork(idelsn),rwork(1),rwork(itwou),&
               rwork(ifouru),start,phase1,nornd,stiff,intout,iwork(26),&
               iwork(27),iwork(28),iwork(29),iwork(30),iwork(31),&
               iwork(32),iwork(33),iwork(34),iwork(35),iwork(45))

      iwork(21) = -1
      if (start) iwork(21) = 1
      iwork(22) = -1
      if (phase1) iwork(22) = 1
      iwork(23) = -1
      if (nornd) iwork(23) = 1
      iwork(24) = -1
      if (stiff) iwork(24) = 1
      iwork(25) = -1
      if (intout) iwork(25) = 1

      if (idid /= (-2)) iwork(liw) = iwork(liw) + 1
      if (t /= rwork(itstar)) iwork(liw) = 0

!*****************************************************************************************    
	end subroutine ddeabm
!*****************************************************************************************    
      
!*****************************************************************************************    
	subroutine ddes (me, neq, t, y, tout, info, rtol, atol, idid,&
         ypout, yp, yy, wt, p, phi, alpha, beta, psi, v, w, sig, g, gi,&
         h, eps, x, xold, hold, told, delsgn, tstop, twou, fouru, start,&
         phase1, nornd, stiff, intout, ns, kord, kold, init, ksteps,&
         kle4, iquit, kprev, ivc, iv, kgi)
!*****************************************************************************************    
!****f* ddeabm_module/ddeabm
!
!  NAME
!    ddeabm
!
!  DESCRIPTION
!   ddeabm merely allocates storage for ddes to relieve the user of the
!   inconvenience of a long call list.  consequently ddes is used as
!   described in the comments for ddeabm.
!
!  AUTHOR
!	watts, h. a., (snla)
!
!  SEE ALSO
!	ddeabm
!
!  HISTORY
!   820301  date written
!   890531  changed all specific intrinsics to generic.  (wrb)
!   890831  modified array declarations.  (wrb)
!   891214  prologue converted to version 4.0 format.  (bab)
!   900328  added type section.  (wrb)
!   900510  convert xerrwv calls to xermsg calls, cvt gotos to
!           if-then-else.  (rwc)
!   910722  updated author section.  (als)
!
!*****************************************************************************************    

	implicit none
	
	class(ddeabm_class),intent(inout)			:: me
	integer,intent(in)							:: neq
	real(wp),intent(inout)						:: t
	real(wp),dimension(neq),intent(inout)		:: y
	real(wp),intent(in)							:: tout
	integer,dimension(4),intent(inout)			:: info
	real(wp),dimension(:),intent(inout)			:: rtol
	real(wp),dimension(:),intent(inout)			:: atol
	integer,intent(inout)						:: idid
	real(wp),dimension(neq),intent(inout)		:: ypout
	real(wp),dimension(neq),intent(inout)		:: yp
	real(wp),dimension(neq),intent(inout)		:: yy
	real(wp),dimension(neq),intent(inout)		:: wt
	real(wp),dimension(neq),intent(inout)		:: p
	real(wp),dimension(neq,16),intent(inout)	:: phi
	real(wp),dimension(12),intent(inout)		:: alpha
	real(wp),dimension(12),intent(inout)		:: beta
	real(wp),dimension(12),intent(inout)		:: psi
	real(wp),dimension(12),intent(inout)		:: v
	real(wp),dimension(12),intent(inout)		:: w
	real(wp),dimension(13),intent(inout)		:: sig
	real(wp),dimension(13),intent(inout)		:: g
	real(wp),dimension(11),intent(inout)		:: gi
	real(wp),intent(inout)						:: h
	real(wp),intent(inout)						:: eps
	real(wp),intent(inout)						:: x
	real(wp),intent(inout)						:: xold
	real(wp),intent(inout)						:: hold
	real(wp),intent(inout)						:: told
	real(wp),intent(inout)						:: delsgn
	real(wp),intent(inout)						:: tstop
	real(wp),intent(inout)						:: twou
	real(wp),intent(inout)						:: fouru
	logical,intent(inout)						:: start
	logical,intent(inout)						:: phase1
	logical,intent(inout)						:: nornd
	logical,intent(inout)						:: stiff
	logical,intent(inout)						:: intout
	integer,intent(inout)						:: ns
	integer,intent(inout)						:: kord
	integer,intent(inout)						:: kold
	integer,intent(inout)						:: init
	integer,intent(inout)						:: ksteps
	integer,intent(inout)						:: kle4
	integer,intent(inout)						:: iquit
	integer,intent(inout)						:: kprev
	integer,intent(inout)						:: ivc
	integer,dimension(10),intent(inout)			:: iv
	integer,intent(inout)						:: kgi
	
	integer k, l, ltol, natolp, nrtolp
	real(wp) a, absdel, del, dt, ha, u
	logical :: crash
	character(len=8) :: xern1
	character(len=16) :: xern3, xern4
	
!.......................................................................
!
! on the first call , perform initialization

	if (info(1) == 0) then
		u=d1mach4		! machine unit roundoff quantity
		twou=2.d0*u   	! set associated machine dependent parameters
		fouru=4.d0*u	!
		iquit=0			! set termination flag
		init=0			! set initialization indicator
		ksteps=0		! set counter for attempted steps
		intout= .false.	! set indicator for intermediate-output
		stiff= .false.	! set indicator for stiffness detection
		kle4=0			! set step counter for stiffness detection
		start= .true. 	! set indicators for steps code
		phase1= .true. 	!
		nornd= .true. 	!
		info(1)=1   	! reset info(1) for subsequent calls
	end if

!.......................................................................
!
!      check validity of input parameters on each entry
!
      if (info(1) /= 0  .and.  info(1) /= 1) then
         write (xern1, '(i8)') info(1)
         call xermsg ('slatec', 'ddes', 'in ddeabm, info(1) must be ' //&
            'set to 0 for the start of a new problem, and must be ' //&
            'set to 1 following an interrupted task.  you are ' //&
            'attempting to continue the integration illegally by ' //&
            'calling the code with info(1) = ' // xern1, 3, 1)
         idid=-33
      end if

      if (info(2) /= 0  .and.  info(2) /= 1) then
         write (xern1, '(i8)') info(2)
         call xermsg ('slatec', 'ddes', 'in ddeabm, info(2) must be ' //&
            '0 or 1 indicating scalar and vector error tolerances, ' //&
            'respectively.  you have called the code with info(2) = ' //&
            xern1, 4, 1)
         idid=-33
      end if

      if (info(3) /= 0  .and.  info(3) /= 1) then
         write (xern1, '(i8)') info(3)
         call xermsg ('slatec', 'ddes', 'in ddeabm, info(3) must be ' //&
            '0 or 1 indicating the interval or intermediate-output ' //&
            'mode of integration, respectively.  you have called ' //&
            'the code with  info(3) = ' // xern1, 5, 1)
         idid=-33
      end if

      if (info(4) /= 0  .and.  info(4) /= 1) then
         write (xern1, '(i8)') info(4)
         call xermsg ('slatec', 'ddes', 'in ddeabm, info(4) must be ' //&
            '0 or 1 indicating whether or not the integration ' //&
            'interval is to be restricted by a point tstop.  you ' //&
            'have called the code with info(4) = ' // xern1, 14, 1)
         idid=-33
      end if

      if (neq < 1) then
         write (xern1, '(i8)') neq
         call xermsg ('slatec', 'ddes', 'in ddeabm,  the number of ' //&
            'equations neq must be a positive integer.  you have ' //&
            'called the code with  neq = ' // xern1, 6, 1)
         idid=-33
      end if

      nrtolp = 0
      natolp = 0
      do 90 k=1,neq
         if (nrtolp == 0 .and. rtol(k) < 0.d0) then
            write (xern1, '(i8)') k
            write (xern3, '(1pe15.6)') rtol(k)
            call xermsg ('slatec', 'ddes', 'in ddeabm, the relative ' //&
               'error tolerances rtol must be non-negative.  you ' //&
               'have called the code with  rtol(' // xern1 // ') = ' //&
               xern3 // '.  in the case of vector error tolerances, ' //&
               'no further checking of rtol components is done.', 7, 1)
            idid = -33
            nrtolp = 1
         end if

         if (natolp == 0 .and. atol(k) < 0.d0) then
            write (xern1, '(i8)') k
            write (xern3, '(1pe15.6)') atol(k)
            call xermsg ('slatec', 'ddes', 'in ddeabm, the absolute ' //&
               'error tolerances atol must be non-negative.  you ' //&
               'have called the code with  atol(' // xern1 // ') = ' //&
               xern3 // '.  in the case of vector error tolerances, ' //&
               'no further checking of atol components is done.', 8, 1)
            idid = -33
            natolp = 1
         end if

         if (info(2) == 0) go to 100
         if (natolp>0 .and. nrtolp>0) go to 100
   90 continue

  100 if (info(4) == 1) then
         if (sign(1.d0,tout-t) /= sign(1.d0,tstop-t)&
            .or. abs(tout-t) > abs(tstop-t)) then
            write (xern3, '(1pe15.6)') tout
            write (xern4, '(1pe15.6)') tstop
            call xermsg ('slatec', 'ddes', 'in ddeabm, you have ' //&
               'called the code with  tout = ' // xern3 // ' but ' //&
               'you have also told the code (info(4) = 1) not to ' //&
               'integrate past the point tstop = ' // xern4 //&
               ' these instructions conflict.', 14, 1)
            idid=-33
         end if
      end if

!     check some continuation possibilities

      if (init /= 0) then
         if (t == tout) then
            write (xern3, '(1pe15.6)') t
            call xermsg ('slatec', 'ddes', 'in ddeabm, you have ' //&
               'called the code with  t = tout = ' // xern3 //&
               '$$this is not allowed on continuation calls.', 9, 1)
            idid=-33
         end if

         if (t /= told) then
            write (xern3, '(1pe15.6)') told
            write (xern4, '(1pe15.6)') t
            call xermsg ('slatec', 'ddes', 'in ddeabm, you have ' //&
               'changed the value of t from ' // xern3 // ' to ' //&
               xern4 //'  this is not allowed on continuation calls.',&
               10, 1)
            idid=-33
         end if

         if (init /= 1) then
            if (delsgn*(tout-t) < 0.d0) then
               write (xern3, '(1pe15.6)') tout
               call xermsg ('slatec', 'ddes', 'in ddeabm, by ' //&
                  'calling the code with tout = ' // xern3 //&
                  ' you are attempting to change the direction of ' //&
                  'integration.$$this is not allowed without ' //&
                  'restarting.', 11, 1)
               idid=-33
            end if
         end if
      end if

!     invalid input detected

      if (idid == -33) then
         if (iquit /= -33) then
            iquit = -33
            info(1) = -1
         else
            call xermsg ('slatec', 'ddes', 'in ddeabm, invalid ' //&
               'input was detected on successive entries.  it is ' //&
               'impossible to proceed because you have not ' //&
               'corrected the problem, so execution is being ' //&
               'terminated.', 12, 2)
         end if
         return
      end if

!.......................................................................
!
!     rtol = atol = 0. is allowed as valid input and interpreted as
!     asking for the most accurate solution possible. in this case,
!     the relative error tolerance rtol is reset to the smallest value
!     fouru which is likely to be reasonable for this method and machine
!
      do k=1,neq
        if (rtol(k)+atol(k) <= 0.d0) then
        	rtol(k)=fouru
        	idid=-2
        end if
        if (info(2) == 0) exit
      end do

      if (idid == -2) then
		! rtol=atol=0 on input, so rtol is changed to a small positive value
      	info(1)=-1
      	return    
      end if

!     branch on status of initialization indicator
!            init=0 means initial derivatives and nominal step size
!                   and direction not yet set
!            init=1 means nominal step size and direction not yet set
!            init=2 means no further initialization required
!
      if (init == 0) go to 210
      if (init == 1) go to 220
      go to 240
!
!.......................................................................
!
!     more initialization --
!                         -- evaluate initial derivatives
!
  210 init=1
      a=t
      call me%df(a,y,yp)
      if (t /= tout) go to 220
      idid=2
      do l = 1,neq
         ypout(l) = yp(l)
      end do
      told=t
      return
!
!                         -- set independent and dependent variables
!                                              x and yy(*) for steps
!                         -- set sign of integration direction
!                         -- initialize the step size
!
  220 init = 2
      x = t
      do l = 1,neq
        yy(l) = y(l)
      end do
      delsgn = sign(1.0d0,tout-t)
      h = sign(max(fouru*abs(x),abs(tout-x)),tout-x)
!
!.......................................................................
!
!   on each call set information which determines the allowed interval
!   of integration before returning with an answer at tout
!
  240 del = tout - t
      absdel = abs(del)
!
!.......................................................................
!
!   if already past output point, interpolate and return
!
  250 if (abs(x-t) >= absdel) then
		  call dintp(x,yy,tout,y,ypout,neq,kold,phi,ivc,iv,kgi,gi,alpha,g,w,xold,p)
		  idid = 3
		  if (x == tout) then
			 idid = 2
			 intout = .false.
		  end if
		  t = tout
		  told = t
		  return	  
      end if
      
!   if cannot go past tstop and sufficiently close, extrapolate and return

      if (info(4)==1 .and. abs(tstop-x)<fouru*abs(x)) then
		  dt = tout - x
		  do l = 1,neq
			y(l) = yy(l) + dt*yp(l)
		  end do
		  call me%df(tout,y,ypout)
		  idid = 3
		  t = tout
		  told = t
		  return
      end if
      
    if (info(3) == 0  .or.  .not.intout) go to 300
!
!   intermediate-output mode
!
      idid = 1
      do l = 1,neq
        y(l)=yy(l)
        ypout(l) = yp(l)
      end do
      t = x
      told = t
      intout = .false.
      return
!
!.......................................................................
!
!     monitor number of steps attempted
!
  300 if (ksteps > me%maxnum) then
  
	! a significant amount of work has been expended
      idid=-1
      ksteps=0
      if (stiff) then
		! problem appears to be stiff
      	idid=-4
      	stiff= .false.
      	kle4=0
      end if

      do l = 1,neq
        y(l) = yy(l)
        ypout(l) = yp(l)
      end do
      t = x
      told = t
      info(1) = -1
      intout = .false.
      return
      
    end if
!
!.......................................................................
!
!   limit step size, set weight vector and take a step
!
  330 ha = abs(h)
      if (info(4) /= 1) go to 340
      ha = min(ha,abs(tstop-x))
  340 h = sign(ha,h)
      eps = 1.0d0
      ltol = 1
      do 350 l = 1,neq
        if (info(2) == 1) ltol = l
        wt(l) = rtol(ltol)*abs(yy(l)) + atol(ltol)
        if (wt(l) <= 0.0d0) go to 360
  350   continue
      go to 380
!
!                       relative error criterion inappropriate
  360 idid = -3
      do 370 l = 1,neq
        y(l) = yy(l)
  370   ypout(l) = yp(l)
      t = x
      told = t
      info(1) = -1
      intout = .false.
      return
!
  380 call me%dsteps(neq,yy,x,h,eps,wt,start,hold,kord,kold,crash,phi,p,&
                 yp,psi,alpha,beta,sig,v,w,g,phase1,ns,nornd,ksteps,&
                 twou,fouru,xold,kprev,ivc,iv,kgi,gi)
!
!.......................................................................
!
      if (.not.crash) go to 420
!
!                       tolerances too small
      idid = -2
      rtol(1) = eps*rtol(1)
      atol(1) = eps*atol(1)
      if (info(2) == 0) go to 400
      do 390 l = 2,neq
        rtol(l) = eps*rtol(l)
  390   atol(l) = eps*atol(l)
  400 do 410 l = 1,neq
        y(l) = yy(l)
  410   ypout(l) = yp(l)
      t = x
      told = t
      info(1) = -1
      intout = .false.
      return
!
!   (stiffness test) count number of consecutive steps taken with the
!   order of the method being less or equal to four
!
  420 kle4 = kle4 + 1
      if (kold > 4) kle4 = 0
      if (kle4 >= 50) stiff = .true.
      intout = .true.
      go to 250
      
!*****************************************************************************************    
	end subroutine ddes
!*****************************************************************************************    
      
!*****************************************************************************************    
	subroutine dhstrt (me,neq,a,b,y,yprime,etol,morder,small,big,spy,pv,yp,sf,h)
!*****************************************************************************************    
!***begin prologue  dhstrt
!***subsidiary
!***purpose  subsidiary to ddeabm, ddebdf and dderkf
!***library   slatec
!***type      double precision (hstart-s, dhstrt-d)
!***author  watts, h. a., (snla)
!***description
!
!   dhstrt computes a starting step size to be used in solving initial
!   value problems in ordinary differential equations.
!
! **********************************************************************
!  abstract
!
!     subroutine dhstrt computes a starting step size to be used by an
!     initial value method in solving ordinary differential equations.
!     it is based on an estimate of the local lipschitz constant for the
!     differential equation   (lower bound on a norm of the jacobian) ,
!     a bound on the differential equation  (first derivative) , and
!     a bound on the partial derivative of the equation with respect to
!     the independent variable.
!     (all approximated near the initial point a)
!
!     subroutine dhstrt uses a function subprogram dhvnrm for computing
!     a vector norm. the maximum norm is presently utilized though it
!     can easily be replaced by any other vector norm. it is presumed
!     that any replacement norm routine would be carefully coded to
!     prevent unnecessary underflows or overflows from occurring, and
!     also, would not alter the vector or number of components.
!
! **********************************************************************
!  on input you must provide the following
!
!      df -- this is a subroutine of the form
!                               df(x,u,uprime)
!             which defines the system of first order differential
!             equations to be solved. for the given values of x and the
!             vector  u(*)=(u(1),u(2),...,u(neq)) , the subroutine must
!             evaluate the neq components of the system of differential
!             equations  du/dx=df(x,u)  and store the derivatives in the
!             array uprime(*), that is,  uprime(i) = * du(i)/dx *  for
!             equations i=1,...,neq.
!
!             subroutine df must not alter x or u(*). you must declare
!             the name df in an external statement in your program that
!             calls dhstrt. you must dimension u and uprime in df.
!
!      neq -- this is the number of (first order) differential equations
!             to be integrated.
!
!      a -- this is the initial point of integration.
!
!      b -- this is a value of the independent variable used to define
!             the direction of integration. a reasonable choice is to
!             set  b  to the first point at which a solution is desired.
!             you can also use  b, if necessary, to restrict the length
!             of the first integration step because the algorithm will
!             not compute a starting step length which is bigger than
!             abs(b-a), unless  b  has been chosen too close to  a.
!             (it is presumed that dhstrt has been called with  b
!             different from  a  on the machine being used. also see the
!             discussion about the parameter  small.)
!
!      y(*) -- this is the vector of initial values of the neq solution
!             components at the initial point  a.
!
!      yprime(*) -- this is the vector of derivatives of the neq
!             solution components at the initial point  a.
!             (defined by the differential equations in subroutine df)
!
!      etol -- this is the vector of error tolerances corresponding to
!             the neq solution components. it is assumed that all
!             elements are positive. following the first integration
!             step, the tolerances are expected to be used by the
!             integrator in an error test which roughly requires that
!                        abs(local error)  <=  etol
!             for each vector component.
!
!      morder -- this is the order of the formula which will be used by
!             the initial value method for taking the first integration
!             step.
!
!      small -- this is a small positive machine dependent constant
!             which is used for protecting against computations with
!             numbers which are too small relative to the precision of
!             floating point arithmetic.  small  should be set to
!             (approximately) the smallest positive double precision
!             number such that  (1.+small) > 1.  on the machine being
!             used. the quantity  small**(3/8)  is used in computing
!             increments of variables for approximating derivatives by
!             differences.  also the algorithm will not compute a
!             starting step length which is smaller than
!             100*small*abs(a).
!
!      big -- this is a large positive machine dependent constant which
!             is used for preventing machine overflows. a reasonable
!             choice is to set big to (approximately) the square root of
!             the largest double precision number which can be held in
!             the machine.
!
!      spy(*),pv(*),yp(*),sf(*) -- these are double precision work
!             arrays of length neq which provide the routine with needed
!             storage space.
!
! **********************************************************************
!  on output  (after the return from dhstrt),
!
!      h -- is an appropriate starting step size to be attempted by the
!             differential equation method.
!
!           all parameters in the call list remain unchanged except for
!           the working arrays spy(*),pv(*),yp(*), and sf(*).
!
! **********************************************************************
!
!***see also  ddeabm, ddebdf, dderkf
!***routines called  dhvnrm
!***revision history  (yymmdd)
!   820301  date written
!   890531  changed all specific intrinsics to generic.  (wrb)
!   890831  modified array declarations.  (wrb)
!   890911  removed unnecessary intrinsics.  (wrb)
!   891024  changed references from dvnorm to dhvnrm.  (wrb)
!   891214  prologue converted to version 4.0 format.  (bab)
!   900328  added type section.  (wrb)
!   910722  updated author section.  (als)
! 
!*****************************************************************************************    

	class(ddeabm_class),intent(inout) :: me	
	integer,intent(in) :: neq
	real(wp),intent(in) :: a
	real(wp),intent(in) :: b	
    real(wp),intent(in),dimension(neq) :: y
    real(wp),intent(in),dimension(neq) :: yprime
    real(wp),intent(in),dimension(neq) :: etol
    integer,intent(in) :: morder
    real(wp),intent(in) :: small
    real(wp),intent(in) :: big
    real(wp),intent(out) :: h
    real(wp),dimension(neq),intent(inout) :: spy
    real(wp),dimension(neq),intent(inout) :: pv
    real(wp),dimension(neq),intent(inout) :: yp
    real(wp),dimension(neq),intent(inout) :: sf
	
	integer j, k, lk
	real(wp) absdx, da, delf, dely,&
		dfdub, dfdxb,&
		dx, dy, fbnd, relper,&
		srydpb, tolexp, tolmin, tolp, tolsum, ydpb
	
!     ..................................................................
!
!     begin block permitting ...exits to 160
!***first executable statement  dhstrt
         dx = b - a
         absdx = abs(dx)
         relper = small**0.375d0
!
!        ...............................................................
!
!             compute an approximate bound (dfdxb) on the partial
!             derivative of the equation with respect to the
!             independent variable. protect against an overflow.
!             also compute a bound (fbnd) on the first derivative
!             locally.
!
         da = sign(max(min(relper*abs(a),absdx),&
                          100.0d0*small*abs(a)),dx)
         if (da == 0.0d0) da = relper*dx
         call me%df(a+da,y,sf)
         do 10 j = 1, neq
            yp(j) = sf(j) - yprime(j)
   10    continue
         delf = dhvnrm(yp,neq)
         dfdxb = big
         if (delf < big*abs(da)) dfdxb = delf/abs(da)
         fbnd = dhvnrm(sf,neq)
!
!        ...............................................................
!
!             compute an estimate (dfdub) of the local lipschitz
!             constant for the system of differential equations. this
!             also represents an estimate of the norm of the jacobian
!             locally.  three iterations (two when neq=1) are used to
!             estimate the lipschitz constant by numerical differences.
!             the first perturbation vector is based on the initial
!             derivatives and direction of integration. the second
!             perturbation vector is formed using another evaluation of
!             the differential equation.  the third perturbation vector
!             is formed using perturbations based only on the initial
!             values. components that are zero are always changed to
!             non-zero values (except on the first iteration). when
!             information is available, care is taken to ensure that
!             components of the perturbation vector have signs which are
!             consistent with the slopes of local solution curves.
!             also choose the largest bound (fbnd) for the first
!             derivative.
!
!                               perturbation vector size is held
!                               constant for all iterations. compute
!                               this change from the
!                                       size of the vector of initial
!                                       values.
         dely = relper*dhvnrm(y,neq)
         if (dely == 0.0d0) dely = relper
         dely = sign(dely,dx)
         delf = dhvnrm(yprime,neq)
         fbnd = max(fbnd,delf)
         if (delf == 0.0d0) go to 30
!           use initial derivatives for first perturbation
            do 20 j = 1, neq
               spy(j) = yprime(j)
               yp(j) = yprime(j)
   20       continue
         go to 50
   30    continue
!           cannot have a null perturbation vector
            do 40 j = 1, neq
               spy(j) = 0.0d0
               yp(j) = 1.0d0
   40       continue
            delf = dhvnrm(yp,neq)
   50    continue
!
         dfdub = 0.0d0
         lk = min(neq+1,3)
         do 140 k = 1, lk
!           define perturbed vector of initial values
            do 60 j = 1, neq
               pv(j) = y(j) + dely*(yp(j)/delf)
   60       continue
            if (k == 2) go to 80
!              evaluate derivatives associated with perturbed
!              vector  and  compute corresponding differences
               call me%df(a,pv,yp)
               do 70 j = 1, neq
                  pv(j) = yp(j) - yprime(j)
   70          continue
            go to 100
   80       continue
!              use a shifted value of the independent variable
!                                    in computing one estimate
               call me%df(a+da,pv,yp)
               do 90 j = 1, neq
                  pv(j) = yp(j) - sf(j)
   90          continue
  100       continue
!           choose largest bounds on the first derivative
!                          and a local lipschitz constant
            fbnd = max(fbnd,dhvnrm(yp,neq))
            delf = dhvnrm(pv,neq)
!        ...exit
            if (delf >= big*abs(dely)) go to 150
            dfdub = max(dfdub,delf/abs(dely))
!     ......exit
            if (k == lk) go to 160
!           choose next perturbation vector
            if (delf == 0.0d0) delf = 1.0d0
            do 130 j = 1, neq
               if (k == 2) go to 110
                  dy = abs(pv(j))
                  if (dy == 0.0d0) dy = delf
               go to 120
  110          continue
                  dy = y(j)
                  if (dy == 0.0d0) dy = dely/relper
  120          continue
               if (spy(j) == 0.0d0) spy(j) = yp(j)
               if (spy(j) /= 0.0d0) dy = sign(dy,spy(j))
               yp(j) = dy
  130       continue
            delf = dhvnrm(yp,neq)
  140    continue
  150    continue
!
!        protect against an overflow
         dfdub = big
  160 continue
!
!     ..................................................................
!
!          compute a bound (ydpb) on the norm of the second derivative
!
      ydpb = dfdxb + dfdub*fbnd
!
!     ..................................................................
!
!          define the tolerance parameter upon which the starting step
!          size is to be based.  a value in the middle of the error
!          tolerance range is selected.
!
      tolmin = big
      tolsum = 0.0d0
      do 170 k = 1, neq
         tolexp = log10(etol(k))
         tolmin = min(tolmin,tolexp)
         tolsum = tolsum + tolexp
  170 continue
      tolp = 10.0d0**(0.5d0*(tolsum/neq + tolmin)/(morder+1))
!
!     ..................................................................
!
!          compute a starting step size based on the above first and
!          second derivative information
!
!                            restrict the step length to be not bigger
!                            than abs(b-a).   (unless  b  is too close
!                            to  a)
      h = absdx
!
      if (ydpb /= 0.0d0 .or. fbnd /= 0.0d0) go to 180
!
!        both first derivative term (fbnd) and second
!                     derivative term (ydpb) are zero
         if (tolp < 1.0d0) h = absdx*tolp
      go to 200
  180 continue
!
      if (ydpb /= 0.0d0) go to 190
!
!        only second derivative term (ydpb) is zero
         if (tolp < fbnd*absdx) h = tolp/fbnd
      go to 200
  190 continue
!
!        second derivative term (ydpb) is non-zero
         srydpb = sqrt(0.5d0*ydpb)
         if (tolp < srydpb*absdx) h = tolp/srydpb
  200 continue
!
!     further restrict the step length to be not
!                               bigger than  1/dfdub
      if (h*dfdub > 1.0d0) h = 1.0d0/dfdub
!
!     finally, restrict the step length to be not
!     smaller than  100*small*abs(a).  however, if
!     a=0. and the computed h underflowed to zero,
!     the algorithm returns  small*abs(b)  for the
!                                     step length.
      h = max(h,100.0d0*small*abs(a))
      if (h == 0.0d0) h = small*abs(b)
!
!     now set direction of integration
      h = sign(h,dx)

!*****************************************************************************************    
	end subroutine dhstrt
!*****************************************************************************************    
      
!*****************************************************************************************    
	function dhvnrm (v, n) result(m)
!*****************************************************************************************
!
!  NAME
!	dhvnrm
!	
!  DESCRIPTION
!	Compute the maximum norm of the vector v of length n
!
!  HISTORY
!	JW : 7/1/2014 : replace original routine
!
!*****************************************************************************************    

	implicit none
	
	real(wp)		 					:: m
	integer,intent(in) 					:: n
	real(wp),dimension(:),intent(in) 	:: v
   
	m = maxval(abs(v(1:n)))

!*****************************************************************************************    
	end function dhvnrm
!*****************************************************************************************    
      
!*****************************************************************************************    
	subroutine dintp (x,y,xout,yout,ypout,neqn,kold,phi,ivc,iv,kgi,gi,alpha,og,ow,ox,oy)
!*****************************************************************************************    
!***purpose  approximate the solution at xout by evaluating the
!            polynomial computed in dsteps at xout.  must be used in
!            conjunction with dsteps.
!***library   slatec (depac)
!***category  i1a1b
!***type      double precision (sintrp-s, dintp-d)
!***keywords  adams method, depac, initial value problems, ode,
!             ordinary differential equations, predictor-corrector,
!             smooth interpolant
!***author  watts, h. a., (snla)
!***description
!
!   the methods in subroutine  dsteps  approximate the solution near  x
!   by a polynomial.  subroutine  dintp  approximates the solution at
!   xout  by evaluating the polynomial there.  information defining this
!   polynomial is passed from  dsteps  so  dintp  cannot be used alone.
!
!   subroutine dsteps is completely explained and documented in the text
!   "computer solution of ordinary differential equations, the initial
!   value problem"  by l. f. shampine and m. k. gordon.
!
!   input to dintp --
!
!   the user provides storage in the calling program for the arrays in
!   the call list
!      dimension y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),oy(neqn)
!                and alpha(12),og(13),ow(12),gi(11),iv(10)
!   and defines
!      xout -- point at which solution is desired.
!   the remaining parameters are defined in  dsteps  and passed to
!   dintp  from that subroutine
!
!   output from  dintp --
!
!      yout(*) -- solution at  xout
!      ypout(*) -- derivative of solution at  xout
!   the remaining parameters are returned unaltered from their input
!   values.  integration with  dsteps  may be continued.
!
!***references  h. a. watts, a smoother interpolant for de/step, intrp
!                 ii, report sand84-0293, sandia laboratories, 1984.
!***routines called  (none)
!***revision history  (yymmdd)
!   840201  date written
!   890831  modified array declarations.  (wrb)
!   890831  revision date from version 3.2
!   891214  prologue converted to version 4.0 format.  (bab)
!   920501  reformatted the references section.  (wrb)
!
!*****************************************************************************************    

	implicit none
	
	integer,intent(in)						:: neqn
	real(wp),intent(in)						:: x
	real(wp),dimension(neqn),intent(in)		:: y
	real(wp),intent(in)						:: xout
	integer,intent(in)						:: kold
	real(wp),dimension(neqn,16),intent(in)	:: phi
	integer,intent(in)						:: ivc
	integer,dimension(10),intent(in)		:: iv
	integer,intent(in)						:: kgi
	real(wp),dimension(11),intent(in)		:: gi
	real(wp),dimension(12),intent(in)		:: alpha
	real(wp),dimension(13),intent(in)		:: og
	real(wp),dimension(12),intent(in)		:: ow
	real(wp),dimension(neqn),intent(in)		:: oy
	real(wp),intent(in)						:: ox	
	real(wp),dimension(neqn),intent(out)	:: yout		! solution at xout
	real(wp),dimension(neqn),intent(out)	:: ypout	! derivative of solution at xout
	
	!local variables:
	integer :: i, iq, iw, j, jq, kp1, kp2, l, m
	real(wp) :: alp, c(13), g(13), gdi, gdif, gamma, h, hi,&
						 hmu, rmu, sigma, temp1, temp2, temp3,&
						 w(13), xi, xim1, xiq

      kp1 = kold + 1
      kp2 = kold + 2

      hi = xout - ox
      h = x - ox
      xi = hi/h
      xim1 = xi - 1.d0

!   initialize w(*) for computing g(*)

      xiq = xi
      do iq = 1,kp1
        xiq = xi*xiq
        temp1 = iq*(iq+1)
        w(iq) = xiq/temp1
      end do
!
!   compute the double integral term gdi
!
      if (kold <= kgi) go to 50
      if (ivc > 0) go to 20
      gdi = 1.0d0/temp1
      m = 2
      go to 30
 20   iw = iv(ivc)
      gdi = ow(iw)
      m = kold - iw + 3
 30   if (m > kold) go to 60
      do 40 i = m,kold
 40     gdi = ow(kp2-i) - alpha(i)*gdi
      go to 60
 50   gdi = gi(kold)
!
!   compute g(*) and c(*)
!
 60   g(1) = xi
      g(2) = 0.5d0*xi*xi
      c(1) = 1.0d0
      c(2) = xi
      if (kold < 2) go to 90
      do 80 i = 2,kold
        alp = alpha(i)
        gamma = 1.0d0 + xim1*alp
        l = kp2 - i
        do 70 jq = 1,l
 70       w(jq) = gamma*w(jq) - alp*w(jq+1)
        g(i+1) = w(1)
 80     c(i+1) = gamma*c(i)
!
!   define interpolation parameters
!
 90   sigma = (w(2) - xim1*w(1))/gdi
      rmu = xim1*c(kp1)/gdi
      hmu = rmu/h
!
!   interpolate for the solution -- yout
!   and for the derivative of the solution -- ypout
!
      do 100 l = 1,neqn
        yout(l) = 0.0d0
 100    ypout(l) = 0.0d0
      do 120 j = 1,kold
        i = kp2 - j
        gdif = og(i) - og(i-1)
        temp2 = (g(i) - g(i-1)) - sigma*gdif
        temp3 = (c(i) - c(i-1)) + rmu*gdif
        do 110 l = 1,neqn
          yout(l) = yout(l) + temp2*phi(l,i)
 110      ypout(l) = ypout(l) + temp3*phi(l,i)
 120    continue
      do 130 l = 1,neqn
        yout(l) = ((1.0d0 - sigma)*oy(l) + sigma*y(l)) +&
                   h*(yout(l) + (g(1) - sigma*og(1))*phi(l,1))
 130    ypout(l) = hmu*(oy(l) - y(l)) +&
                      (ypout(l) + (c(1) + rmu*og(1))*phi(l,1))

!*****************************************************************************************    
	end subroutine dintp
!*****************************************************************************************    
      
!*****************************************************************************************    
	subroutine dsteps (me, neqn, y, x, h, eps, wt, start, hold, k,&
         kold, crash, phi, p, yp, psi, alpha, beta, sig, v, w, g,&
         phase1, ns, nornd, ksteps, twou, fouru, xold, kprev, ivc, iv,&
         kgi, gi)
!*****************************************************************************************    
!***purpose  integrate a system of first order ordinary differential
!            equations one step.
!***library   slatec (depac)
!***category  i1a1b
!***type      double precision (steps-s, dsteps-d)
!***keywords  adams method, depac, initial value problems, ode,
!             ordinary differential equations, predictor-corrector
!***author  shampine, l. f., (snla)
!           gordon, m. k., (snla)
!             modified by h.a. watts
!***description
!
!   written by l. f. shampine and m. k. gordon
!
!   abstract
!
!   subroutine  dsteps  is normally used indirectly through subroutine
!   ddeabm .  because ddeabm suffices for most problems and is much
!   easier to use, using it should be considered before using  dsteps
!   alone.
!
!   subroutine dsteps integrates a system of  neqn  first order ordinary
!   differential equations one step, normally from x to x+h, using a
!   modified divided difference form of the adams pece formulas.  local
!   extrapolation is used to improve absolute stability and accuracy.
!   the code adjusts its order and step size to control the local error
!   per unit step in a generalized sense.  special devices are included
!   to control roundoff error and to detect when the user is requesting
!   too much accuracy.
!
!   this code is completely explained and documented in the text,
!   computer solution of ordinary differential equations, the initial
!   value problem  by l. f. shampine and m. k. gordon.
!   further details on use of this code are available in "solving
!   ordinary differential equations with ode, step, and intrp",
!   by l. f. shampine and m. k. gordon, sla-73-1060.
!
!
!   the parameters represent --
!      df -- subroutine to evaluate derivatives
!      neqn -- number of equations to be integrated
!      y(*) -- solution vector at x
!      x -- independent variable
!      h -- appropriate step size for next step.  normally determined by
!           code
!      eps -- local error tolerance
!      wt(*) -- vector of weights for error criterion
!      start -- logical variable set .true. for first step,  .false.
!           otherwise
!      hold -- step size used for last successful step
!      k -- appropriate order for next step (determined by code)
!      kold -- order used for last successful step
!      crash -- logical variable set .true. when no step can be taken,
!           .false. otherwise.
!      yp(*) -- derivative of solution vector at  x  after successful
!           step
!      ksteps -- counter on attempted steps
!      twou -- 2.*u where u is machine unit roundoff quantity
!      fouru -- 4.*u where u is machine unit roundoff quantity
!   the variables x,xold,kold,kgi and ivc and the arrays y,phi,alpha,g,
!   w,p,iv and gi are required for the interpolation subroutine sintrp.
!   the remaining variables and arrays are included in the call list
!   only to eliminate local retention of variables between calls.
!
!   input to dsteps
!
!      first call --
!
!   the user must provide storage in his calling program for all arrays
!   in the call list, namely
!
!     dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12),
!    1  alpha(12),beta(12),sig(13),v(12),w(12),g(13),gi(11),iv(10)
!
!    **note**
!
!   the user must also declare  start ,  crash ,  phase1  and  nornd
!   logical variables and  df  an external subroutine, supply the
!   subroutine  df(x,y,yp)  to evaluate
!      dy(i)/dx = yp(i) = df(x,y(1),y(2),...,y(neqn))
!   and initialize only the following parameters.
!      neqn -- number of equations to be integrated
!      y(*) -- vector of initial values of dependent variables
!      x -- initial value of the independent variable
!      h -- nominal step size indicating direction of integration
!           and maximum size of step.  must be variable
!      eps -- local error tolerance per step.  must be variable
!      wt(*) -- vector of non-zero weights for error criterion
!      start -- .true.
!      yp(*) -- vector of initial derivative values
!      ksteps -- set ksteps to zero
!      twou -- 2.*u where u is machine unit roundoff quantity
!      fouru -- 4.*u where u is machine unit roundoff quantity
!   define u to be the machine unit roundoff quantity by calling
!   the function routine  d1mach,  u = d1mach(4), or by
!   computing u so that u is the smallest positive number such
!   that 1.0+u > 1.0.
!
!   dsteps  requires that the l2 norm of the vector with components
!   local error(l)/wt(l)  be less than  eps  for a successful step.  the
!   array  wt  allows the user to specify an error test appropriate
!   for his problem.  for example,
!      wt(l) = 1.0  specifies absolute error,
!            = abs(y(l))  error relative to the most recent value of the
!                 l-th component of the solution,
!            = abs(yp(l))  error relative to the most recent value of
!                 the l-th component of the derivative,
!            = max(wt(l),abs(y(l)))  error relative to the largest
!                 magnitude of l-th component obtained so far,
!            = abs(y(l))*relerr/eps + abserr/eps  specifies a mixed
!                 relative-absolute test where  relerr  is relative
!                 error,  abserr  is absolute error and  eps =
!                 max(relerr,abserr) .
!
!      subsequent calls --
!
!   subroutine  dsteps  is designed so that all information needed to
!   continue the integration, including the step size  h  and the order
!   k , is returned with each step.  with the exception of the step
!   size, the error tolerance, and the weights, none of the parameters
!   should be altered.  the array  wt  must be updated after each step
!   to maintain relative error tests like those above.  normally the
!   integration is continued just beyond the desired endpoint and the
!   solution interpolated there with subroutine  sintrp .  if it is
!   impossible to integrate beyond the endpoint, the step size may be
!   reduced to hit the endpoint since the code will not take a step
!   larger than the  h  input.  changing the direction of integration,
!   i.e., the sign of  h , requires the user set  start = .true. before
!   calling  dsteps  again.  this is the only situation in which  start
!   should be altered.
!
!   output from dsteps
!
!      successful step --
!
!   the subroutine returns after each successful step with  start  and
!   crash  set .false. .  x  represents the independent variable
!   advanced one step of length  hold  from its value on input and  y
!   the solution vector at the new value of  x .  all other parameters
!   represent information corresponding to the new  x  needed to
!   continue the integration.
!
!      unsuccessful step --
!
!   when the error tolerance is too small for the machine precision,
!   the subroutine returns without taking a step and  crash = .true. .
!   an appropriate step size and error tolerance for continuing are
!   estimated and all other information is restored as upon input
!   before returning.  to continue with the larger tolerance, the user
!   just calls the code again.  a restart is neither required nor
!   desirable.
!
!***references  l. f. shampine and m. k. gordon, solving ordinary
!                 differential equations with ode, step, and intrp,
!                 report sla-73-1060, sandia laboratories, 1973.
!***routines called  d1mach, dhstrt
!***revision history  (yymmdd)
!   740101  date written
!   890531  changed all specific intrinsics to generic.  (wrb)
!   890831  modified array declarations.  (wrb)
!   890831  revision date from version 3.2
!   891214  prologue converted to version 4.0 format.  (bab)
!   920501  reformatted the references section.  (wrb)
!
!*****************************************************************************************    

	implicit none
	
	class(ddeabm_class),intent(inout) :: me	
    real(wp),intent(inout) :: x
	real(wp),intent(inout) :: h
	real(wp),intent(inout) :: eps
	real(wp),intent(inout) :: hold
	integer,intent(inout)  :: k
    integer,intent(inout)  :: kold 
    integer,intent(inout)  :: ns
    integer,intent(inout)  :: ksteps
    real(wp),intent(inout) :: twou
    real(wp),intent(inout) :: fouru
    real(wp),intent(inout) :: xold
    integer,intent(inout)  :: kprev
    integer,intent(inout)  :: ivc
    integer,intent(inout)  :: kgi
    integer,intent(in)	   :: neqn
    logical,intent(inout)  :: start,crash,phase1,nornd
    real(wp),dimension(neqn),intent(inout) :: y
    real(wp),dimension(neqn),intent(inout) :: wt
    real(wp),dimension(neqn,16),intent(inout) :: phi
    real(wp),dimension(neqn),intent(inout) :: p
    real(wp),dimension(neqn),intent(inout) :: yp
    real(wp),dimension(12),intent(inout) :: psi
	real(wp),dimension(12),intent(inout) :: alpha
	real(wp),dimension(12),intent(inout) :: beta
	real(wp),dimension(13),intent(inout) :: sig
	real(wp),dimension(12),intent(inout) :: v
	real(wp),dimension(12),intent(inout) :: w
	real(wp),dimension(13),intent(inout) :: g
	real(wp),dimension(11),intent(inout) :: gi
	integer, dimension(10),intent(inout)  :: iv

	integer i, ifail, im1, ip1, iq, j, km1, km2, knew,&
		kp1, kp2, l, limit1, limit2, nsm2,&
		nsp1, nsp2, jv
	real(wp) absh, big,&
		erk, erkm1, erkm2, erkp1, err,&
		hnew, p5eps, r,&
		reali, realns, rho, round, tau, temp1,&
		temp2, temp3, temp4, temp5, temp6, u
	  
	!parameters:            
	real(wp),dimension(13),parameter :: two = [	2.0_wp, &
												4.0_wp, &
												8.0_wp, &
												16.0_wp, &
												32.0_wp, &
												64.0_wp, &
												128.0_wp, &
												256.0_wp, &
												512.0_wp, &
												1024.0_wp, &
												2048.0_wp, &
												4096.0_wp, &
												8192.0_wp]
 
 	!note: this is a modification of the original code.
 	! The full-precision coefficients are used here, instead 
 	!	of the less precise ones in the original.
 	! These were computed from the equation on p. 159 of Shampine/Gordon, 
 	!	"Computer Solution of Ordinary Differential Equations", 1975.
  	real(wp),dimension(13),parameter :: gstr = [0.5000000000000000E+00_wp, &
  												0.8333333333333331E-01_wp, &
  												0.4166666666666669E-01_wp, &
  												0.2638888888888891E-01_wp, &
  												0.1874999999999996E-01_wp, &
  												0.1426917989417992E-01_wp, &
  												0.1136739417989419E-01_wp, &
  												0.9356536596119916E-02_wp, &
  												0.7892554012345690E-02_wp, &
  												0.6785849984634704E-02_wp, &
  												0.5924056412337661E-02_wp, &
  												0.5236693257950287E-02_wp, &
  												0.4677498407042263E-02_wp]
 												         
!
!       ***     begin block 0     ***
!   check if step size or error tolerance is too small for machine
!   precision.  if first step, initialize phi array and estimate a
!   starting step size.
!                   ***

!   if step size is too small, determine an acceptable one
      crash = .true.
      if (abs(h) < fouru*abs(x)) then
      	h = sign(fouru*abs(x),h)
      	return
      end if
      
      p5eps = 0.5d0*eps

!   if error tolerance is too small, increase it to an acceptable value
      round = 0.0d0
      do l = 1,neqn
        round = round + (y(l)/wt(l))**2
 	  end do
      round = twou*sqrt(round)
      if (p5eps < round) then
      	eps = 2.0d0*round*(1.0d0 + fouru)
      	return
      end if
      
      crash = .false.
      g(1) = 1.0d0
      g(2) = 0.5d0
      sig(1) = 1.0d0
      
      if (start) then

	!   initialize.  compute appropriate step size for first step
	!     call me%df(x,y,yp)
	!     sum = 0.0
		  do l = 1,neqn
			phi(l,1) = yp(l)
			phi(l,2) = 0.0d0
		  end do
	!       sum = sum + (yp(l)/wt(l))**2
	!     end do
	!     sum = sqrt(sum)
	!     absh = abs(h)
	!     if (eps < 16.0*sum*h*h) absh = 0.25*sqrt(eps/sum)
	!     h = sign(max(absh,fouru*abs(x)),h)
	!
		  u = d1mach4
		  big = sqrt(d1mach2)
		  call me%dhstrt(neqn,x,x+h,y,yp,wt,1,u,big,&
					   phi(1,3),phi(1,4),phi(1,5),phi(1,6),h)

		  hold = 0.0d0
		  k = 1
		  kold = 0
		  kprev = 0
		  start = .false.
		  phase1 = .true.
		  nornd = .true.
		  if (p5eps <= 100.0d0*round) then
			  nornd = .false.
			  do l = 1,neqn
				phi(l,15) = 0.0d0
			  end do
		  end if
	  
	  end if
	  
      ifail = 0
 
!       ***     end block 0     ***

!       ***     begin block 1     ***
!   compute coefficients of formulas for this step.  avoid computing
!   those quantities not changed when step size is not changed.
!                   ***
!
 100  kp1 = k+1
      kp2 = k+2
      km1 = k-1
      km2 = k-2
!
!   ns is the number of dsteps taken with size h, including the current
!   one.  when k<ns, no coefficients change
!
      if (h /= hold) ns = 0
      if (ns<=kold) ns = ns+1
      nsp1 = ns+1
      if (k < ns) go to 199
!
!   compute those components of alpha(*),beta(*),psi(*),sig(*) which
!   are changed
!
      beta(ns) = 1.0d0
      realns = ns
      alpha(ns) = 1.0d0/realns
      temp1 = h*realns
      sig(nsp1) = 1.0d0
      if (k < nsp1) go to 110
      do 105 i = nsp1,k
        im1 = i-1
        temp2 = psi(im1)
        psi(im1) = temp1
        beta(i) = beta(im1)*psi(im1)/temp2
        temp1 = temp2 + h
        alpha(i) = h/temp1
        reali = i
 105    sig(i+1) = reali*alpha(i)*sig(i)
 110  psi(k) = temp1
!
!   compute coefficients g(*)
!
!   initialize v(*) and set w(*).
!
      if (ns > 1) go to 120
      do 115 iq = 1,k
        temp3 = iq*(iq+1)
        v(iq) = 1.0d0/temp3
 115    w(iq) = v(iq)
      ivc = 0
      kgi = 0
      if (k == 1) go to 140
      kgi = 1
      gi(1) = w(2)
      go to 140
!
!   if order was raised, update diagonal part of v(*)
!
 120  if (k <= kprev) go to 130
      if (ivc == 0) go to 122
      jv = kp1 - iv(ivc)
      ivc = ivc - 1
      go to 123
 122  jv = 1
      temp4 = k*kp1
      v(k) = 1.0d0/temp4
      w(k) = v(k)
      if (k /= 2) go to 123
      kgi = 1
      gi(1) = w(2)
 123  nsm2 = ns-2
      if (nsm2 < jv) go to 130
      do 125 j = jv,nsm2
        i = k-j
        v(i) = v(i) - alpha(j+1)*v(i+1)
 125    w(i) = v(i)
      if (i /= 2) go to 130
      kgi = ns - 1
      gi(kgi) = w(2)
!
!   update v(*) and set w(*)
!
 130  limit1 = kp1 - ns
      temp5 = alpha(ns)
      do 135 iq = 1,limit1
        v(iq) = v(iq) - temp5*v(iq+1)
 135    w(iq) = v(iq)
      g(nsp1) = w(1)
      if (limit1 == 1) go to 137
      kgi = ns
      gi(kgi) = w(2)
 137  w(limit1+1) = v(limit1+1)
      if (k >= kold) go to 140
      ivc = ivc + 1
      iv(ivc) = limit1 + 2
!
!   compute the g(*) in the work vector w(*)
!
 140  nsp2 = ns + 2
      kprev = k
      if (kp1 < nsp2) go to 199
      do 150 i = nsp2,kp1
        limit2 = kp2 - i
        temp6 = alpha(i-1)
        do 145 iq = 1,limit2
 145      w(iq) = w(iq) - temp6*w(iq+1)
 150    g(i) = w(1)
 199    continue
!       ***     end block 1     ***
!
!       ***     begin block 2     ***
!   predict a solution p(*), evaluate derivatives using predicted
!   solution, estimate local error at order k and errors at orders k,
!   k-1, k-2 as if constant step size were used.
!                   ***
!
!   increment counter on attempted dsteps
!
      ksteps = ksteps + 1
!
!   change phi to phi star
!
      if (k < nsp1) go to 215
      do 210 i = nsp1,k
        temp1 = beta(i)
        do 205 l = 1,neqn
 205      phi(l,i) = temp1*phi(l,i)
 210    continue
!
!   predict solution and differences
!
 215  do 220 l = 1,neqn
        phi(l,kp2) = phi(l,kp1)
        phi(l,kp1) = 0.0d0
 220    p(l) = 0.0d0
      do 230 j = 1,k
        i = kp1 - j
        ip1 = i+1
        temp2 = g(i)
        do 225 l = 1,neqn
          p(l) = p(l) + temp2*phi(l,i)
 225      phi(l,i) = phi(l,i) + phi(l,ip1)
 230    continue
      if (nornd) go to 240
      do 235 l = 1,neqn
        tau = h*p(l) - phi(l,15)
        p(l) = y(l) + tau
 235    phi(l,16) = (p(l) - y(l)) - tau
      go to 250
 240  do 245 l = 1,neqn
 245    p(l) = y(l) + h*p(l)
 250  xold = x
      x = x + h
      absh = abs(h)
      call me%df(x,p,yp)
!
!   estimate errors at orders k,k-1,k-2
!
      erkm2 = 0.0d0
      erkm1 = 0.0d0
      erk = 0.0d0
      do 265 l = 1,neqn
        temp3 = 1.0d0/wt(l)
        temp4 = yp(l) - phi(l,1)
        if (km2)265,260,255
 255    erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
 260    erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
 265    erk = erk + (temp4*temp3)**2
      if (km2)280,275,270
 270  erkm2 = absh*sig(km1)*gstr(km2)*sqrt(erkm2)
 275  erkm1 = absh*sig(k)*gstr(km1)*sqrt(erkm1)
 280  temp5 = absh*sqrt(erk)
      err = temp5*(g(k)-g(kp1))
      erk = temp5*sig(kp1)*gstr(k)
      knew = k
!
!   test if order should be lowered
!
      if (km2)299,290,285
 285  if (max(erkm1,erkm2) <= erk) knew = km1
      go to 299
 290  if (erkm1 <= 0.5d0*erk) knew = km1
!
!   test if step successful
!
 299  if (err <= eps) go to 400
!       ***     end block 2     ***
!
!       ***     begin block 3     ***
!   the step is unsuccessful.  restore  x, phi(*,*), psi(*) .
!   if third consecutive failure, set order to one.  if step fails more
!   than three times, consider an optimal step size.  double error
!   tolerance and return if estimated step size is too small for machine
!   precision.
!                   ***
!
!   restore x, phi(*,*) and psi(*)
!
      phase1 = .false.
      x = xold
      do 310 i = 1,k
        temp1 = 1.0d0/beta(i)
        ip1 = i+1
        do 305 l = 1,neqn
 305      phi(l,i) = temp1*(phi(l,i) - phi(l,ip1))
 310    continue
      if (k < 2) go to 320
      do 315 i = 2,k
 315    psi(i-1) = psi(i) - h
!
!   on third failure, set order to one.  thereafter, use optimal step
!   size
!
 320  ifail = ifail + 1
      temp2 = 0.5d0
      if (ifail - 3) 335,330,325
 325  if (p5eps < 0.25d0*erk) temp2 = sqrt(p5eps/erk)
 330  knew = 1
 335  h = temp2*h
      k = knew
      ns = 0
      if (abs(h) >= fouru*abs(x)) go to 340
      crash = .true.
      h = sign(fouru*abs(x),h)
      eps = eps + eps
      return
 340  go to 100
!       ***     end block 3     ***
!
!       ***     begin block 4     ***
!   the step is successful.  correct the predicted solution, evaluate
!   the derivatives using the corrected solution and update the
!   differences.  determine best order and step size for next step.
!                   ***
 400  kold = k
      hold = h
!
!   correct and evaluate
!
      temp1 = h*g(kp1)
      if (nornd) go to 410
      do 405 l = 1,neqn
        temp3 = y(l)
        rho = temp1*(yp(l) - phi(l,1)) - phi(l,16)
        y(l) = p(l) + rho
        phi(l,15) = (y(l) - p(l)) - rho
 405    p(l) = temp3
      go to 420
 410  do 415 l = 1,neqn
        temp3 = y(l)
        y(l) = p(l) + temp1*(yp(l) - phi(l,1))
 415    p(l) = temp3
 420  call me%df(x,y,yp)
!
!   update differences for next step
!
      do 425 l = 1,neqn
        phi(l,kp1) = yp(l) - phi(l,1)
 425    phi(l,kp2) = phi(l,kp1) - phi(l,kp2)
      do 435 i = 1,k
        do 430 l = 1,neqn
 430      phi(l,i) = phi(l,i) + phi(l,kp1)
 435    continue
!
!   estimate error at order k+1 unless:
!     in first phase when always raise order,
!     already decided to lower order,
!     step size not constant so estimate unreliable
!
      erkp1 = 0.0d0
      if (knew == km1  .or.  k == 12) phase1 = .false.
      if (phase1) go to 450
      if (knew == km1) go to 455
      if (kp1 > ns) go to 460
      do 440 l = 1,neqn
 440    erkp1 = erkp1 + (phi(l,kp2)/wt(l))**2
      erkp1 = absh*gstr(kp1)*sqrt(erkp1)
!
!   using estimated error at order k+1, determine appropriate order
!   for next step
!
      if (k > 1) go to 445
      if (erkp1 >= 0.5d0*erk) go to 460
      go to 450
 445  if (erkm1 <= min(erk,erkp1)) go to 455
      if (erkp1 >= erk  .or.  k == 12) go to 460
!
!   here erkp1 < erk < max(erkm1,erkm2) else order would have
!   been lowered in block 2.  thus order is to be raised
!
!   raise order
!
 450  k = kp1
      erk = erkp1
      go to 460
!
!   lower order
!
 455  k = km1
      erk = erkm1
!
!   with new order determine appropriate step size for next step
!
 460  hnew = h + h
      if (phase1) go to 465
      if (p5eps >= erk*two(k+1)) go to 465
      hnew = h
      if (p5eps >= erk) go to 465
      temp2 = k+1
      r = (p5eps/erk)**(1.0d0/temp2)
      hnew = absh*max(0.5d0,min(0.9d0,r))
      hnew = sign(max(hnew,fouru*abs(x)),h)
 465  h = hnew
!       ***     end block 4     ***

!*****************************************************************************************    
	end subroutine dsteps
!*****************************************************************************************    

!*****************************************************************************************    
    subroutine xermsg (librar, subrou, messg, nerr, level)
!*****************************************************************************************    

    implicit none

    character(len=*),intent(in) :: librar, subrou, messg
    integer,intent(in) :: nerr, level

    write(*,'(A)') 'Error in '//trim(subrou)//': '//trim(messg)      

!*****************************************************************************************    
    end subroutine xermsg
!*****************************************************************************************    

!*****************************************************************************************    
    end module ddeabm_module
!*****************************************************************************************    
