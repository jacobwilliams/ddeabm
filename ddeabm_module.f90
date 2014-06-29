    module ddeabm_module
    
    use, intrinsic :: iso_fortran_env, wp=>real64    !double precision
    
    private
    
    !parameters:    
    real(wp),parameter :: d1mach2 = huge(1.0_wp)                        ! the largest magnitude
    real(wp),parameter :: d1mach4 = RADIX(1.0_wp)**(1-DIGITS(1.0_wp))   ! the largest relative spacing

    !public routines:
    public :: ddeabm
    
    contains

!DECK DDEABM
      SUBROUTINE DDEABM (DF, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID,&
         RWORK, LRW, IWORK, LIW, RPAR, IPAR)
!***BEGIN PROLOGUE  DDEABM
!***PURPOSE  Solve an initial value problem in ordinary differential
!            equations using an Adams-Bashforth method.
!***LIBRARY   SLATEC (DEPAC)
!***CATEGORY  I1A1B
!***TYPE      DOUBLE PRECISION (DEABM-S, DDEABM-D)
!***KEYWORDS  ADAMS-BASHFORTH METHOD, DEPAC, INITIAL VALUE PROBLEMS,
!             ODE, ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR
!***AUTHOR  Shampine, L. F., (SNLA)
!           Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   This is the Adams code in the package of differential equation
!   solvers DEPAC, consisting of the codes DDERKF, DDEABM, and DDEBDF.
!   Design of the package was by L. F. Shampine and H. A. Watts.
!   It is documented in
!        SAND79-2374 , DEPAC - Design of a User Oriented Package of ODE
!                              Solvers.
!   DDEABM is a driver for a modification of the code ODE written by
!             L. F. Shampine and M. K. Gordon
!             Sandia Laboratories
!             Albuquerque, New Mexico 87185
!
! **********************************************************************
! * ABSTRACT *
! ************
!
!   Subroutine DDEABM uses the Adams-Bashforth-Moulton
!   Predictor-Corrector formulas of orders one through twelve to
!   integrate a system of NEQ first order ordinary differential
!   equations of the form
!                         DU/DX = DF(X,U)
!   when the vector Y(*) of initial values for U(*) at X=T is given.
!   The subroutine integrates from T to TOUT. It is easy to continue the
!   integration to get results at additional TOUT.  This is the interval
!   mode of operation.  It is also easy for the routine to return with
!   the solution at each intermediate step on the way to TOUT.  This is
!   the intermediate-output mode of operation.
!
!   DDEABM uses subprograms DDES, DSTEPS, DINTP, DHSTRT, DHVNRM,
!   D1MACH, and the error handling routine XERMSG.  The only machine
!   dependent parameters to be assigned appear in D1MACH.
!
! **********************************************************************
! * Description of The Arguments To DDEABM (An Overview) *
! **********************************************************************
!
!   The Parameters are
!
!      DF -- This is the name of a subroutine which you provide to
!             define the differential equations.
!
!      NEQ -- This is the number of (first order) differential
!             equations to be integrated.
!
!      T -- This is a DOUBLE PRECISION value of the independent
!           variable.
!
!      Y(*) -- This DOUBLE PRECISION array contains the solution
!              components at T.
!
!      TOUT -- This is a DOUBLE PRECISION point at which a solution is
!              desired.
!
!      INFO(*) -- The basic task of the code is to integrate the
!             differential equations from T to TOUT and return an
!             answer at TOUT.  INFO(*) is an INTEGER array which is used
!             to communicate exactly how you want this task to be
!             carried out.
!
!      RTOL, ATOL -- These DOUBLE PRECISION quantities represent
!                    relative and absolute error tolerances which you
!                    provide to indicate how accurately you wish the
!                    solution to be computed.  You may choose them to be
!                    both scalars or else both vectors.
!
!      IDID -- This scalar quantity is an indicator reporting what
!             the code did.  You must monitor this INTEGER variable to
!             decide what action to take next.
!
!      RWORK(*), LRW -- RWORK(*) is a DOUBLE PRECISION work array of
!             length LRW which provides the code with needed storage
!             space.
!
!      IWORK(*), LIW -- IWORK(*) is an INTEGER work array of length LIW
!             which provides the code with needed storage space and an
!             across call flag.
!
!      RPAR, IPAR -- These are DOUBLE PRECISION and INTEGER parameter
!             arrays which you can use for communication between your
!             calling program and the DF subroutine.
!
!  Quantities which are used as input items are
!             NEQ, T, Y(*), TOUT, INFO(*),
!             RTOL, ATOL, RWORK(1), LRW and LIW.
!
!  Quantities which may be altered by the code are
!             T, Y(*), INFO(1), RTOL, ATOL,
!             IDID, RWORK(*) and IWORK(*).
!
! **********************************************************************
! * INPUT -- What To Do On The First Call To DDEABM *
! **********************************************************************
!
!   The first call of the code is defined to be the start of each new
!   problem.  Read through the descriptions of all the following items,
!   provide sufficient storage space for designated arrays, set
!   appropriate variables for the initialization of the problem, and
!   give information about how you want the problem to be solved.
!
!
!      DF -- Provide a subroutine of the form
!                               DF(X,U,UPRIME,RPAR,IPAR)
!             to define the system of first order differential equations
!             which is to be solved.  For the given values of X and the
!             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
!             evaluate the NEQ components of the system of differential
!             equations  DU/DX=DF(X,U)  and store the derivatives in the
!             array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for
!             equations I=1,...,NEQ.
!
!             Subroutine DF must NOT alter X or U(*).  You must declare
!             the name df in an external statement in your program that
!             calls DDEABM.  You must dimension U and UPRIME in DF.
!
!             RPAR and IPAR are DOUBLE PRECISION and INTEGER parameter
!             arrays which you can use for communication between your
!             calling program and subroutine DF. They are not used or
!             altered by DDEABM.  If you do not need RPAR or IPAR,
!             ignore these parameters by treating them as dummy
!             arguments. If you do choose to use them, dimension them in
!             your calling program and in DF as arrays of appropriate
!             length.
!
!      NEQ -- Set it to the number of differential equations.
!             (NEQ >= 1)
!
!      T -- Set it to the initial point of the integration.
!             You must use a program variable for T because the code
!             changes its value.
!
!      Y(*) -- Set this vector to the initial values of the NEQ solution
!             components at the initial point.  You must dimension Y at
!             least NEQ in your calling program.
!
!      TOUT -- Set it to the first point at which a solution
!             is desired.  You can take TOUT = T, in which case the code
!             will evaluate the derivative of the solution at T and
!             return. Integration either forward in T  (TOUT > T)  or
!             backward in T  (TOUT < T)  is permitted.
!
!             The code advances the solution from T to TOUT using
!             step sizes which are automatically selected so as to
!             achieve the desired accuracy.  If you wish, the code will
!             return with the solution and its derivative following
!             each intermediate step (intermediate-output mode) so that
!             you can monitor them, but you still must provide TOUT in
!             accord with the basic aim of the code.
!
!             The first step taken by the code is a critical one
!             because it must reflect how fast the solution changes near
!             the initial point.  The code automatically selects an
!             initial step size which is practically always suitable for
!             the problem. By using the fact that the code will not step
!             past TOUT in the first step, you could, if necessary,
!             restrict the length of the initial step size.
!
!             For some problems it may not be permissible to integrate
!             past a point TSTOP because a discontinuity occurs there
!             or the solution or its derivative is not defined beyond
!             TSTOP.  When you have declared a TSTOP point (see INFO(4)
!             and RWORK(1)), you have told the code not to integrate
!             past TSTOP.  In this case any TOUT beyond TSTOP is invalid
!             input.
!
!      INFO(*) -- Use the INFO array to give the code more details about
!             how you want your problem solved.  This array should be
!             dimensioned of length 15 to accommodate other members of
!             DEPAC or possible future extensions, though DDEABM uses
!             only the first four entries.  You must respond to all of
!             the following items which are arranged as questions.  The
!             simplest use of the code corresponds to answering all
!             questions as YES ,i.e. setting ALL entries of INFO to 0.
!
!        INFO(1) -- This parameter enables the code to initialize
!               itself.  You must set it to indicate the start of every
!               new problem.
!
!            **** Is this the first call for this problem ...
!                  YES -- set INFO(1) = 0
!                   NO -- not applicable here.
!                         See below for continuation calls.  ****
!
!        INFO(2) -- How much accuracy you want of your solution
!               is specified by the error tolerances RTOL and ATOL.
!               The simplest use is to take them both to be scalars.
!               To obtain more flexibility, they can both be vectors.
!               The code must be told your choice.
!
!            **** Are both error tolerances RTOL, ATOL scalars ...
!                  YES -- set INFO(2) = 0
!                         and input scalars for both RTOL and ATOL
!                   NO -- set INFO(2) = 1
!                         and input arrays for both RTOL and ATOL ****
!
!        INFO(3) -- The code integrates from T in the direction
!               of TOUT by steps.  If you wish, it will return the
!               computed solution and derivative at the next
!               intermediate step (the intermediate-output mode) or
!               TOUT, whichever comes first.  This is a good way to
!               proceed if you want to see the behavior of the solution.
!               If you must have solutions at a great many specific
!               TOUT points, this code will compute them efficiently.
!
!            **** Do you want the solution only at
!                 TOUT (and not at the next intermediate step) ...
!                  YES -- set INFO(3) = 0
!                   NO -- set INFO(3) = 1 ****
!
!        INFO(4) -- To handle solutions at a great many specific
!               values TOUT efficiently, this code may integrate past
!               TOUT and interpolate to obtain the result at TOUT.
!               Sometimes it is not possible to integrate beyond some
!               point TSTOP because the equation changes there or it is
!               not defined past TSTOP.  Then you must tell the code
!               not to go past.
!
!            **** Can the integration be carried out without any
!                 Restrictions on the independent variable T ...
!                  YES -- set INFO(4)=0
!                   NO -- set INFO(4)=1
!                         and define the stopping point TSTOP by
!                         setting RWORK(1)=TSTOP ****
!
!      RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL)
!             error tolerances to tell the code how accurately you want
!             the solution to be computed.  They must be defined as
!             program variables because the code may change them.  You
!             have two choices --
!                  Both RTOL and ATOL are scalars. (INFO(2)=0)
!                  Both RTOL and ATOL are vectors. (INFO(2)=1)
!             In either case all components must be non-negative.
!
!             The tolerances are used by the code in a local error test
!             at each step which requires roughly that
!                     ABS(LOCAL ERROR) <= RTOL*ABS(Y)+ATOL
!             for each vector component.
!             (More specifically, a Euclidean norm is used to measure
!             the size of vectors, and the error test uses the magnitude
!             of the solution at the beginning of the step.)
!
!             The true (global) error is the difference between the true
!             solution of the initial value problem and the computed
!             approximation.  Practically all present day codes,
!             including this one, control the local error at each step
!             and do not even attempt to control the global error
!             directly.  Roughly speaking, they produce a solution Y(T)
!             which satisfies the differential equations with a
!             residual R(T),    DY(T)/DT = DF(T,Y(T)) + R(T)   ,
!             and, almost always, R(T) is bounded by the error
!             tolerances.  Usually, but not always, the true accuracy of
!             the computed Y is comparable to the error tolerances. This
!             code will usually, but not always, deliver a more accurate
!             solution if you reduce the tolerances and integrate again.
!             By comparing two such solutions you can get a fairly
!             reliable idea of the true error in the solution at the
!             bigger tolerances.
!
!             Setting ATOL=0.D0 results in a pure relative error test on
!             that component. Setting RTOL=0. results in a pure absolute
!             error test on that component.  A mixed test with non-zero
!             RTOL and ATOL corresponds roughly to a relative error
!             test when the solution component is much bigger than ATOL
!             and to an absolute error test when the solution component
!             is smaller than the threshold ATOL.
!
!             Proper selection of the absolute error control parameters
!             ATOL  requires you to have some idea of the scale of the
!             solution components.  To acquire this information may mean
!             that you will have to solve the problem more than once. In
!             the absence of scale information, you should ask for some
!             relative accuracy in all the components (by setting  RTOL
!             values non-zero) and perhaps impose extremely small
!             absolute error tolerances to protect against the danger of
!             a solution component becoming zero.
!
!             The code will not attempt to compute a solution at an
!             accuracy unreasonable for the machine being used.  It will
!             advise you if you ask for too much accuracy and inform
!             you as to the maximum accuracy it believes possible.
!
!      RWORK(*) -- Dimension this DOUBLE PRECISION work array of length
!             LRW in your calling program.
!
!      RWORK(1) -- If you have set INFO(4)=0, you can ignore this
!             optional input parameter.  Otherwise you must define a
!             stopping point TSTOP by setting   RWORK(1) = TSTOP.
!             (for some problems it may not be permissible to integrate
!             past a point TSTOP because a discontinuity occurs there
!             or the solution or its derivative is not defined beyond
!             TSTOP.)
!
!      LRW -- Set it to the declared length of the RWORK array.
!             You must have  LRW >= 130+21*NEQ
!
!      IWORK(*) -- Dimension this INTEGER work array of length LIW in
!             your calling program.
!
!      LIW -- Set it to the declared length of the IWORK array.
!             You must have  LIW >= 51
!
!      RPAR, IPAR -- These are parameter arrays, of DOUBLE PRECISION and
!             INTEGER type, respectively.  You can use them for
!             communication between your program that calls DDEABM and
!             the  DF subroutine.  They are not used or altered by
!             DDEABM.  If you do not need RPAR or IPAR, ignore these
!             parameters by treating them as dummy arguments.  If you do
!             choose to use them, dimension them in your calling program
!             and in DF as arrays of appropriate length.
!
! **********************************************************************
! * OUTPUT -- After Any Return From DDEABM *
! **********************************************************************
!
!   The principal aim of the code is to return a computed solution at
!   TOUT, although it is also possible to obtain intermediate results
!   along the way.  To find out whether the code achieved its goal
!   or if the integration process was interrupted before the task was
!   completed, you must check the IDID parameter.
!
!
!      T -- The solution was successfully advanced to the
!             output value of T.
!
!      Y(*) -- Contains the computed solution approximation at T.
!             You may also be interested in the approximate derivative
!             of the solution at T.  It is contained in
!             RWORK(21),...,RWORK(20+NEQ).
!
!      IDID -- Reports what the code did
!
!                         *** Task Completed ***
!                   Reported by positive values of IDID
!
!             IDID = 1 -- A step was successfully taken in the
!                       intermediate-output mode.  The code has not
!                       yet reached TOUT.
!
!             IDID = 2 -- The integration to TOUT was successfully
!                       completed (T=TOUT) by stepping exactly to TOUT.
!
!             IDID = 3 -- The integration to TOUT was successfully
!                       completed (T=TOUT) by stepping past TOUT.
!                       Y(*) is obtained by interpolation.
!
!                         *** Task Interrupted ***
!                   Reported by negative values of IDID
!
!             IDID = -1 -- A large amount of work has been expended.
!                       (500 steps attempted)
!
!             IDID = -2 -- The error tolerances are too stringent.
!
!             IDID = -3 -- The local error test cannot be satisfied
!                       because you specified a zero component in ATOL
!                       and the corresponding computed solution
!                       component is zero.  Thus, a pure relative error
!                       test is impossible for this component.
!
!             IDID = -4 -- The problem appears to be stiff.
!
!             IDID = -5,-6,-7,..,-32  -- Not applicable for this code
!                       but used by other members of DEPAC or possible
!                       future extensions.
!
!                         *** Task Terminated ***
!                   Reported by the value of IDID=-33
!
!             IDID = -33 -- The code has encountered trouble from which
!                       it cannot recover.  A message is printed
!                       explaining the trouble and control is returned
!                       to the calling program. For example, this occurs
!                       when invalid input is detected.
!
!      RTOL, ATOL -- These quantities remain unchanged except when
!             IDID = -2.  In this case, the error tolerances have been
!             increased by the code to values which are estimated to be
!             appropriate for continuing the integration.  However, the
!             reported solution at T was obtained using the input values
!             of RTOL and ATOL.
!
!      RWORK, IWORK -- Contain information which is usually of no
!             interest to the user but necessary for subsequent calls.
!             However, you may find use for
!
!             RWORK(11)--which contains the step size H to be
!                        attempted on the next step.
!
!             RWORK(12)--if the tolerances have been increased by the
!                        code (IDID = -2) , they were multiplied by the
!                        value in RWORK(12).
!
!             RWORK(13)--Which contains the current value of the
!                        independent variable, i.e. the farthest point
!                        integration has reached. This will be different
!                        from T only when interpolation has been
!                        performed (IDID=3).
!
!             RWORK(20+I)--Which contains the approximate derivative
!                        of the solution component Y(I).  In DDEABM, it
!                        is obtained by calling subroutine DF to
!                        evaluate the differential equation using T and
!                        Y(*) when IDID=1 or 2, and by interpolation
!                        when IDID=3.
!
! **********************************************************************
! * INPUT -- What To Do To Continue The Integration *
! *             (calls after the first)             *
! **********************************************************************
!
!        This code is organized so that subsequent calls to continue the
!        integration involve little (if any) additional effort on your
!        part. You must monitor the IDID parameter in order to determine
!        what to do next.
!
!        Recalling that the principal task of the code is to integrate
!        from T to TOUT (the interval mode), usually all you will need
!        to do is specify a new TOUT upon reaching the current TOUT.
!
!        Do not alter any quantity not specifically permitted below,
!        in particular do not alter NEQ, T, Y(*), RWORK(*), IWORK(*) or
!        the differential equation in subroutine DF. Any such alteration
!        constitutes a new problem and must be treated as such, i.e.
!        you must start afresh.
!
!        You cannot change from vector to scalar error control or vice
!        versa (INFO(2)) but you can change the size of the entries of
!        RTOL, ATOL.  Increasing a tolerance makes the equation easier
!        to integrate.  Decreasing a tolerance will make the equation
!        harder to integrate and should generally be avoided.
!
!        You can switch from the intermediate-output mode to the
!        interval mode (INFO(3)) or vice versa at any time.
!
!        If it has been necessary to prevent the integration from going
!        past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
!        code will not integrate to any TOUT beyond the currently
!        specified TSTOP.  Once TSTOP has been reached you must change
!        the value of TSTOP or set INFO(4)=0.  You may change INFO(4)
!        or TSTOP at any time but you must supply the value of TSTOP in
!        RWORK(1) whenever you set INFO(4)=1.
!
!        The parameter INFO(1) is used by the code to indicate the
!        beginning of a new problem and to indicate whether integration
!        is to be continued.  You must input the value  INFO(1) = 0
!        when starting a new problem.  You must input the value
!        INFO(1) = 1  if you wish to continue after an interrupted task.
!        Do not set  INFO(1) = 0  on a continuation call unless you
!        want the code to restart at the current T.
!
!                         *** Following A Completed Task ***
!         If
!             IDID = 1, call the code again to continue the integration
!                     another step in the direction of TOUT.
!
!             IDID = 2 or 3, define a new TOUT and call the code again.
!                     TOUT must be different from T. You cannot change
!                     the direction of integration without restarting.
!
!                         *** Following An Interrupted Task ***
!                     To show the code that you realize the task was
!                     interrupted and that you want to continue, you
!                     must take appropriate action and reset INFO(1) = 1
!         If
!             IDID = -1, the code has attempted 500 steps.
!                     If you want to continue, set INFO(1) = 1 and
!                     call the code again. An additional 500 steps
!                     will be allowed.
!
!             IDID = -2, the error tolerances RTOL, ATOL have been
!                     increased to values the code estimates appropriate
!                     for continuing.  You may want to change them
!                     yourself.  If you are sure you want to continue
!                     with relaxed error tolerances, set INFO(1)=1 and
!                     call the code again.
!
!             IDID = -3, a solution component is zero and you set the
!                     corresponding component of ATOL to zero.  If you
!                     are sure you want to continue, you must first
!                     alter the error criterion to use positive values
!                     for those components of ATOL corresponding to zero
!                     solution components, then set INFO(1)=1 and call
!                     the code again.
!
!             IDID = -4, the problem appears to be stiff.  It is very
!                     inefficient to solve such problems with DDEABM.
!                     The code DDEBDF in DEPAC handles this task
!                     efficiently.  If you are absolutely sure you want
!                     to continue with DDEABM, set INFO(1)=1 and call
!                     the code again.
!
!             IDID = -5,-6,-7,..,-32  --- cannot occur with this code
!                     but used by other members of DEPAC or possible
!                     future extensions.
!
!                         *** Following A Terminated Task ***
!         If
!             IDID = -33, you cannot continue the solution of this
!                     problem.  An attempt to do so will result in your
!                     run being terminated.
!
! **********************************************************************
! *Long Description:
!
! **********************************************************************
! *             DEPAC Package Overview           *
! **********************************************************************
!
! ....   You have a choice of three differential equation solvers from
! ....   DEPAC. The following brief descriptions are meant to aid you in
! ....   choosing the most appropriate code for your problem.
!
! ....   DDERKF is a fifth order Runge-Kutta code. It is the simplest of
! ....   the three choices, both algorithmically and in the use of the
! ....   code. DDERKF is primarily designed to solve non-stiff and
! ....   mildly stiff differential equations when derivative evaluations
! ....   are not expensive. It should generally not be used to get high
! ....   accuracy results nor answers at a great many specific points.
! ....   Because DDERKF has very low overhead costs, it will usually
! ....   result in the least expensive integration when solving
! ....   problems requiring a modest amount of accuracy and having
! ....   equations that are not costly to evaluate. DDERKF attempts to
! ....   discover when it is not suitable for the task posed.
!
! ....   DDEABM is a variable order (one through twelve) Adams code.
! ....   Its complexity lies somewhere between that of DDERKF and
! ....   DDEBDF.  DDEABM is primarily designed to solve non-stiff and
! ....   mildly stiff differential equations when derivative evaluations
! ....   are expensive, high accuracy results are needed or answers at
! ....   many specific points are required. DDEABM attempts to discover
! ....   when it is not suitable for the task posed.
!
! ....   DDEBDF is a variable order (one through five) backward
! ....   differentiation formula code. it is the most complicated of
! ....   the three choices. DDEBDF is primarily designed to solve stiff
! ....   differential equations at crude to moderate tolerances.
! ....   If the problem is very stiff at all, DDERKF and DDEABM will be
! ....   quite inefficient compared to DDEBDF. However, DDEBDF will be
! ....   inefficient compared to DDERKF and DDEABM on non-stiff problems
! ....   because it uses much more storage, has a much larger overhead,
! ....   and the low order formulas will not give high accuracies
! ....   efficiently.
!
! ....   The concept of stiffness cannot be described in a few words.
! ....   If you do not know the problem to be stiff, try either DDERKF
! ....   or DDEABM. Both of these codes will inform you of stiffness
! ....   when the cost of solving such problems becomes important.
!
! *********************************************************************
!
!***REFERENCES  L. F. Shampine and H. A. Watts, DEPAC - design of a user
!                 oriented package of ODE solvers, Report SAND79-2374,
!                 Sandia Laboratories, 1979.
!***ROUTINES CALLED  DDES, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891024  Changed references from DVNORM to DHVNRM.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DDEABM
!
      INTEGER IALPHA, IBETA, IDELSN, IDID, IFOURU, IG, IHOLD,&
            INFO, IP, IPAR, IPHI, IPSI, ISIG, ITOLD, ITSTAR, ITWOU,&
            IV, IW, IWORK, IWT, IYP, IYPOUT, IYY, LIW, LRW, NEQ
      DOUBLE PRECISION ATOL, RPAR, RTOL, RWORK, T, TOUT, Y
      LOGICAL START,PHASE1,NORND,STIFF,INTOUT
!
      DIMENSION Y(*),INFO(15),RTOL(*),ATOL(*),RWORK(*),IWORK(*),&
                RPAR(*),IPAR(*)
!
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3
!
      EXTERNAL DF
!
!     CHECK FOR AN APPARENT INFINITE LOOP
!
!***FIRST EXECUTABLE STATEMENT  DDEABM
      IF ( INFO(1) == 0 ) IWORK(LIW) = 0
      IF (IWORK(LIW) >= 5) THEN
         IF (T == RWORK(21 + NEQ)) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DDEABM',&
               'AN APPARENT INFINITE LOOP HAS BEEN DETECTED.$$' //&
               'YOU HAVE MADE REPEATED CALLS AT T = ' // XERN3 //&
               ' AND THE INTEGRATION HAS NOT ADVANCED.  CHECK THE ' //&
               'WAY YOU HAVE SET PARAMETERS FOR THE CALL TO THE ' //&
               'CODE, PARTICULARLY INFO(1).', 13, 2)
            RETURN
         ENDIF
      ENDIF
!
!     CHECK LRW AND LIW FOR SUFFICIENT STORAGE ALLOCATION
!
      IDID=0
      IF (LRW < 130+21*NEQ) THEN
         WRITE (XERN1, '(I8)') LRW
         CALL XERMSG ('SLATEC', 'DDEABM', 'THE LENGTH OF THE RWORK ' //&
            'ARRAY MUST BE AT LEAST 130 + 21*NEQ.$$' //&
            'YOU HAVE CALLED THE CODE WITH LRW = ' // XERN1, 1, 1)
         IDID=-33
      ENDIF
!
      IF (LIW < 51) THEN
         WRITE (XERN1, '(I8)') LIW
         CALL XERMSG ('SLATEC', 'DDEABM', 'THE LENGTH OF THE IWORK ' //&
            'ARRAY MUST BE AT LEAST 51.$$YOU HAVE CALLED THE CODE ' //&
            'WITH LIW = ' // XERN1, 2, 1)
         IDID=-33
      ENDIF
!
!     COMPUTE THE INDICES FOR THE ARRAYS TO BE STORED IN THE WORK ARRAY
!
      IYPOUT = 21
      ITSTAR = NEQ + 21
      IYP = 1 + ITSTAR
      IYY = NEQ + IYP
      IWT = NEQ + IYY
      IP = NEQ + IWT
      IPHI = NEQ + IP
      IALPHA = (NEQ*16) + IPHI
      IBETA = 12 + IALPHA
      IPSI = 12 + IBETA
      IV = 12 + IPSI
      IW = 12 + IV
      ISIG = 12 + IW
      IG = 13 + ISIG
      IGI = 13 + IG
      IXOLD = 11 + IGI
      IHOLD = 1 + IXOLD
      ITOLD = 1 + IHOLD
      IDELSN = 1 + ITOLD
      ITWOU = 1 + IDELSN
      IFOURU = 1 + ITWOU
!
      RWORK(ITSTAR) = T
      IF (INFO(1) == 0) GO TO 50
      START = IWORK(21) /= (-1)
      PHASE1 = IWORK(22) /= (-1)
      NORND = IWORK(23) /= (-1)
      STIFF = IWORK(24) /= (-1)
      INTOUT = IWORK(25) /= (-1)
!
 50   CALL DDES(DF,NEQ,T,Y,TOUT,INFO,RTOL,ATOL,IDID,RWORK(IYPOUT),&
               RWORK(IYP),RWORK(IYY),RWORK(IWT),RWORK(IP),RWORK(IPHI),&
               RWORK(IALPHA),RWORK(IBETA),RWORK(IPSI),RWORK(IV),&
               RWORK(IW),RWORK(ISIG),RWORK(IG),RWORK(IGI),RWORK(11),&
               RWORK(12),RWORK(13),RWORK(IXOLD),RWORK(IHOLD),&
               RWORK(ITOLD),RWORK(IDELSN),RWORK(1),RWORK(ITWOU),&
               RWORK(IFOURU),START,PHASE1,NORND,STIFF,INTOUT,IWORK(26),&
               IWORK(27),IWORK(28),IWORK(29),IWORK(30),IWORK(31),&
               IWORK(32),IWORK(33),IWORK(34),IWORK(35),IWORK(45),&
               RPAR,IPAR)
!
      IWORK(21) = -1
      IF (START) IWORK(21) = 1
      IWORK(22) = -1
      IF (PHASE1) IWORK(22) = 1
      IWORK(23) = -1
      IF (NORND) IWORK(23) = 1
      IWORK(24) = -1
      IF (STIFF) IWORK(24) = 1
      IWORK(25) = -1
      IF (INTOUT) IWORK(25) = 1
!
      IF (IDID /= (-2)) IWORK(LIW) = IWORK(LIW) + 1
      IF (T /= RWORK(ITSTAR)) IWORK(LIW) = 0
!
      RETURN
      END SUBROUTINE DDEABM
      
!DECK DDES
      SUBROUTINE DDES (DF, NEQ, T, Y, TOUT, INFO, RTOL, ATOL, IDID,&
         YPOUT, YP, YY, WT, P, PHI, ALPHA, BETA, PSI, V, W, SIG, G, GI,&
         H, EPS, X, XOLD, HOLD, TOLD, DELSGN, TSTOP, TWOU, FOURU, START,&
         PHASE1, NORND, STIFF, INTOUT, NS, KORD, KOLD, INIT, KSTEPS,&
         KLE4, IQUIT, KPREV, IVC, IV, KGI, RPAR, IPAR)
!***BEGIN PROLOGUE  DDES
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DDEABM
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (DES-S, DDES-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   DDEABM merely allocates storage for DDES to relieve the user of the
!   inconvenience of a long call list.  Consequently  DDES  is used as
!   described in the comments for  DDEABM .
!
!***SEE ALSO  DDEABM
!***ROUTINES CALLED  D1MACH, DINTP, DSTEPS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls, cvt GOTOs to
!           IF-THEN-ELSE.  (RWC)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DDES
!
      INTEGER IDID, INFO, INIT, IPAR, IQUIT, IV, IVC, K, KGI, KLE4,&
            KOLD, KORD, KPREV, KSTEPS, L, LTOL, MAXNUM, NATOLP, NEQ,&
            NRTOLP, NS
      DOUBLE PRECISION A, ABSDEL, ALPHA, ATOL, BETA,&
            DEL, DELSGN, DT, EPS, FOURU, G, GI, H,&
            HA, HOLD, P, PHI, PSI, RPAR, RTOL, SIG, T, TOLD, TOUT,&
            TSTOP, TWOU, U, V, W, WT, X, XOLD, Y, YP, YPOUT, YY
      LOGICAL STIFF,CRASH,START,PHASE1,NORND,INTOUT
!
      DIMENSION Y(*),YY(*),WT(*),PHI(NEQ,16),P(*),YP(*),&
        YPOUT(*),PSI(12),ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),&
        GI(11),IV(10),INFO(15),RTOL(*),ATOL(*),RPAR(*),IPAR(*)
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3, XERN4
!
      EXTERNAL DF
!
!.......................................................................
!
!  THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
!  NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE COUNTER
!  IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE
!  WORK.
!
      SAVE MAXNUM
      DATA MAXNUM/500/
!
!.......................................................................
!
!***FIRST EXECUTABLE STATEMENT  DDES
      IF (INFO(1) == 0) THEN
!
! ON THE FIRST CALL , PERFORM INITIALIZATION --
!        DEFINE THE MACHINE UNIT ROUNDOFF QUANTITY  U  BY CALLING THE
!        FUNCTION ROUTINE  D1MACH. THE USER MUST MAKE SURE THAT THE
!        VALUES SET IN D1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
!
         U=d1mach4
!                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS
         TWOU=2.D0*U
         FOURU=4.D0*U
!                       -- SET TERMINATION FLAG
         IQUIT=0
!                       -- SET INITIALIZATION INDICATOR
         INIT=0
!                       -- SET COUNTER FOR ATTEMPTED STEPS
         KSTEPS=0
!                       -- SET INDICATOR FOR INTERMEDIATE-OUTPUT
         INTOUT= .FALSE.
!                       -- SET INDICATOR FOR STIFFNESS DETECTION
         STIFF= .FALSE.
!                       -- SET STEP COUNTER FOR STIFFNESS DETECTION
         KLE4=0
!                       -- SET INDICATORS FOR STEPS CODE
         START= .TRUE.
         PHASE1= .TRUE.
         NORND= .TRUE.
!                       -- RESET INFO(1) FOR SUBSEQUENT CALLS
         INFO(1)=1
      ENDIF
!
!.......................................................................
!
!      CHECK VALIDITY OF INPUT PARAMETERS ON EACH ENTRY
!
      IF (INFO(1) /= 0  .AND.  INFO(1) /= 1) THEN
         WRITE (XERN1, '(I8)') INFO(1)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(1) MUST BE ' //&
            'SET TO 0 FOR THE START OF A NEW PROBLEM, AND MUST BE ' //&
            'SET TO 1 FOLLOWING AN INTERRUPTED TASK.  YOU ARE ' //&
            'ATTEMPTING TO CONTINUE THE INTEGRATION ILLEGALLY BY ' //&
            'CALLING THE CODE WITH INFO(1) = ' // XERN1, 3, 1)
         IDID=-33
      ENDIF
!
      IF (INFO(2) /= 0  .AND.  INFO(2) /= 1) THEN
         WRITE (XERN1, '(I8)') INFO(2)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(2) MUST BE ' //&
            '0 OR 1 INDICATING SCALAR AND VECTOR ERROR TOLERANCES, ' //&
            'RESPECTIVELY.  YOU HAVE CALLED THE CODE WITH INFO(2) = ' //&
            XERN1, 4, 1)
         IDID=-33
      ENDIF
!
      IF (INFO(3) /= 0  .AND.  INFO(3) /= 1) THEN
         WRITE (XERN1, '(I8)') INFO(3)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(3) MUST BE ' //&
            '0 OR 1 INDICATING THE INTERVAL OR INTERMEDIATE-OUTPUT ' //&
            'MODE OF INTEGRATION, RESPECTIVELY.  YOU HAVE CALLED ' //&
            'THE CODE WITH  INFO(3) = ' // XERN1, 5, 1)
         IDID=-33
      ENDIF
!
      IF (INFO(4) /= 0  .AND.  INFO(4) /= 1) THEN
         WRITE (XERN1, '(I8)') INFO(4)
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INFO(4) MUST BE ' //&
            '0 OR 1 INDICATING WHETHER OR NOT THE INTEGRATION ' //&
            'INTERVAL IS TO BE RESTRICTED BY A POINT TSTOP.  YOU ' //&
            'HAVE CALLED THE CODE WITH INFO(4) = ' // XERN1, 14, 1)
         IDID=-33
      ENDIF
!
      IF (NEQ < 1) THEN
         WRITE (XERN1, '(I8)') NEQ
         CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM,  THE NUMBER OF ' //&
            'EQUATIONS NEQ MUST BE A POSITIVE INTEGER.  YOU HAVE ' //&
            'CALLED THE CODE WITH  NEQ = ' // XERN1, 6, 1)
         IDID=-33
      ENDIF
!
      NRTOLP = 0
      NATOLP = 0
      DO 90 K=1,NEQ
         IF (NRTOLP == 0 .AND. RTOL(K) < 0.D0) THEN
            WRITE (XERN1, '(I8)') K
            WRITE (XERN3, '(1PE15.6)') RTOL(K)
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, THE RELATIVE ' //&
               'ERROR TOLERANCES RTOL MUST BE NON-NEGATIVE.  YOU ' //&
               'HAVE CALLED THE CODE WITH  RTOL(' // XERN1 // ') = ' //&
               XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' //&
               'NO FURTHER CHECKING OF RTOL COMPONENTS IS DONE.', 7, 1)
            IDID = -33
            NRTOLP = 1
         ENDIF
!
         IF (NATOLP == 0 .AND. ATOL(K) < 0.D0) THEN
            WRITE (XERN1, '(I8)') K
            WRITE (XERN3, '(1PE15.6)') ATOL(K)
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, THE ABSOLUTE ' //&
               'ERROR TOLERANCES ATOL MUST BE NON-NEGATIVE.  YOU ' //&
               'HAVE CALLED THE CODE WITH  ATOL(' // XERN1 // ') = ' //&
               XERN3 // '.  IN THE CASE OF VECTOR ERROR TOLERANCES, ' //&
               'NO FURTHER CHECKING OF ATOL COMPONENTS IS DONE.', 8, 1)
            IDID = -33
            NATOLP = 1
         ENDIF
!
         IF (INFO(2) == 0) GO TO 100
         IF (NATOLP>0 .AND. NRTOLP>0) GO TO 100
   90 CONTINUE
!
  100 IF (INFO(4) == 1) THEN
         IF (SIGN(1.D0,TOUT-T) /= SIGN(1.D0,TSTOP-T)&
            .OR. ABS(TOUT-T) > ABS(TSTOP-T)) THEN
            WRITE (XERN3, '(1PE15.6)') TOUT
            WRITE (XERN4, '(1PE15.6)') TSTOP
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, YOU HAVE ' //&
               'CALLED THE CODE WITH  TOUT = ' // XERN3 // ' BUT ' //&
               'YOU HAVE ALSO TOLD THE CODE (INFO(4) = 1) NOT TO ' //&
               'INTEGRATE PAST THE POINT TSTOP = ' // XERN4 //&
               ' THESE INSTRUCTIONS CONFLICT.', 14, 1)
            IDID=-33
         ENDIF
      ENDIF
!
!     CHECK SOME CONTINUATION POSSIBILITIES
!
      IF (INIT /= 0) THEN
         IF (T == TOUT) THEN
            WRITE (XERN3, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, YOU HAVE ' //&
               'CALLED THE CODE WITH  T = TOUT = ' // XERN3 //&
               '$$THIS IS NOT ALLOWED ON CONTINUATION CALLS.', 9, 1)
            IDID=-33
         ENDIF
!
         IF (T /= TOLD) THEN
            WRITE (XERN3, '(1PE15.6)') TOLD
            WRITE (XERN4, '(1PE15.6)') T
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, YOU HAVE ' //&
               'CHANGED THE VALUE OF T FROM ' // XERN3 // ' TO ' //&
               XERN4 //'  THIS IS NOT ALLOWED ON CONTINUATION CALLS.',&
               10, 1)
            IDID=-33
         ENDIF
!
         IF (INIT /= 1) THEN
            IF (DELSGN*(TOUT-T) < 0.D0) THEN
               WRITE (XERN3, '(1PE15.6)') TOUT
               CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, BY ' //&
                  'CALLING THE CODE WITH TOUT = ' // XERN3 //&
                  ' YOU ARE ATTEMPTING TO CHANGE THE DIRECTION OF ' //&
                  'INTEGRATION.$$THIS IS NOT ALLOWED WITHOUT ' //&
                  'RESTARTING.', 11, 1)
               IDID=-33
            ENDIF
         ENDIF
      ENDIF
!
!     INVALID INPUT DETECTED
!
      IF (IDID == (-33)) THEN
         IF (IQUIT /= (-33)) THEN
            IQUIT = -33
            INFO(1) = -1
         ELSE
            CALL XERMSG ('SLATEC', 'DDES', 'IN DDEABM, INVALID ' //&
               'INPUT WAS DETECTED ON SUCCESSIVE ENTRIES.  IT IS ' //&
               'IMPOSSIBLE TO PROCEED BECAUSE YOU HAVE NOT ' //&
               'CORRECTED THE PROBLEM, SO EXECUTION IS BEING ' //&
               'TERMINATED.', 12, 2)
         ENDIF
         RETURN
      ENDIF
!
!.......................................................................
!
!     RTOL = ATOL = 0. IS ALLOWED AS VALID INPUT AND INTERPRETED AS
!     ASKING FOR THE MOST ACCURATE SOLUTION POSSIBLE. IN THIS CASE,
!     THE RELATIVE ERROR TOLERANCE RTOL IS RESET TO THE SMALLEST VALUE
!     FOURU WHICH IS LIKELY TO BE REASONABLE FOR THIS METHOD AND MACHINE
!
      DO 180 K=1,NEQ
        IF (RTOL(K)+ATOL(K) > 0.D0) GO TO 170
        RTOL(K)=FOURU
        IDID=-2
  170   IF (INFO(2) == 0) GO TO 190
  180   CONTINUE
!
  190 IF (IDID /= (-2)) GO TO 200
!                       RTOL=ATOL=0 ON INPUT, SO RTOL IS CHANGED TO A
!                                                SMALL POSITIVE VALUE
      INFO(1)=-1
      RETURN
!
!     BRANCH ON STATUS OF INITIALIZATION INDICATOR
!            INIT=0 MEANS INITIAL DERIVATIVES AND NOMINAL STEP SIZE
!                   AND DIRECTION NOT YET SET
!            INIT=1 MEANS NOMINAL STEP SIZE AND DIRECTION NOT YET SET
!            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED
!
  200 IF (INIT == 0) GO TO 210
      IF (INIT == 1) GO TO 220
      GO TO 240
!
!.......................................................................
!
!     MORE INITIALIZATION --
!                         -- EVALUATE INITIAL DERIVATIVES
!
  210 INIT=1
      A=T
      CALL DF(A,Y,YP,RPAR,IPAR)
      IF (T /= TOUT) GO TO 220
      IDID=2
      DO 215 L = 1,NEQ
  215    YPOUT(L) = YP(L)
      TOLD=T
      RETURN
!
!                         -- SET INDEPENDENT AND DEPENDENT VARIABLES
!                                              X AND YY(*) FOR STEPS
!                         -- SET SIGN OF INTEGRATION DIRECTION
!                         -- INITIALIZE THE STEP SIZE
!
  220 INIT = 2
      X = T
      DO 230 L = 1,NEQ
  230   YY(L) = Y(L)
      DELSGN = SIGN(1.0D0,TOUT-T)
      H = SIGN(MAX(FOURU*ABS(X),ABS(TOUT-X)),TOUT-X)
!
!.......................................................................
!
!   ON EACH CALL SET INFORMATION WHICH DETERMINES THE ALLOWED INTERVAL
!   OF INTEGRATION BEFORE RETURNING WITH AN ANSWER AT TOUT
!
  240 DEL = TOUT - T
      ABSDEL = ABS(DEL)
!
!.......................................................................
!
!   IF ALREADY PAST OUTPUT POINT, INTERPOLATE AND RETURN
!
  250 IF(ABS(X-T) < ABSDEL) GO TO 260
      CALL DINTP(X,YY,TOUT,Y,YPOUT,NEQ,KOLD,PHI,IVC,IV,KGI,GI,&
                                              ALPHA,G,W,XOLD,P)
      IDID = 3
      IF (X /= TOUT) GO TO 255
      IDID = 2
      INTOUT = .FALSE.
  255 T = TOUT
      TOLD = T
      RETURN
!
!   IF CANNOT GO PAST TSTOP AND SUFFICIENTLY CLOSE,
!   EXTRAPOLATE AND RETURN
!
  260 IF (INFO(4) /= 1) GO TO 280
      IF (ABS(TSTOP-X) >= FOURU*ABS(X)) GO TO 280
      DT = TOUT - X
      DO 270 L = 1,NEQ
  270   Y(L) = YY(L) + DT*YP(L)
      CALL DF(TOUT,Y,YPOUT,RPAR,IPAR)
      IDID = 3
      T = TOUT
      TOLD = T
      RETURN
!
  280 IF (INFO(3) == 0  .OR.  .NOT.INTOUT) GO TO 300
!
!   INTERMEDIATE-OUTPUT MODE
!
      IDID = 1
      DO 290 L = 1,NEQ
        Y(L)=YY(L)
  290   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INTOUT = .FALSE.
      RETURN
!
!.......................................................................
!
!     MONITOR NUMBER OF STEPS ATTEMPTED
!
  300 IF (KSTEPS <= MAXNUM) GO TO 330
!
!                       A SIGNIFICANT AMOUNT OF WORK HAS BEEN EXPENDED
      IDID=-1
      KSTEPS=0
      IF (.NOT. STIFF) GO TO 310
!
!                       PROBLEM APPEARS TO BE STIFF
      IDID=-4
      STIFF= .FALSE.
      KLE4=0
!
  310 DO 320 L = 1,NEQ
        Y(L) = YY(L)
  320   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INFO(1) = -1
      INTOUT = .FALSE.
      RETURN
!
!.......................................................................
!
!   LIMIT STEP SIZE, SET WEIGHT VECTOR AND TAKE A STEP
!
  330 HA = ABS(H)
      IF (INFO(4) /= 1) GO TO 340
      HA = MIN(HA,ABS(TSTOP-X))
  340 H = SIGN(HA,H)
      EPS = 1.0D0
      LTOL = 1
      DO 350 L = 1,NEQ
        IF (INFO(2) == 1) LTOL = L
        WT(L) = RTOL(LTOL)*ABS(YY(L)) + ATOL(LTOL)
        IF (WT(L) <= 0.0D0) GO TO 360
  350   CONTINUE
      GO TO 380
!
!                       RELATIVE ERROR CRITERION INAPPROPRIATE
  360 IDID = -3
      DO 370 L = 1,NEQ
        Y(L) = YY(L)
  370   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INFO(1) = -1
      INTOUT = .FALSE.
      RETURN
!
  380 CALL DSTEPS(DF,NEQ,YY,X,H,EPS,WT,START,HOLD,KORD,KOLD,CRASH,PHI,P,&
                 YP,PSI,ALPHA,BETA,SIG,V,W,G,PHASE1,NS,NORND,KSTEPS,&
                 TWOU,FOURU,XOLD,KPREV,IVC,IV,KGI,GI,RPAR,IPAR)
!
!.......................................................................
!
      IF(.NOT.CRASH) GO TO 420
!
!                       TOLERANCES TOO SMALL
      IDID = -2
      RTOL(1) = EPS*RTOL(1)
      ATOL(1) = EPS*ATOL(1)
      IF (INFO(2) == 0) GO TO 400
      DO 390 L = 2,NEQ
        RTOL(L) = EPS*RTOL(L)
  390   ATOL(L) = EPS*ATOL(L)
  400 DO 410 L = 1,NEQ
        Y(L) = YY(L)
  410   YPOUT(L) = YP(L)
      T = X
      TOLD = T
      INFO(1) = -1
      INTOUT = .FALSE.
      RETURN
!
!   (STIFFNESS TEST) COUNT NUMBER OF CONSECUTIVE STEPS TAKEN WITH THE
!   ORDER OF THE METHOD BEING LESS OR EQUAL TO FOUR
!
  420 KLE4 = KLE4 + 1
      IF(KOLD > 4) KLE4 = 0
      IF(KLE4 >= 50) STIFF = .TRUE.
      INTOUT = .TRUE.
      GO TO 250
      END SUBROUTINE DDES
      
!DECK DHSTRT
      SUBROUTINE DHSTRT (DF, NEQ, A, B, Y, YPRIME, ETOL, MORDER, SMALL,&
         BIG, SPY, PV, YP, SF, RPAR, IPAR, H)
!***BEGIN PROLOGUE  DHSTRT
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DDEABM, DDEBDF and DDERKF
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (HSTART-S, DHSTRT-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   DHSTRT computes a starting step size to be used in solving initial
!   value problems in ordinary differential equations.
!
! **********************************************************************
!  ABSTRACT
!
!     Subroutine DHSTRT computes a starting step size to be used by an
!     initial value method in solving ordinary differential equations.
!     It is based on an estimate of the local Lipschitz constant for the
!     differential equation   (lower bound on a norm of the Jacobian) ,
!     a bound on the differential equation  (first derivative) , and
!     a bound on the partial derivative of the equation with respect to
!     the independent variable.
!     (all approximated near the initial point A)
!
!     Subroutine DHSTRT uses a function subprogram DHVNRM for computing
!     a vector norm. The maximum norm is presently utilized though it
!     can easily be replaced by any other vector norm. It is presumed
!     that any replacement norm routine would be carefully coded to
!     prevent unnecessary underflows or overflows from occurring, and
!     also, would not alter the vector or number of components.
!
! **********************************************************************
!  On input you must provide the following
!
!      DF -- This is a subroutine of the form
!                               DF(X,U,UPRIME,RPAR,IPAR)
!             which defines the system of first order differential
!             equations to be solved. For the given values of X and the
!             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
!             evaluate the NEQ components of the system of differential
!             equations  DU/DX=DF(X,U)  and store the derivatives in the
!             array UPRIME(*), that is,  UPRIME(I) = * DU(I)/DX *  for
!             equations I=1,...,NEQ.
!
!             Subroutine DF must not alter X or U(*). You must declare
!             the name DF in an external statement in your program that
!             calls DHSTRT. You must dimension U and UPRIME in DF.
!
!             RPAR and IPAR are DOUBLE PRECISION and INTEGER parameter
!             arrays which you can use for communication between your
!             program and subroutine DF. They are not used or altered by
!             DHSTRT. If you do not need RPAR or IPAR, ignore these
!             parameters by treating them as dummy arguments. If you do
!             choose to use them, dimension them in your program and in
!             DF as arrays of appropriate length.
!
!      NEQ -- This is the number of (first order) differential equations
!             to be integrated.
!
!      A -- This is the initial point of integration.
!
!      B -- This is a value of the independent variable used to define
!             the direction of integration. A reasonable choice is to
!             set  B  to the first point at which a solution is desired.
!             You can also use  B, if necessary, to restrict the length
!             of the first integration step because the algorithm will
!             not compute a starting step length which is bigger than
!             ABS(B-A), unless  B  has been chosen too close to  A.
!             (it is presumed that DHSTRT has been called with  B
!             different from  A  on the machine being used. Also see the
!             discussion about the parameter  SMALL.)
!
!      Y(*) -- This is the vector of initial values of the NEQ solution
!             components at the initial point  A.
!
!      YPRIME(*) -- This is the vector of derivatives of the NEQ
!             solution components at the initial point  A.
!             (defined by the differential equations in subroutine DF)
!
!      ETOL -- This is the vector of error tolerances corresponding to
!             the NEQ solution components. It is assumed that all
!             elements are positive. Following the first integration
!             step, the tolerances are expected to be used by the
!             integrator in an error test which roughly requires that
!                        ABS(LOCAL ERROR)  <=  ETOL
!             for each vector component.
!
!      MORDER -- This is the order of the formula which will be used by
!             the initial value method for taking the first integration
!             step.
!
!      SMALL -- This is a small positive machine dependent constant
!             which is used for protecting against computations with
!             numbers which are too small relative to the precision of
!             floating point arithmetic.  SMALL  should be set to
!             (approximately) the smallest positive DOUBLE PRECISION
!             number such that  (1.+SMALL) > 1.  on the machine being
!             used. The quantity  SMALL**(3/8)  is used in computing
!             increments of variables for approximating derivatives by
!             differences.  Also the algorithm will not compute a
!             starting step length which is smaller than
!             100*SMALL*ABS(A).
!
!      BIG -- This is a large positive machine dependent constant which
!             is used for preventing machine overflows. A reasonable
!             choice is to set big to (approximately) the square root of
!             the largest DOUBLE PRECISION number which can be held in
!             the machine.
!
!      SPY(*),PV(*),YP(*),SF(*) -- These are DOUBLE PRECISION work
!             arrays of length NEQ which provide the routine with needed
!             storage space.
!
!      RPAR,IPAR -- These are parameter arrays, of DOUBLE PRECISION and
!             INTEGER type, respectively, which can be used for
!             communication between your program and the DF subroutine.
!             They are not used or altered by DHSTRT.
!
! **********************************************************************
!  On Output  (after the return from DHSTRT),
!
!      H -- is an appropriate starting step size to be attempted by the
!             differential equation method.
!
!           All parameters in the call list remain unchanged except for
!           the working arrays SPY(*),PV(*),YP(*), and SF(*).
!
! **********************************************************************
!
!***SEE ALSO  DDEABM, DDEBDF, DDERKF
!***ROUTINES CALLED  DHVNRM
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891024  Changed references from DVNORM to DHVNRM.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DHSTRT
!
      INTEGER IPAR, J, K, LK, MORDER, NEQ
      DOUBLE PRECISION A, ABSDX, B, BIG, DA, DELF, DELY,&
            DFDUB, DFDXB,&
            DX, DY, ETOL, FBND, H, PV, RELPER, RPAR, SF, SMALL, SPY,&
            SRYDPB, TOLEXP, TOLMIN, TOLP, TOLSUM, Y, YDPB, YP, YPRIME
      DIMENSION Y(*),YPRIME(*),ETOL(*),SPY(*),PV(*),YP(*),&
                SF(*),RPAR(*),IPAR(*)
      EXTERNAL DF
!
!     ..................................................................
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 160
!***FIRST EXECUTABLE STATEMENT  DHSTRT
         DX = B - A
         ABSDX = ABS(DX)
         RELPER = SMALL**0.375D0
!
!        ...............................................................
!
!             COMPUTE AN APPROXIMATE BOUND (DFDXB) ON THE PARTIAL
!             DERIVATIVE OF THE EQUATION WITH RESPECT TO THE
!             INDEPENDENT VARIABLE. PROTECT AGAINST AN OVERFLOW.
!             ALSO COMPUTE A BOUND (FBND) ON THE FIRST DERIVATIVE
!             LOCALLY.
!
         DA = SIGN(MAX(MIN(RELPER*ABS(A),ABSDX),&
                          100.0D0*SMALL*ABS(A)),DX)
         IF (DA == 0.0D0) DA = RELPER*DX
         CALL DF(A+DA,Y,SF,RPAR,IPAR)
         DO 10 J = 1, NEQ
            YP(J) = SF(J) - YPRIME(J)
   10    CONTINUE
         DELF = DHVNRM(YP,NEQ)
         DFDXB = BIG
         IF (DELF < BIG*ABS(DA)) DFDXB = DELF/ABS(DA)
         FBND = DHVNRM(SF,NEQ)
!
!        ...............................................................
!
!             COMPUTE AN ESTIMATE (DFDUB) OF THE LOCAL LIPSCHITZ
!             CONSTANT FOR THE SYSTEM OF DIFFERENTIAL EQUATIONS. THIS
!             ALSO REPRESENTS AN ESTIMATE OF THE NORM OF THE JACOBIAN
!             LOCALLY.  THREE ITERATIONS (TWO WHEN NEQ=1) ARE USED TO
!             ESTIMATE THE LIPSCHITZ CONSTANT BY NUMERICAL DIFFERENCES.
!             THE FIRST PERTURBATION VECTOR IS BASED ON THE INITIAL
!             DERIVATIVES AND DIRECTION OF INTEGRATION. THE SECOND
!             PERTURBATION VECTOR IS FORMED USING ANOTHER EVALUATION OF
!             THE DIFFERENTIAL EQUATION.  THE THIRD PERTURBATION VECTOR
!             IS FORMED USING PERTURBATIONS BASED ONLY ON THE INITIAL
!             VALUES. COMPONENTS THAT ARE ZERO ARE ALWAYS CHANGED TO
!             NON-ZERO VALUES (EXCEPT ON THE FIRST ITERATION). WHEN
!             INFORMATION IS AVAILABLE, CARE IS TAKEN TO ENSURE THAT
!             COMPONENTS OF THE PERTURBATION VECTOR HAVE SIGNS WHICH ARE
!             CONSISTENT WITH THE SLOPES OF LOCAL SOLUTION CURVES.
!             ALSO CHOOSE THE LARGEST BOUND (FBND) FOR THE FIRST
!             DERIVATIVE.
!
!                               PERTURBATION VECTOR SIZE IS HELD
!                               CONSTANT FOR ALL ITERATIONS. COMPUTE
!                               THIS CHANGE FROM THE
!                                       SIZE OF THE VECTOR OF INITIAL
!                                       VALUES.
         DELY = RELPER*DHVNRM(Y,NEQ)
         IF (DELY == 0.0D0) DELY = RELPER
         DELY = SIGN(DELY,DX)
         DELF = DHVNRM(YPRIME,NEQ)
         FBND = MAX(FBND,DELF)
         IF (DELF == 0.0D0) GO TO 30
!           USE INITIAL DERIVATIVES FOR FIRST PERTURBATION
            DO 20 J = 1, NEQ
               SPY(J) = YPRIME(J)
               YP(J) = YPRIME(J)
   20       CONTINUE
         GO TO 50
   30    CONTINUE
!           CANNOT HAVE A NULL PERTURBATION VECTOR
            DO 40 J = 1, NEQ
               SPY(J) = 0.0D0
               YP(J) = 1.0D0
   40       CONTINUE
            DELF = DHVNRM(YP,NEQ)
   50    CONTINUE
!
         DFDUB = 0.0D0
         LK = MIN(NEQ+1,3)
         DO 140 K = 1, LK
!           DEFINE PERTURBED VECTOR OF INITIAL VALUES
            DO 60 J = 1, NEQ
               PV(J) = Y(J) + DELY*(YP(J)/DELF)
   60       CONTINUE
            IF (K == 2) GO TO 80
!              EVALUATE DERIVATIVES ASSOCIATED WITH PERTURBED
!              VECTOR  AND  COMPUTE CORRESPONDING DIFFERENCES
               CALL DF(A,PV,YP,RPAR,IPAR)
               DO 70 J = 1, NEQ
                  PV(J) = YP(J) - YPRIME(J)
   70          CONTINUE
            GO TO 100
   80       CONTINUE
!              USE A SHIFTED VALUE OF THE INDEPENDENT VARIABLE
!                                    IN COMPUTING ONE ESTIMATE
               CALL DF(A+DA,PV,YP,RPAR,IPAR)
               DO 90 J = 1, NEQ
                  PV(J) = YP(J) - SF(J)
   90          CONTINUE
  100       CONTINUE
!           CHOOSE LARGEST BOUNDS ON THE FIRST DERIVATIVE
!                          AND A LOCAL LIPSCHITZ CONSTANT
            FBND = MAX(FBND,DHVNRM(YP,NEQ))
            DELF = DHVNRM(PV,NEQ)
!        ...EXIT
            IF (DELF >= BIG*ABS(DELY)) GO TO 150
            DFDUB = MAX(DFDUB,DELF/ABS(DELY))
!     ......EXIT
            IF (K == LK) GO TO 160
!           CHOOSE NEXT PERTURBATION VECTOR
            IF (DELF == 0.0D0) DELF = 1.0D0
            DO 130 J = 1, NEQ
               IF (K == 2) GO TO 110
                  DY = ABS(PV(J))
                  IF (DY == 0.0D0) DY = DELF
               GO TO 120
  110          CONTINUE
                  DY = Y(J)
                  IF (DY == 0.0D0) DY = DELY/RELPER
  120          CONTINUE
               IF (SPY(J) == 0.0D0) SPY(J) = YP(J)
               IF (SPY(J) /= 0.0D0) DY = SIGN(DY,SPY(J))
               YP(J) = DY
  130       CONTINUE
            DELF = DHVNRM(YP,NEQ)
  140    CONTINUE
  150    CONTINUE
!
!        PROTECT AGAINST AN OVERFLOW
         DFDUB = BIG
  160 CONTINUE
!
!     ..................................................................
!
!          COMPUTE A BOUND (YDPB) ON THE NORM OF THE SECOND DERIVATIVE
!
      YDPB = DFDXB + DFDUB*FBND
!
!     ..................................................................
!
!          DEFINE THE TOLERANCE PARAMETER UPON WHICH THE STARTING STEP
!          SIZE IS TO BE BASED.  A VALUE IN THE MIDDLE OF THE ERROR
!          TOLERANCE RANGE IS SELECTED.
!
      TOLMIN = BIG
      TOLSUM = 0.0D0
      DO 170 K = 1, NEQ
         TOLEXP = LOG10(ETOL(K))
         TOLMIN = MIN(TOLMIN,TOLEXP)
         TOLSUM = TOLSUM + TOLEXP
  170 CONTINUE
      TOLP = 10.0D0**(0.5D0*(TOLSUM/NEQ + TOLMIN)/(MORDER+1))
!
!     ..................................................................
!
!          COMPUTE A STARTING STEP SIZE BASED ON THE ABOVE FIRST AND
!          SECOND DERIVATIVE INFORMATION
!
!                            RESTRICT THE STEP LENGTH TO BE NOT BIGGER
!                            THAN ABS(B-A).   (UNLESS  B  IS TOO CLOSE
!                            TO  A)
      H = ABSDX
!
      IF (YDPB /= 0.0D0 .OR. FBND /= 0.0D0) GO TO 180
!
!        BOTH FIRST DERIVATIVE TERM (FBND) AND SECOND
!                     DERIVATIVE TERM (YDPB) ARE ZERO
         IF (TOLP < 1.0D0) H = ABSDX*TOLP
      GO TO 200
  180 CONTINUE
!
      IF (YDPB /= 0.0D0) GO TO 190
!
!        ONLY SECOND DERIVATIVE TERM (YDPB) IS ZERO
         IF (TOLP < FBND*ABSDX) H = TOLP/FBND
      GO TO 200
  190 CONTINUE
!
!        SECOND DERIVATIVE TERM (YDPB) IS NON-ZERO
         SRYDPB = SQRT(0.5D0*YDPB)
         IF (TOLP < SRYDPB*ABSDX) H = TOLP/SRYDPB
  200 CONTINUE
!
!     FURTHER RESTRICT THE STEP LENGTH TO BE NOT
!                               BIGGER THAN  1/DFDUB
      IF (H*DFDUB > 1.0D0) H = 1.0D0/DFDUB
!
!     FINALLY, RESTRICT THE STEP LENGTH TO BE NOT
!     SMALLER THAN  100*SMALL*ABS(A).  HOWEVER, IF
!     A=0. AND THE COMPUTED H UNDERFLOWED TO ZERO,
!     THE ALGORITHM RETURNS  SMALL*ABS(B)  FOR THE
!                                     STEP LENGTH.
      H = MAX(H,100.0D0*SMALL*ABS(A))
      IF (H == 0.0D0) H = SMALL*ABS(B)
!
!     NOW SET DIRECTION OF INTEGRATION
      H = SIGN(H,DX)
!
      RETURN
      END SUBROUTINE DHSTRT
      
!DECK DHVNRM
      DOUBLE PRECISION FUNCTION DHVNRM (V, NCOMP)
!***BEGIN PROLOGUE  DHVNRM
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DDEABM, DDEBDF and DDERKF
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (HVNRM-S, DHVNRM-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     Compute the maximum norm of the vector V(*) of length NCOMP and
!     return the result as DHVNRM
!
!***SEE ALSO  DDEABM, DDEBDF, DDERKF
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891024  Changed references from DVNORM to DHVNRM.  (WRB)
!   891024  Changed routine name from DVNORM to DHVNRM.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DHVNRM
!
      INTEGER K, NCOMP
      DOUBLE PRECISION V
      DIMENSION V(*)
!***FIRST EXECUTABLE STATEMENT  DHVNRM
      DHVNRM = 0.0D0
      DO 10 K = 1, NCOMP
         DHVNRM = MAX(DHVNRM,ABS(V(K)))
   10 CONTINUE
      RETURN
      END FUNCTION DHVNRM
      
!DECK DINTP
      SUBROUTINE DINTP (X, Y, XOUT, YOUT, YPOUT, NEQN, KOLD, PHI, IVC,&
         IV, KGI, GI, ALPHA, OG, OW, OX, OY)
!***BEGIN PROLOGUE  DINTP
!***PURPOSE  Approximate the solution at XOUT by evaluating the
!            polynomial computed in DSTEPS at XOUT.  Must be used in
!            conjunction with DSTEPS.
!***LIBRARY   SLATEC (DEPAC)
!***CATEGORY  I1A1B
!***TYPE      DOUBLE PRECISION (SINTRP-S, DINTP-D)
!***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE,
!             ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR,
!             SMOOTH INTERPOLANT
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   The methods in subroutine  DSTEPS  approximate the solution near  X
!   by a polynomial.  Subroutine  DINTP  approximates the solution at
!   XOUT  by evaluating the polynomial there.  Information defining this
!   polynomial is passed from  DSTEPS  so  DINTP  cannot be used alone.
!
!   Subroutine DSTEPS is completely explained and documented in the text
!   "Computer Solution of Ordinary Differential Equations, the Initial
!   Value Problem"  by L. F. Shampine and M. K. Gordon.
!
!   Input to DINTP --
!
!   The user provides storage in the calling program for the arrays in
!   the call list
!      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),OY(NEQN)
!                AND ALPHA(12),OG(13),OW(12),GI(11),IV(10)
!   and defines
!      XOUT -- point at which solution is desired.
!   The remaining parameters are defined in  DSTEPS  and passed to
!   DINTP  from that subroutine
!
!   Output from  DINTP --
!
!      YOUT(*) -- solution at  XOUT
!      YPOUT(*) -- derivative of solution at  XOUT
!   The remaining parameters are returned unaltered from their input
!   values.  Integration with  DSTEPS  may be continued.
!
!***REFERENCES  H. A. Watts, A smoother interpolant for DE/STEP, INTRP
!                 II, Report SAND84-0293, Sandia Laboratories, 1984.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   840201  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DINTP
!
      INTEGER I, IQ, IV, IVC, IW, J, JQ, KGI, KOLD, KP1, KP2,&
              L, M, NEQN
      DOUBLE PRECISION ALP, ALPHA, C, G, GDI, GDIF, GI, GAMMA, H, HI,&
             HMU, OG, OW, OX, OY, PHI, RMU, SIGMA, TEMP1, TEMP2, TEMP3,&
             W, X, XI, XIM1, XIQ, XOUT, Y, YOUT, YPOUT
!
      DIMENSION Y(*),YOUT(*),YPOUT(*),PHI(NEQN,16),OY(*)
      DIMENSION G(13),C(13),W(13),OG(13),OW(12),ALPHA(12),GI(11),IV(10)
!
!***FIRST EXECUTABLE STATEMENT  DINTP
      KP1 = KOLD + 1
      KP2 = KOLD + 2
!
      HI = XOUT - OX
      H = X - OX
      XI = HI/H
      XIM1 = XI - 1.D0
!
!   INITIALIZE W(*) FOR COMPUTING G(*)
!
      XIQ = XI
      DO 10 IQ = 1,KP1
        XIQ = XI*XIQ
        TEMP1 = IQ*(IQ+1)
 10     W(IQ) = XIQ/TEMP1
!
!   COMPUTE THE DOUBLE INTEGRAL TERM GDI
!
      IF (KOLD <= KGI) GO TO 50
      IF (IVC > 0) GO TO 20
      GDI = 1.0D0/TEMP1
      M = 2
      GO TO 30
 20   IW = IV(IVC)
      GDI = OW(IW)
      M = KOLD - IW + 3
 30   IF (M > KOLD) GO TO 60
      DO 40 I = M,KOLD
 40     GDI = OW(KP2-I) - ALPHA(I)*GDI
      GO TO 60
 50   GDI = GI(KOLD)
!
!   COMPUTE G(*) AND C(*)
!
 60   G(1) = XI
      G(2) = 0.5D0*XI*XI
      C(1) = 1.0D0
      C(2) = XI
      IF (KOLD < 2) GO TO 90
      DO 80 I = 2,KOLD
        ALP = ALPHA(I)
        GAMMA = 1.0D0 + XIM1*ALP
        L = KP2 - I
        DO 70 JQ = 1,L
 70       W(JQ) = GAMMA*W(JQ) - ALP*W(JQ+1)
        G(I+1) = W(1)
 80     C(I+1) = GAMMA*C(I)
!
!   DEFINE INTERPOLATION PARAMETERS
!
 90   SIGMA = (W(2) - XIM1*W(1))/GDI
      RMU = XIM1*C(KP1)/GDI
      HMU = RMU/H
!
!   INTERPOLATE FOR THE SOLUTION -- YOUT
!   AND FOR THE DERIVATIVE OF THE SOLUTION -- YPOUT
!
      DO 100 L = 1,NEQN
        YOUT(L) = 0.0D0
 100    YPOUT(L) = 0.0D0
      DO 120 J = 1,KOLD
        I = KP2 - J
        GDIF = OG(I) - OG(I-1)
        TEMP2 = (G(I) - G(I-1)) - SIGMA*GDIF
        TEMP3 = (C(I) - C(I-1)) + RMU*GDIF
        DO 110 L = 1,NEQN
          YOUT(L) = YOUT(L) + TEMP2*PHI(L,I)
 110      YPOUT(L) = YPOUT(L) + TEMP3*PHI(L,I)
 120    CONTINUE
      DO 130 L = 1,NEQN
        YOUT(L) = ((1.0D0 - SIGMA)*OY(L) + SIGMA*Y(L)) +&
                   H*(YOUT(L) + (G(1) - SIGMA*OG(1))*PHI(L,1))
 130    YPOUT(L) = HMU*(OY(L) - Y(L)) +&
                      (YPOUT(L) + (C(1) + RMU*OG(1))*PHI(L,1))
!
      RETURN
      END SUBROUTINE DINTP
      
!DECK DSTEPS
      SUBROUTINE DSTEPS (DF, NEQN, Y, X, H, EPS, WT, START, HOLD, K,&
         KOLD, CRASH, PHI, P, YP, PSI, ALPHA, BETA, SIG, V, W, G,&
         PHASE1, NS, NORND, KSTEPS, TWOU, FOURU, XOLD, KPREV, IVC, IV,&
         KGI, GI, RPAR, IPAR)
!***BEGIN PROLOGUE  DSTEPS
!***PURPOSE  Integrate a system of first order ordinary differential
!            equations one step.
!***LIBRARY   SLATEC (DEPAC)
!***CATEGORY  I1A1B
!***TYPE      DOUBLE PRECISION (STEPS-S, DSTEPS-D)
!***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE,
!             ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR
!***AUTHOR  Shampine, L. F., (SNLA)
!           Gordon, M. K., (SNLA)
!             MODIFIED BY H.A. WATTS
!***DESCRIPTION
!
!   Written by L. F. Shampine and M. K. Gordon
!
!   Abstract
!
!   Subroutine  DSTEPS  is normally used indirectly through subroutine
!   DDEABM .  Because  DDEABM  suffices for most problems and is much
!   easier to use, using it should be considered before using  DSTEPS
!   alone.
!
!   Subroutine DSTEPS integrates a system of  NEQN  first order ordinary
!   differential equations one step, normally from X to X+H, using a
!   modified divided difference form of the Adams Pece formulas.  Local
!   extrapolation is used to improve absolute stability and accuracy.
!   The code adjusts its order and step size to control the local error
!   per unit step in a generalized sense.  Special devices are included
!   to control roundoff error and to detect when the user is requesting
!   too much accuracy.
!
!   This code is completely explained and documented in the text,
!   Computer Solution of Ordinary Differential Equations, The Initial
!   Value Problem  by L. F. Shampine and M. K. Gordon.
!   Further details on use of this code are available in "Solving
!   Ordinary Differential Equations with ODE, STEP, and INTRP",
!   by L. F. Shampine and M. K. Gordon, SLA-73-1060.
!
!
!   The parameters represent --
!      DF -- subroutine to evaluate derivatives
!      NEQN -- number of equations to be integrated
!      Y(*) -- solution vector at X
!      X -- independent variable
!      H -- appropriate step size for next step.  Normally determined by
!           code
!      EPS -- local error tolerance
!      WT(*) -- vector of weights for error criterion
!      START -- logical variable set .TRUE. for first step,  .FALSE.
!           otherwise
!      HOLD -- step size used for last successful step
!      K -- appropriate order for next step (determined by code)
!      KOLD -- order used for last successful step
!      CRASH -- logical variable set .TRUE. when no step can be taken,
!           .FALSE. otherwise.
!      YP(*) -- derivative of solution vector at  X  after successful
!           step
!      KSTEPS -- counter on attempted steps
!      TWOU -- 2.*U where U is machine unit roundoff quantity
!      FOURU -- 4.*U where U is machine unit roundoff quantity
!      RPAR,IPAR -- parameter arrays which you may choose to use
!            for communication between your program and subroutine F.
!            They are not altered or used by DSTEPS.
!   The variables X,XOLD,KOLD,KGI and IVC and the arrays Y,PHI,ALPHA,G,
!   W,P,IV and GI are required for the interpolation subroutine SINTRP.
!   The remaining variables and arrays are included in the call list
!   only to eliminate local retention of variables between calls.
!
!   Input to DSTEPS
!
!      First call --
!
!   The user must provide storage in his calling program for all arrays
!   in the call list, namely
!
!     DIMENSION Y(NEQN),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN),PSI(12),
!    1  ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10),
!    2  RPAR(*),IPAR(*)
!
!    **Note**
!
!   The user must also declare  START ,  CRASH ,  PHASE1  and  NORND
!   logical variables and  DF  an EXTERNAL subroutine, supply the
!   subroutine  DF(X,Y,YP)  to evaluate
!      DY(I)/DX = YP(I) = DF(X,Y(1),Y(2),...,Y(NEQN))
!   and initialize only the following parameters.
!      NEQN -- number of equations to be integrated
!      Y(*) -- vector of initial values of dependent variables
!      X -- initial value of the independent variable
!      H -- nominal step size indicating direction of integration
!           and maximum size of step.  Must be variable
!      EPS -- local error tolerance per step.  Must be variable
!      WT(*) -- vector of non-zero weights for error criterion
!      START -- .TRUE.
!      YP(*) -- vector of initial derivative values
!      KSTEPS -- set KSTEPS to zero
!      TWOU -- 2.*U where U is machine unit roundoff quantity
!      FOURU -- 4.*U where U is machine unit roundoff quantity
!   Define U to be the machine unit roundoff quantity by calling
!   the function routine  D1MACH,  U = D1MACH(4), or by
!   computing U so that U is the smallest positive number such
!   that 1.0+U > 1.0.
!
!   DSTEPS  requires that the L2 norm of the vector with components
!   LOCAL ERROR(L)/WT(L)  be less than  EPS  for a successful step.  The
!   array  WT  allows the user to specify an error test appropriate
!   for his problem.  For example,
!      WT(L) = 1.0  specifies absolute error,
!            = ABS(Y(L))  error relative to the most recent value of the
!                 L-th component of the solution,
!            = ABS(YP(L))  error relative to the most recent value of
!                 the L-th component of the derivative,
!            = MAX(WT(L),ABS(Y(L)))  error relative to the largest
!                 magnitude of L-th component obtained so far,
!            = ABS(Y(L))*RELERR/EPS + ABSERR/EPS  specifies a mixed
!                 relative-absolute test where  RELERR  is relative
!                 error,  ABSERR  is absolute error and  EPS =
!                 MAX(RELERR,ABSERR) .
!
!      Subsequent calls --
!
!   Subroutine  DSTEPS  is designed so that all information needed to
!   continue the integration, including the step size  H  and the order
!   K , is returned with each step.  With the exception of the step
!   size, the error tolerance, and the weights, none of the parameters
!   should be altered.  The array  WT  must be updated after each step
!   to maintain relative error tests like those above.  Normally the
!   integration is continued just beyond the desired endpoint and the
!   solution interpolated there with subroutine  SINTRP .  If it is
!   impossible to integrate beyond the endpoint, the step size may be
!   reduced to hit the endpoint since the code will not take a step
!   larger than the  H  input.  Changing the direction of integration,
!   i.e., the sign of  H , requires the user set  START = .TRUE. before
!   calling  DSTEPS  again.  This is the only situation in which  START
!   should be altered.
!
!   Output from DSTEPS
!
!      Successful Step --
!
!   The subroutine returns after each successful step with  START  and
!   CRASH  set .FALSE. .  X  represents the independent variable
!   advanced one step of length  HOLD  from its value on input and  Y
!   the solution vector at the new value of  X .  All other parameters
!   represent information corresponding to the new  X  needed to
!   continue the integration.
!
!      Unsuccessful Step --
!
!   When the error tolerance is too small for the machine precision,
!   the subroutine returns without taking a step and  CRASH = .TRUE. .
!   An appropriate step size and error tolerance for continuing are
!   estimated and all other information is restored as upon input
!   before returning.  To continue with the larger tolerance, the user
!   just calls the code again.  A restart is neither required nor
!   desirable.
!
!***REFERENCES  L. F. Shampine and M. K. Gordon, Solving ordinary
!                 differential equations with ODE, STEP, and INTRP,
!                 Report SLA-73-1060, Sandia Laboratories, 1973.
!***ROUTINES CALLED  D1MACH, DHSTRT
!***REVISION HISTORY  (YYMMDD)
!   740101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DSTEPS
!
      INTEGER I, IFAIL, IM1, IP1, IPAR, IQ, J, K, KM1, KM2, KNEW,&
            KOLD, KP1, KP2, KSTEPS, L, LIMIT1, LIMIT2, NEQN, NS, NSM2,&
            NSP1, NSP2
      DOUBLE PRECISION ABSH, ALPHA, BETA, BIG,&
            EPS, ERK, ERKM1, ERKM2, ERKP1, ERR,&
            FOURU, G, GI, GSTR, H, HNEW, HOLD, P, P5EPS, PHI, PSI, R,&
            REALI, REALNS, RHO, ROUND, RPAR, SIG, TAU, TEMP1,&
            TEMP2, TEMP3, TEMP4, TEMP5, TEMP6, TWO, TWOU, U, V, W, WT,&
            X, XOLD, Y, YP
      LOGICAL START,CRASH,PHASE1,NORND
      DIMENSION Y(*),WT(*),PHI(NEQN,16),P(*),YP(*),PSI(12),&
        ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10),&
        RPAR(*),IPAR(*)
      DIMENSION TWO(13),GSTR(13)
      EXTERNAL DF
      SAVE TWO, GSTR
!
      DATA TWO(1),TWO(2),TWO(3),TWO(4),TWO(5),TWO(6),TWO(7),TWO(8),&
           TWO(9),TWO(10),TWO(11),TWO(12),TWO(13)&
           /2.0D0,4.0D0,8.0D0,16.0D0,32.0D0,64.0D0,128.0D0,256.0D0,&
            512.0D0,1024.0D0,2048.0D0,4096.0D0,8192.0D0/
      DATA GSTR(1),GSTR(2),GSTR(3),GSTR(4),GSTR(5),GSTR(6),GSTR(7),&
           GSTR(8),GSTR(9),GSTR(10),GSTR(11),GSTR(12),GSTR(13)&
           /0.5D0,0.0833D0,0.0417D0,0.0264D0,0.0188D0,0.0143D0,0.0114D0,&
            0.00936D0,0.00789D0,0.00679D0,0.00592D0,0.00524D0,0.00468D0/
!
!       ***     BEGIN BLOCK 0     ***
!   CHECK IF STEP SIZE OR ERROR TOLERANCE IS TOO SMALL FOR MACHINE
!   PRECISION.  IF FIRST STEP, INITIALIZE PHI ARRAY AND ESTIMATE A
!   STARTING STEP SIZE.
!                   ***
!
!   IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE
!
!***FIRST EXECUTABLE STATEMENT  DSTEPS
      CRASH = .TRUE.
      IF(ABS(H) >= FOURU*ABS(X)) GO TO 5
      H = SIGN(FOURU*ABS(X),H)
      RETURN
 5    P5EPS = 0.5D0*EPS
!
!   IF ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE
!
      ROUND = 0.0D0
      DO 10 L = 1,NEQN
 10     ROUND = ROUND + (Y(L)/WT(L))**2
      ROUND = TWOU*SQRT(ROUND)
      IF(P5EPS >= ROUND) GO TO 15
      EPS = 2.0D0*ROUND*(1.0D0 + FOURU)
      RETURN
 15   CRASH = .FALSE.
      G(1) = 1.0D0
      G(2) = 0.5D0
      SIG(1) = 1.0D0
      IF(.NOT.START) GO TO 99
!
!   INITIALIZE.  COMPUTE APPROPRIATE STEP SIZE FOR FIRST STEP
!
!     CALL DF(X,Y,YP,RPAR,IPAR)
!     SUM = 0.0
      DO 20 L = 1,NEQN
        PHI(L,1) = YP(L)
   20   PHI(L,2) = 0.0D0
!20     SUM = SUM + (YP(L)/WT(L))**2
!     SUM = SQRT(SUM)
!     ABSH = ABS(H)
!     IF(EPS < 16.0*SUM*H*H) ABSH = 0.25*SQRT(EPS/SUM)
!     H = SIGN(MAX(ABSH,FOURU*ABS(X)),H)
!
      U = d1mach4
      BIG = SQRT(d1mach2)
      CALL DHSTRT(DF,NEQN,X,X+H,Y,YP,WT,1,U,BIG,&
                   PHI(1,3),PHI(1,4),PHI(1,5),PHI(1,6),RPAR,IPAR,H)
!
      HOLD = 0.0D0
      K = 1
      KOLD = 0
      KPREV = 0
      START = .FALSE.
      PHASE1 = .TRUE.
      NORND = .TRUE.
      IF(P5EPS > 100.0D0*ROUND) GO TO 99
      NORND = .FALSE.
      DO 25 L = 1,NEQN
 25     PHI(L,15) = 0.0D0
 99   IFAIL = 0
!       ***     END BLOCK 0     ***
!
!       ***     BEGIN BLOCK 1     ***
!   COMPUTE COEFFICIENTS OF FORMULAS FOR THIS STEP.  AVOID COMPUTING
!   THOSE QUANTITIES NOT CHANGED WHEN STEP SIZE IS NOT CHANGED.
!                   ***
!
 100  KP1 = K+1
      KP2 = K+2
      KM1 = K-1
      KM2 = K-2
!
!   NS IS THE NUMBER OF DSTEPS TAKEN WITH SIZE H, INCLUDING THE CURRENT
!   ONE.  WHEN K<NS, NO COEFFICIENTS CHANGE
!
      IF(H /= HOLD) NS = 0
      IF (NS<=KOLD) NS = NS+1
      NSP1 = NS+1
      IF (K < NS) GO TO 199
!
!   COMPUTE THOSE COMPONENTS OF ALPHA(*),BETA(*),PSI(*),SIG(*) WHICH
!   ARE CHANGED
!
      BETA(NS) = 1.0D0
      REALNS = NS
      ALPHA(NS) = 1.0D0/REALNS
      TEMP1 = H*REALNS
      SIG(NSP1) = 1.0D0
      IF(K < NSP1) GO TO 110
      DO 105 I = NSP1,K
        IM1 = I-1
        TEMP2 = PSI(IM1)
        PSI(IM1) = TEMP1
        BETA(I) = BETA(IM1)*PSI(IM1)/TEMP2
        TEMP1 = TEMP2 + H
        ALPHA(I) = H/TEMP1
        REALI = I
 105    SIG(I+1) = REALI*ALPHA(I)*SIG(I)
 110  PSI(K) = TEMP1
!
!   COMPUTE COEFFICIENTS G(*)
!
!   INITIALIZE V(*) AND SET W(*).
!
      IF(NS > 1) GO TO 120
      DO 115 IQ = 1,K
        TEMP3 = IQ*(IQ+1)
        V(IQ) = 1.0D0/TEMP3
 115    W(IQ) = V(IQ)
      IVC = 0
      KGI = 0
      IF (K == 1) GO TO 140
      KGI = 1
      GI(1) = W(2)
      GO TO 140
!
!   IF ORDER WAS RAISED, UPDATE DIAGONAL PART OF V(*)
!
 120  IF(K <= KPREV) GO TO 130
      IF (IVC == 0) GO TO 122
      JV = KP1 - IV(IVC)
      IVC = IVC - 1
      GO TO 123
 122  JV = 1
      TEMP4 = K*KP1
      V(K) = 1.0D0/TEMP4
      W(K) = V(K)
      IF (K /= 2) GO TO 123
      KGI = 1
      GI(1) = W(2)
 123  NSM2 = NS-2
      IF(NSM2 < JV) GO TO 130
      DO 125 J = JV,NSM2
        I = K-J
        V(I) = V(I) - ALPHA(J+1)*V(I+1)
 125    W(I) = V(I)
      IF (I /= 2) GO TO 130
      KGI = NS - 1
      GI(KGI) = W(2)
!
!   UPDATE V(*) AND SET W(*)
!
 130  LIMIT1 = KP1 - NS
      TEMP5 = ALPHA(NS)
      DO 135 IQ = 1,LIMIT1
        V(IQ) = V(IQ) - TEMP5*V(IQ+1)
 135    W(IQ) = V(IQ)
      G(NSP1) = W(1)
      IF (LIMIT1 == 1) GO TO 137
      KGI = NS
      GI(KGI) = W(2)
 137  W(LIMIT1+1) = V(LIMIT1+1)
      IF (K >= KOLD) GO TO 140
      IVC = IVC + 1
      IV(IVC) = LIMIT1 + 2
!
!   COMPUTE THE G(*) IN THE WORK VECTOR W(*)
!
 140  NSP2 = NS + 2
      KPREV = K
      IF(KP1 < NSP2) GO TO 199
      DO 150 I = NSP2,KP1
        LIMIT2 = KP2 - I
        TEMP6 = ALPHA(I-1)
        DO 145 IQ = 1,LIMIT2
 145      W(IQ) = W(IQ) - TEMP6*W(IQ+1)
 150    G(I) = W(1)
 199    CONTINUE
!       ***     END BLOCK 1     ***
!
!       ***     BEGIN BLOCK 2     ***
!   PREDICT A SOLUTION P(*), EVALUATE DERIVATIVES USING PREDICTED
!   SOLUTION, ESTIMATE LOCAL ERROR AT ORDER K AND ERRORS AT ORDERS K,
!   K-1, K-2 AS IF CONSTANT STEP SIZE WERE USED.
!                   ***
!
!   INCREMENT COUNTER ON ATTEMPTED DSTEPS
!
      KSTEPS = KSTEPS + 1
!
!   CHANGE PHI TO PHI STAR
!
      IF(K < NSP1) GO TO 215
      DO 210 I = NSP1,K
        TEMP1 = BETA(I)
        DO 205 L = 1,NEQN
 205      PHI(L,I) = TEMP1*PHI(L,I)
 210    CONTINUE
!
!   PREDICT SOLUTION AND DIFFERENCES
!
 215  DO 220 L = 1,NEQN
        PHI(L,KP2) = PHI(L,KP1)
        PHI(L,KP1) = 0.0D0
 220    P(L) = 0.0D0
      DO 230 J = 1,K
        I = KP1 - J
        IP1 = I+1
        TEMP2 = G(I)
        DO 225 L = 1,NEQN
          P(L) = P(L) + TEMP2*PHI(L,I)
 225      PHI(L,I) = PHI(L,I) + PHI(L,IP1)
 230    CONTINUE
      IF(NORND) GO TO 240
      DO 235 L = 1,NEQN
        TAU = H*P(L) - PHI(L,15)
        P(L) = Y(L) + TAU
 235    PHI(L,16) = (P(L) - Y(L)) - TAU
      GO TO 250
 240  DO 245 L = 1,NEQN
 245    P(L) = Y(L) + H*P(L)
 250  XOLD = X
      X = X + H
      ABSH = ABS(H)
      CALL DF(X,P,YP,RPAR,IPAR)
!
!   ESTIMATE ERRORS AT ORDERS K,K-1,K-2
!
      ERKM2 = 0.0D0
      ERKM1 = 0.0D0
      ERK = 0.0D0
      DO 265 L = 1,NEQN
        TEMP3 = 1.0D0/WT(L)
        TEMP4 = YP(L) - PHI(L,1)
        IF(KM2)265,260,255
 255    ERKM2 = ERKM2 + ((PHI(L,KM1)+TEMP4)*TEMP3)**2
 260    ERKM1 = ERKM1 + ((PHI(L,K)+TEMP4)*TEMP3)**2
 265    ERK = ERK + (TEMP4*TEMP3)**2
      IF(KM2)280,275,270
 270  ERKM2 = ABSH*SIG(KM1)*GSTR(KM2)*SQRT(ERKM2)
 275  ERKM1 = ABSH*SIG(K)*GSTR(KM1)*SQRT(ERKM1)
 280  TEMP5 = ABSH*SQRT(ERK)
      ERR = TEMP5*(G(K)-G(KP1))
      ERK = TEMP5*SIG(KP1)*GSTR(K)
      KNEW = K
!
!   TEST IF ORDER SHOULD BE LOWERED
!
      IF(KM2)299,290,285
 285  IF(MAX(ERKM1,ERKM2) <= ERK) KNEW = KM1
      GO TO 299
 290  IF(ERKM1 <= 0.5D0*ERK) KNEW = KM1
!
!   TEST IF STEP SUCCESSFUL
!
 299  IF(ERR <= EPS) GO TO 400
!       ***     END BLOCK 2     ***
!
!       ***     BEGIN BLOCK 3     ***
!   THE STEP IS UNSUCCESSFUL.  RESTORE  X, PHI(*,*), PSI(*) .
!   IF THIRD CONSECUTIVE FAILURE, SET ORDER TO ONE.  IF STEP FAILS MORE
!   THAN THREE TIMES, CONSIDER AN OPTIMAL STEP SIZE.  DOUBLE ERROR
!   TOLERANCE AND RETURN IF ESTIMATED STEP SIZE IS TOO SMALL FOR MACHINE
!   PRECISION.
!                   ***
!
!   RESTORE X, PHI(*,*) AND PSI(*)
!
      PHASE1 = .FALSE.
      X = XOLD
      DO 310 I = 1,K
        TEMP1 = 1.0D0/BETA(I)
        IP1 = I+1
        DO 305 L = 1,NEQN
 305      PHI(L,I) = TEMP1*(PHI(L,I) - PHI(L,IP1))
 310    CONTINUE
      IF(K < 2) GO TO 320
      DO 315 I = 2,K
 315    PSI(I-1) = PSI(I) - H
!
!   ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP
!   SIZE
!
 320  IFAIL = IFAIL + 1
      TEMP2 = 0.5D0
      IF(IFAIL - 3) 335,330,325
 325  IF(P5EPS < 0.25D0*ERK) TEMP2 = SQRT(P5EPS/ERK)
 330  KNEW = 1
 335  H = TEMP2*H
      K = KNEW
      NS = 0
      IF(ABS(H) >= FOURU*ABS(X)) GO TO 340
      CRASH = .TRUE.
      H = SIGN(FOURU*ABS(X),H)
      EPS = EPS + EPS
      RETURN
 340  GO TO 100
!       ***     END BLOCK 3     ***
!
!       ***     BEGIN BLOCK 4     ***
!   THE STEP IS SUCCESSFUL.  CORRECT THE PREDICTED SOLUTION, EVALUATE
!   THE DERIVATIVES USING THE CORRECTED SOLUTION AND UPDATE THE
!   DIFFERENCES.  DETERMINE BEST ORDER AND STEP SIZE FOR NEXT STEP.
!                   ***
 400  KOLD = K
      HOLD = H
!
!   CORRECT AND EVALUATE
!
      TEMP1 = H*G(KP1)
      IF(NORND) GO TO 410
      DO 405 L = 1,NEQN
        TEMP3 = Y(L)
        RHO = TEMP1*(YP(L) - PHI(L,1)) - PHI(L,16)
        Y(L) = P(L) + RHO
        PHI(L,15) = (Y(L) - P(L)) - RHO
 405    P(L) = TEMP3
      GO TO 420
 410  DO 415 L = 1,NEQN
        TEMP3 = Y(L)
        Y(L) = P(L) + TEMP1*(YP(L) - PHI(L,1))
 415    P(L) = TEMP3
 420  CALL DF(X,Y,YP,RPAR,IPAR)
!
!   UPDATE DIFFERENCES FOR NEXT STEP
!
      DO 425 L = 1,NEQN
        PHI(L,KP1) = YP(L) - PHI(L,1)
 425    PHI(L,KP2) = PHI(L,KP1) - PHI(L,KP2)
      DO 435 I = 1,K
        DO 430 L = 1,NEQN
 430      PHI(L,I) = PHI(L,I) + PHI(L,KP1)
 435    CONTINUE
!
!   ESTIMATE ERROR AT ORDER K+1 UNLESS:
!     IN FIRST PHASE WHEN ALWAYS RAISE ORDER,
!     ALREADY DECIDED TO LOWER ORDER,
!     STEP SIZE NOT CONSTANT SO ESTIMATE UNRELIABLE
!
      ERKP1 = 0.0D0
      IF(KNEW == KM1  .OR.  K == 12) PHASE1 = .FALSE.
      IF(PHASE1) GO TO 450
      IF(KNEW == KM1) GO TO 455
      IF(KP1 > NS) GO TO 460
      DO 440 L = 1,NEQN
 440    ERKP1 = ERKP1 + (PHI(L,KP2)/WT(L))**2
      ERKP1 = ABSH*GSTR(KP1)*SQRT(ERKP1)
!
!   USING ESTIMATED ERROR AT ORDER K+1, DETERMINE APPROPRIATE ORDER
!   FOR NEXT STEP
!
      IF(K > 1) GO TO 445
      IF(ERKP1 >= 0.5D0*ERK) GO TO 460
      GO TO 450
 445  IF(ERKM1 <= MIN(ERK,ERKP1)) GO TO 455
      IF(ERKP1 >= ERK  .OR.  K == 12) GO TO 460
!
!   HERE ERKP1 < ERK < MAX(ERKM1,ERKM2) ELSE ORDER WOULD HAVE
!   BEEN LOWERED IN BLOCK 2.  THUS ORDER IS TO BE RAISED
!
!   RAISE ORDER
!
 450  K = KP1
      ERK = ERKP1
      GO TO 460
!
!   LOWER ORDER
!
 455  K = KM1
      ERK = ERKM1
!
!   WITH NEW ORDER DETERMINE APPROPRIATE STEP SIZE FOR NEXT STEP
!
 460  HNEW = H + H
      IF(PHASE1) GO TO 465
      IF(P5EPS >= ERK*TWO(K+1)) GO TO 465
      HNEW = H
      IF(P5EPS >= ERK) GO TO 465
      TEMP2 = K+1
      R = (P5EPS/ERK)**(1.0D0/TEMP2)
      HNEW = ABSH*MAX(0.5D0,MIN(0.9D0,R))
      HNEW = SIGN(MAX(HNEW,FOURU*ABS(X)),H)
 465  H = HNEW
      RETURN
!       ***     END BLOCK 4     ***
      END SUBROUTINE DSTEPS

    subroutine xermsg (librar, subrou, messg, nerr, level)

    implicit none

    character(len=*),intent(in) :: librar, subrou, messg
    integer,intent(in) :: nerr, level

    write(*,'(A)') trim(messg)      

    end subroutine xermsg

    end module ddeabm_module
