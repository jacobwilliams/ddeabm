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
            DFDUB, DFDXB, DHVNRM,&
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
      END
