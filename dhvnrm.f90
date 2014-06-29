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
      END
