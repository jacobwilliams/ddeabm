!DECK FDUMP
      SUBROUTINE FDUMP
!***BEGIN PROLOGUE  FDUMP
!***PURPOSE  Symbolic dump (should be locally written).
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (FDUMP-A)
!***KEYWORDS  ERROR, XERMSG
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!        ***Note*** Machine Dependent Routine
!        FDUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!        XGETUA routine.  See XSETUA and XGETUA for details.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  FDUMP
!***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
